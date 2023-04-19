library(tidyverse)
library(progress)
library(pbapply)
library(pbmcapply)
library(smoother)
library(patchwork)
source("R/DynamicFrame.R")

#' Run a single instance of the individual-based simulation
#' 
#' Only arguments with obvious default values have those set. Can override function parameters by passing in a named list to `params`.
run_simulation = function(i, params,
                          dt = 1,
                          tau = 1,
                          nu = 0,
                          iota_max = 0,
                          forcing=function(x) {1},
                          progress=F,
                          n_trace = 1000) {
  # Set param elements to the current environment
  for (name in names(params)) {
    assign(name, params[[name]])
  }
  
  # Prepare trace/logging variables
  n_trace = min(floor(max_t/dt), n_trace)
  max_t = max_t%/%dt*dt
  t_trace = seq(dt, max_t, length.out=n_trace)
  trace = tibble(
    t = NA_real_,
    blood = NA_integer_,
    liver = NA_integer_,
    foci = NA_integer_,
    .rows = n_trace
  )
  i_trace = 1
  
  # Initialise population
  pop = tibble(
    ID = 1:n,
    blood_stage = runif(n) < incidence,
    liver_stage = blood_stage,
    t_infect = -Inf,
    t_treat = Inf
  )
  
  if (progress) pb = progress_bar$new(total = max_t/dt + 1, format = "[:bar] ETA::eta")
  for (t in seq(dt, max_t, dt)) {
    force_of_infection = forcing(t)*lambda * sum(pop$blood_stage) / n + delta
    
    # Infect people
    to_infect = runif(n) < force_of_infection*dt*(1-secondary_immunity*(pop$t_infect>=0)) & !pop$blood_stage # don't allow concurrent infections
    to_relapse = pop$liver_stage & (runif(n) < f*dt) & !pop$blood_stage
    to_schedule_detect = (to_infect | to_relapse) & (runif(n) < rho)
    
    # Detect people passively
    to_detect = pop$t_treat <= t
    
    # Conduct RACD when treating blood
    case_detections = min(sum(to_detect), iota_max * dt)
    racd_weights = 1 + (tau-1) * (pop$blood_stage | pop$liver_stage)
    # Increasing case detections reduces marginal racd proportion because of overlap
    racd_proportion = 1-exp(-case_detections * nu / n) # 40% faster than 1-(1-nu)^c_d. This includes a deterministic reduction factor
    racd_n = n * racd_proportion # racd_n = rbinom(1, n, racd_proportion)
    to_racd = rep(F, n)
    to_racd[sample.int(n, racd_n, prob=racd_weights)] = T
    blood_detectable = pop$blood_stage & runif(n) < racd_blood_sens
    liver_detectable = pop$liver_stage & runif(n) < racd_liver_sens
    # Add to to_detect
    to_detect[to_racd & (blood_detectable | liver_detectable)] = T
    
    # Process treatment
    to_treat_blood = to_detect & runif(n) < alpha
    to_treat_liver = to_treat_blood & (runif(n) < beta)
    
    # Natural clearance
    to_clear_blood = pop$blood_stage & runif(n) < r * dt
    to_clear_liver = pop$liver_stage & runif(n) < gamma * dt
    
    # Apply changes
    pop$t_infect[to_infect | to_relapse] = t
    pop$blood_stage[to_infect | to_relapse] = T
    pop$blood_stage[to_treat_blood | to_clear_blood] = F
    pop$liver_stage[to_infect] = T
    pop$liver_stage[to_treat_liver | to_clear_liver] = F
    pop$t_treat[to_schedule_detect] = t + 14 # schedule treatment effect
    pop$t_treat[pop$t_treat <= t] = Inf # reset treatment schedule if done
    
    # Log
    if ((t+1e-10) >= t_trace[i_trace]) {
      trace$t[i_trace] = t
      trace$blood[i_trace] = sum(pop$blood_stage) / n
      trace$liver[i_trace] = sum(pop$liver_stage) / n
      trace$foci[i_trace] = case_detections / dt
      i_trace = i_trace + 1
    }
    if (progress) pb$tick()
  }
  trace
}

run_repeats = function(params, n, ..., mc.cores=8) {
  runs = pbmclapply(seq_len(n), run_simulation, params, ..., mc.cores=mc.cores)
  run_summary = bind_rows(runs) %>%
    pivot_longer(cols = -t) %>%
    group_by(t, name) %>%
    summarise(id = params$id,
              LQ = quantile(value, 0.95, na.rm=T),
              UQ = quantile(value, 0.05, na.rm=T),
              median = median(value, na.rm=T),
              mean = mean(value, na.rm=T),
              .groups = "drop")
}


params = list(id = "Default",
              incidence = 0.01,
              alpha = 0.5, 
              beta = 0.8,
              rho = 0.9,
              lambda = 0.02,
              delta = 5e-06,
              f = 1/30, # relapse freq
              r = 1/60,
              gamma = 1/223,
              tau = 2,
              nu = 50, # people per focus
              iota_max = 1, # focii per day
              racd_blood_sens = 0.2,
              racd_liver_sens = 0,
              secondary_immunity = 0.5, # scales the force of infection
              dt = 5
)

# Set up scenarios
# Scenarios for RCD
scenario.A = params
scenario.A$id = "Untargeted case detection"
scenario.A$tau = 1

scenario.B = params
scenario.B$id = "RACD, Blood stage only"

scenario.C = params
scenario.C$id = "80% serological sensitivity"
scenario.C$racd_liver_sens = 0.8

# Scenarios for checking timestep invariance
scenario.D = params
scenario.D$id = "dt=2"
scenario.D$dt = 2

scenario.E = params
scenario.E$id = "dt=0.25"
scenario.E$dt = 0.25

# Run simulations
n = 8000
max_t = 6*365.25
n_repeats = 32
scenarios = list(scenario.A, scenario.B, scenario.C)
# scenarios = list(scenario.D, scenario.E)
summaries = lapply(scenarios, run_repeats, n_repeats, forcing=function(t) {
  0.5 + 0.5*cos(2*pi*t/365.25)
})

run_summary = bind_rows(summaries) %>%
  mutate(id = fct_inorder(id)) %>%
  group_by(name, id) %>%
  mutate_at(vars(mean, median, LQ, UQ), smth.gaussian, window=10, tails=T) %>%
  ungroup() %>%
  mutate(name = fct_recode(name,
                           "Blood stage prevalence" = "blood",
                           "Liver stage prevalence" = "liver",
                           "Daily reactive case detection foci" = "foci"))

p_shared_elements = list(
  geom_ribbon(alpha = 0.3),
  geom_line(aes(y = mean, color = id)),
  scale_x_continuous(expand = c(0, 0), breaks=seq(0, floor(max(run_summary$t)))),
  labs(x = "Year", y = NULL, color = "Intervention", fill = "Intervention"),
  facet_wrap(vars(name), ncol = 1),
  theme(plot.background = element_rect(fill = "transparent", color = NA))
)

(run_summary %>%
    filter(name %>% str_detect(c("Blood|Liver"))) %>%
    ggplot(aes(x = t/365.25, ymin = LQ, ymax = UQ, fill = id)) +
    scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) +
  p_shared_elements) /
  (run_summary %>%
     filter(name %>% str_detect("reactive")) %>%
     ggplot(aes(x = t/365.25, ymin = LQ, ymax = UQ, fill = id)) +
     p_shared_elements +
     scale_y_continuous(limits = c(0, NA))) +
  plot_layout(guides = "collect", heights=c(2,1))


# ggsave("~/Downloads/plot.pdf", width=8, height=3, bg="transparent")
# 
# # Render video
# y_max = max(run_summary$UQ) *1.25
# temp_dir = tempdir()
# file = "~/Downloads/animation.mp4"
# tryCatch({
#   skip_factor = 100 # only render every nth time
#   frames = seq(nrow(run_summary) %% skip_factor, nrow(run_summary), skip_factor)
#   frame_files = file.path(temp_dir, paste0("frame_", frames, ".png"))
#   pb = progress::progress_bar$new(total = length(frame_files), format = "Rendering :current of :total frames [:bar] ETA::eta")
#   for (i in seq_along(frames)) {
#     pb$tick()
#     run_summary[1:frames[i],] %>%
#       ggplot(aes(x = t, ymin = LQ, ymax = UQ, fill = id)) +
#       geom_ribbon(alpha = 0.3) +
#       geom_line(aes(y = median, color = id)) +
#       geom_vline(xintercept = run_summary$t[frames[i]], color="grey") +
#       scale_x_continuous(expand = c(0, 0), limits=c(min(run_summary$t), max(run_summary$t))) +
#       scale_y_continuous(limits = c(0, y_max*1.25), labels = scales::label_percent()) +
#       scale_color_discrete(drop = FALSE) +
#       scale_fill_discrete(drop = FALSE) +
#       facet_wrap(vars(name), ncol = 1) +
#       labs(x = "Days", y = NULL, color = "Intervention", fill = "Intervention")
#     ggsave(frame_files[i], width=8, height=3, dpi=300)
#   }
# 
#   message("Encoding to ", file)
#   av::av_encode_video(frame_files, file, framerate = 30, verbose = F)
# },
# finally = {
#   message("Removing temporary frames from ", temp_dir)
#   # unlink(frame_files)
# })
