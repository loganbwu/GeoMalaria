library(tidyverse)
library(progress)
library(pbapply)
source("R/DynamicFrame.R")



params = list(id = "TSY",
              incidence = 0.05,
              prop_import = 0.35,
              alpha = 0.9, 
              beta = 0.8,
              rho = 0.9,
              h = 2.73972602739726e-05, 
              lambda = 0.04,
              delta = 1.06564364876385e-05,
              f = 1/30,
              r = 1/60,
              gamma = 1/223,
              tau = 1.01,
              nu = 5,
              iota_max = 10,
              eta = 1,
              racd_blood_sens = 0.1,
              racd_liver_sens = 0.8
)
n = 1000
dt = 1
max_t = 1000

simulation = function(i) {
  pop = tibble(
    ID = 1:n,
    blood_stage = runif(n) < params$incidence,
    liver_stage = blood_stage,
    t_treat = ifelse(runif(n) < 0.2, runif(n, 0, 14), Inf)
  )
  
  trace = DynamicFrame$new(
    t = 0,
    blood = sum(pop$blood_stage) / n,
    liver = sum(pop$liver_stage) / n
  )
  
  # pb = progress_bar$new(total = round(max_t / dt) + 1, format = "[:bar] ETA::eta")
  for (t in seq(0, max_t, dt)) {
    I = sum(pop$blood_stage) / n
    FoI = params$lambda * I + params$delta
    
    # Infect people
    to_infect = runif(n) < FoI*dt & !pop$blood_stage # don't allow concurrent infections
    to_relapse = pop$liver_stage & (runif(n) < params$f*dt) & !pop$blood_stage
    to_schedule_detect = (to_infect | to_relapse) & (runif(n) < params$rho)
    
    # Detect people passively
    to_detect = pop$t_treat <= t
    
    # Conduct RACD when treating blood
    iota = min(sum(to_detect), params$iota_max)
    no_tests = iota * params$nu
    racd_test = rep(F, n)
    racd_test[sample.int(n, size=no_tests, prob = 1 + (params$tau-1) * (pop$blood_stage | pop$liver_stage))] = T
    blood_detectable = pop$blood_stage & runif(n) < params$racd_blood_sens
    liver_detectable = pop$liver_stage & runif(n) < params$racd_liver_sens
    # Add to to_detect
    to_detect[racd_test & (blood_detectable | liver_detectable)] = T
    
    # Process treatment
    to_treat_blood = to_detect & runif(n) < params$alpha
    to_treat_liver = to_treat_blood & (runif(n) < params$beta)
    
    # Natural clearance
    to_clear_blood = pop$blood_stage & runif(n) < params$r * dt
    to_clear_liver = pop$liver_stage & runif(n) < params$gamma * dt
    
    # Apply changes
    pop$blood_stage[to_infect | to_relapse] = T
    pop$blood_stage[to_treat_blood | to_clear_blood] = F
    pop$liver_stage[to_infect] = T
    pop$liver_stage[to_treat_liver | to_clear_liver] = F
    pop$t_treat[to_schedule_detect] = t + 14 # if infected, schedule treatment of a proportion
    pop$t_treat[pop$t_treat <= t] = Inf # reset treatment if done
    
    # Log
    if (t %% 5 == 0) {
      trace$append(t = t,
                   blood = sum(pop$blood_stage) / n,
                   liver = sum(pop$liver_stage) / n)
    }
    # pb$tick()
  }
  trace$df
}

runs = pblapply(1:10, simulation)
run_summary = bind_rows(runs) %>%
  pivot_longer(cols = c("blood", "liver")) %>%
  group_by(t, name) %>%
  summarise(LQ = quantile(value, 0.75),
            UQ = quantile(value, 0.25),
            median = median(value),
            .groups = "drop")

ggplot(run_summary, aes(x = t, ymin = LQ, ymax = UQ, fill = name)) +
  geom_ribbon(alpha = 0.5) +
  geom_line(aes(y = median, color = name)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(n.breaks = 11, limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Time", y = "Proportion", color = "Prevalence", fill = "Prevalence")
