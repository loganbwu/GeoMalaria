## code to prepare boundary datasets goes here
crs = 24047 # https://epsg.io/24047

library(sf)
adm0 = st_read("https://github.com/wmgeolab/geoBoundaries/raw/df46a8320bb2f703640664510273bec2be88b180/releaseData/gbOpen/THA/ADM0/geoBoundaries-THA-ADM0_simplified.geojson")
adm1 = st_read("https://github.com/wmgeolab/geoBoundaries/raw/eb8c046c824d34547fbe9a50b7308b443e76c106/releaseData/gbOpen/THA/ADM1/geoBoundaries-THA-ADM1_simplified.geojson")
adm2 = st_read("https://github.com/wmgeolab/geoBoundaries/raw/922e33de54786eed936fc13b7b3fd28c50189d5b/releaseData/gbOpen/THA/ADM2/geoBoundaries-THA-ADM2_simplified.geojson")

thailand = adm0 %>% st_transform(crs)
tak = adm1[adm1$shapeName == "Tak Province",] %>% st_transform(crs)
thasongyang = adm2[adm2$shapeName == "Tha Song Yang",] %>% st_transform(crs)

usethis::use_data(thailand, tak, thasongyang, internal = TRUE, overwrite = TRUE)
