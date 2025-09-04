#### CE script sert juste à récuperer les projections d'incidence d'esca, 
# et de calculer la moyenne annuelle selon l'ensemble des modèles climatiques
# et de calculer la moyenne par période et par modèles climatiques

library(dplyr)

# load("data/Projections/Res/esca_historical_ACCESS-CM2_1995-2014.RData") 
# Récuperer l'ensemble des projections par ordre alphabétique donc par scénario puis par modèle puis par dates
cmip6_files <- sort(list.files("data/Projections/Res/"))


hist <- cmip6_files[1:19]

hist_complet <- data.frame()

for(data in hist){
  load(paste0("data/Projections/Res/", data))
  model <- strsplit(data, "_")[[1]][3]
  hist_complet <- rbind(hist_complet, cbind(projections, model))
  rm(projections); gc()
}

hist_moy <- hist_complet[,1:7] %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, an) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(hist_moy, file = paste0("data/Projections/Res_moy/esca_moy_historical.RData"))

hist_moy <- hist_complet %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, model) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(hist_moy, file = paste0("data/Projections/Res_period/esca_moy_historical.RData"))

scen_245_2046_2065_complet <- data.frame()
scen_245_2081_2100_complet <- data.frame()
scen_585_2046_2065_complet <- data.frame()
scen_585_2081_2100_complet <- data.frame()

for(data in cmip6_files[20:93]){
  scen <- strsplit(data, "_")[[1]][2]
  period <- strsplit(data, "_")[[1]][4]  
  model <- strsplit(data, "_")[[1]][3]
  load(paste0("data/Projections/Res/", data))
  if(scen == "ssp245"){
    if(period == "2046-2065.RData") scen_245_2046_2065_complet <- rbind(scen_245_2046_2065_complet, cbind(projections, model))
    else scen_245_2081_2100_complet <- rbind(scen_245_2081_2100_complet, cbind(projections, model))
  }
  else{
    if(period == "2046-2065.RData") scen_585_2046_2065_complet <- rbind(scen_585_2046_2065_complet, cbind(projections, model))
    else scen_585_2081_2100_complet <- rbind(scen_585_2081_2100_complet, cbind(projections, model))
  }
  rm(projections); gc()
}

scen_245_2046_2065_moy <- scen_245_2046_2065_complet[,1:7] %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, an) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_245_2046_2065_moy, file = paste0("data/Projections/Res_moy/esca_moy_245_2046_2065.RData"))

scen_245_2046_2065_moy <- scen_245_2046_2065_complet %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, model) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_245_2046_2065_moy, file = paste0("data/Projections/Res_period/esca_moy_245_2046_2065.RData"))


scen_245_2081_2100_moy <- scen_245_2081_2100_complet[,1:7] %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, an) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_245_2081_2100_moy, file = paste0("data/Projections/Res_moy/esca_moy_245_2081_2100.RData"))

scen_245_2081_2100_moy <- scen_245_2081_2100_complet %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, model) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_245_2081_2100_moy, file = paste0("data/Projections/Res_period/esca_moy_245_2081_2100.RData"))


scen_585_2046_2065_moy <- scen_585_2046_2065_complet[, 1:7] %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, an) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_585_2046_2065_moy, file = paste0("data/Projections/Res_moy/esca_moy_585_2046_2065.RData"))

scen_585_2046_2065_moy <- scen_585_2046_2065_complet %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, model) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_585_2046_2065_moy, file = paste0("data/Projections/Res_period/esca_moy_585_2046_2065.RData"))


scen_585_2081_2100_moy <- scen_585_2081_2100_complet[,1:7] %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, an) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_585_2081_2100_moy, file = paste0("data/Projections/Res_moy/esca_moy_585_2081_2100.RData"))

scen_585_2081_2100_moy <- scen_585_2081_2100_complet %>% 
  group_by(age_parcelle_estime, esca_categoriel, maille, cepage, model) %>%
  summarise(moy_lm = mean(pred_lm), sd_lm = sd(pred_lm), moy_gp = mean(pred_gp), sd_gp = sd(pred_gp))

save(scen_585_2081_2100_moy, file = paste0("data/Projections/Res_period/esca_moy_585_2081_2100.RData"))


##############################################################################################

cmip6_files <- sort(list.files("data/Projections/Climat/"))

hist <- cmip6_files[1:19]

hist_complet <- data.frame()

for(data in hist){
  load(paste0("data/Projections/Climat/", data))
  model <- strsplit(data, "_")[[1]][3]
  hist_complet <- rbind(hist_complet, cbind(projections_climat, model))
  rm(projections_climat); gc()
}

hist_moy <- hist_complet %>% 
  group_by(model, maille, cepage) %>%
  summarise(debourrement = mean(debourrement), floraison = mean(floraison), tm = mean(tm), rr = mean(rr),
            ftsw = mean(ftsw))

save(hist_moy, file = paste0("data/Projections/Climat_moy/esca_moy_historical.RData"))

scen_245_2046_2065_complet <- data.frame()
scen_245_2081_2100_complet <- data.frame()
scen_585_2046_2065_complet <- data.frame()
scen_585_2081_2100_complet <- data.frame()

for(data in cmip6_files[20:91]){
  scen <- strsplit(data, "_")[[1]][2]
  period <- strsplit(data, "_")[[1]][4]  
  model <- strsplit(data, "_")[[1]][3]
  load(paste0("data/Projections/Climat/", data))
  if(scen == "ssp245"){
    if(period == "2046-2065.RData") scen_245_2046_2065_complet <- rbind(scen_245_2046_2065_complet, cbind(projections_climat, model))
    else scen_245_2081_2100_complet <- rbind(scen_245_2081_2100_complet, cbind(projections_climat, model))
  }
  else{
    if(period == "2046-2065.RData") scen_585_2046_2065_complet <- rbind(scen_585_2046_2065_complet, cbind(projections_climat, model))
    else scen_585_2081_2100_complet <- rbind(scen_585_2081_2100_complet, cbind(projections_climat, model))
  }
  rm(projections_climat); gc()
}

scen_245_2046_2065_moy <- scen_245_2046_2065_complet %>% 
  group_by(model, maille, cepage) %>%
  summarise(debourrement = mean(debourrement), floraison = mean(floraison), tm = mean(tm), rr = mean(rr),
            ftsw = mean(ftsw))

save(scen_245_2046_2065_moy, file = paste0("data/Projections/Climat_moy/esca_moy_245_2046_2065.RData"))

scen_245_2081_2100_moy <- scen_245_2081_2100_complet %>% 
  group_by(model, maille, cepage) %>%
  summarise(debourrement = mean(debourrement), floraison = mean(floraison), tm = mean(tm), rr = mean(rr),
            ftsw = mean(ftsw))

save(scen_245_2081_2100_moy, file = paste0("data/Projections/Climat_moy/esca_moy_245_2081_2100.RData"))

scen_585_2046_2065_moy <- scen_585_2046_2065_complet %>% 
  group_by(model, maille, cepage) %>%
  summarise(debourrement = mean(debourrement), floraison = mean(floraison), tm = mean(tm), rr = mean(rr),
            ftsw = mean(ftsw))

save(scen_585_2046_2065_moy, file = paste0("data/Projections/Climat_moy/esca_moy_585_2046_2065.RData"))

scen_585_2081_2100_moy <- scen_585_2081_2100_complet %>% 
  group_by(model, maille, cepage) %>%
  summarise(debourrement = mean(debourrement), floraison = mean(floraison), tm = mean(tm), rr = mean(rr),
            ftsw = mean(ftsw))

save(scen_585_2081_2100_moy, file = paste0("data/Projections/Climat_moy/esca_moy_585_2081_2100.RData"))

