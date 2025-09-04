# 
# Script Sébastien ZITO ; sebastien.zito@inrae.fr
#        Gabriel MACE ; gabriel.mace@inrae.fr
#        Chloe DELMAS ; chloe.delmas@inrae.fr
# 
#  avril et mai 2025
# 
# Ce script à pour but de simuler la phénologoie et le bilan hydrique correspondant à nos données 
# pour chaque couple cépage maille.
# On simule pour toutes les années, même si on n'a pas d'observations d'ESCA.
# On récupère aussi les tendances climatiques sur les différentes périodes phénologiques.


library(dplyr)
# library(progressr) # Pour voir la progression des simulations
library(parallel)
source("script/Bilan_hydrique/AE.clim.ind_V0.8.3_Gabriel.r") # Fonctions bilan hydrique
load("data/climat/climat_maille_cepage_jour_complet.RData") # Données climatiques par maille / cépage
description_parcelles <- read.csv("data/Reg_viti_water_balance_param.csv") # Pratiques agricoles régionales
description_parcelles[description_parcelles$cepage == "cabernet-sauvignon",]$cepage <- "cabernet sauvignon"

info_mailles <- read.csv("data/geographie/correspondance_maille_commune.csv", sep=",")
info_mailles <- info_mailles %>% 
  select(X_maille, Y_maille, ID_MAILLE, RU, altitude, region_viticole) %>% unique()


complet <- complet %>% filter(DATE > 20020000 & DATE < 20240000) # Que les années à partir de 2002, sans 2024

complet <- complet %>%
  mutate(DATE = as.Date(as.character(DATE), format = "%Y%m%d"), YEAR = format(DATE, "%Y"),
         pr = PRENEI_Q + PRELIQ_Q)

# Simulation bilan hydrique : ----
bilan_hydrique <- function(cep, id_maille, region, ru=150, altitude = 45, daily = FALSE){
  climat <- complet %>% # Récupération du climat de la maille
            filter(ID_MAILLE == id_maille, cepage == cep) %>% 
            select(tasmin = TINF_H_Q, tasmax = TSUP_H_Q, pr, Time = DATE, hu = HU_Q, swi = SWI_Q) # Enlver swi pour les projections
  to_use <- description_parcelles[description_parcelles$region == region,] # Pratiques agricoles régionales
  # Simulations BU :
  res <- AgroEcoclim_index(df = climat, lon = 45, lat = 45, alt = altitude, cepage = cep,
                           mat_sugar_content = to_use$mat_sugar_content, canopy.height = to_use$canopy_height,
                           canopy.width = to_use$canopy_width, distance.between.rows = to_use$rows_distance,
                           row.porosity = to_use$row_porosity, row.azimut = to_use$row_azimut, ru = ru,
                           daily.wb.outputs = daily)
  return(res)
}

param_periodes <- c("an", "longueur_periode", "tn","tx","tm", "rr", "VPD", "hu", "rain.days", "swi",
                    "ftsw","isv", "bh0", "bhv", "bh", "et0", "tv", "auc_isv",
                    "sum.days.isv.faible", "sum.days.isv.fai_mod", "sum.days.isv.mod_sev", "sum.days.isv.sev",
                    "sum.frost.days.0", "sum.heat.days.25", "sum.heat.days.30", "sum.heat.days.35",
                    "isv.faible.seq.5", "isv.faible.seq.10","isv.faible.seq.15",
                    "isv.fai_mod.seq.5", "isv.fai_mod.seq.10","isv.fai_mod.seq.15",
                    "isv.mod_sev.seq.5", "isv.mod_sev.seq.10","isv.mod_sev.seq.15",
                    "isv.sev.seq.5", "isv.sev.seq.10","isv.sev.seq.15")# ru swhc, "diff.ftsw.days"


# handlers("progress") # Pour affiché la barre de progression : 
pre <- Sys.time()
# with_progress({
# progress <- progressr::progressor(steps = length(unique(complet$ID_MAILLE)))

#indicateurs_climatiques <- c()
#for(maille in unique(complet$ID_MAILLE)){ # Pour chaque maille
cl <- makeCluster(8)
clusterExport(cl, varlist = ls(globalenv()))
clusterEvalQ(cl, {c(library(dplyr), library(tis))})

sim.result.list <- parLapply(cl, unique(complet$ID_MAILLE), function(maille){
  res_list <- list()
  for(cepage in unique(complet[complet$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    info_maille <- info_mailles %>% filter(ID_MAILLE == maille)
    res <- bilan_hydrique(cepage, maille ,info_maille$region_viticole, info_maille$RU, info_maille$altitude) 
    bh <- data.frame(maille = maille, cepage = cepage, an = res$ph$an, debourrement = res$ph$C, 
                     floraison = res$ph$I)
    bh <- bh %>% filter(!an == 2002) # Enlever la première année car pas de débourrement bien simulé
    bh <- bh %>% left_join(res[["st.an.pheno"]][,param_periodes], by = "an")
    bh <- bh %>% left_join(res[["st.an"]][,param_periodes], by = "an", suffix = c("",".an"))
    bh <- bh %>% left_join(res[["st.dormance"]][,param_periodes], by = "an", suffix = c("",".dormance"))
    bh <- bh %>% left_join(res[["st.deb.to.flo"]][,param_periodes], by = "an", suffix = c("",".deb_flo"))
    bh <- bh %>% left_join(res[["st.symptomes"]][,param_periodes], by = "an", suffix = c("",".symptomes"))
    bh <- bh %>% left_join(res[["st.deb.to.end"]][,param_periodes], by = "an", suffix = c("",".deb_end"))
    # indicateurs_climatiques <- rbind(indicateurs_climatiques,bh)
    res_list[[length(res_list) + 1]] <- bh
  }
  do.call(rbind, res_list)
  #progress()
#}
})
stopCluster(cl)
indicateurs_climatiques <- do.call(rbind, sim.result.list)
post <- Sys.time()
print(post-pre) # 16.27617 secs, 10 premières mailles : 1.746367 mins, totale = 25.33977 mins
# parallélisé : 10 premières mailles : 41.58513 secs, totale =~10 mins

save(indicateurs_climatiques, file = "data/climat/simulations_pheno_bh_obs.RData")

# Lier aux observations ESCA : 
esca <- read.csv("data/ESCA/esca_filtre.csv") # Données ESCA
esca <- esca %>% select(annee, code_commune_INSEE, cepage, age_parcelle_estime, pourcentage_esca)
maille_commune <- read.csv("data/geographie/correspondance_maille_commune.csv", sep=",")
maille_commune <- maille_commune %>% select(ID_MAILLE, INSEE_COM, region_viticole, RU) %>% unique()
esca <- esca %>% left_join(maille_commune, by = c("code_commune_INSEE" = "INSEE_COM"))

observations <- esca %>% left_join(indicateurs_climatiques, by = c("ID_MAILLE" = "maille", "cepage" = "cepage",
                                                           "annee" = "an"), relationship = "many-to-one")
save(observations, file = "data/modelisation/observations.RData")



### Récupération des valeurs quotidiennes de certaines variables (pour 2019 (sèche) et 2021 (humide)) :

complet <- complet %>% filter(YEAR %in% (2018:2021)) # Que les années de 2018 à 2021


cl <- makeCluster(8)
clusterExport(cl, varlist = ls(globalenv()))
clusterEvalQ(cl, {c(library(dplyr), library(tis))})

sim.result.list <- parLapply(cl, unique(complet$ID_MAILLE), function(maille){
  res_list <- list()
  for(cepage in unique(complet[complet$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    info_maille <- info_mailles %>% filter(ID_MAILLE == maille)
    res <- bilan_hydrique(cepage, maille ,info_maille$region_viticole, info_maille$RU, info_maille$altitude,
                          daily = TRUE) 
    
    bh <- data.frame(maille = maille, cepage = cepage, an = res$ph$an, debourrement = res$ph$C, 
                     floraison = res$ph$I)
    bh <- bh %>% filter(an %in% c(2019, 2021)) # Garder que 2019 et 2021.
    bh <- bh %>% left_join(res[["daily.wb"]], by = "an")
    res_list[[length(res_list) + 1]] <- bh
  }
  do.call(rbind, res_list)
})
stopCluster(cl)
obs_an <- do.call(rbind, sim.result.list)


# Lier aux observations ESCA : 
esca <- esca %>% filter(annee %in% c(2019,2021))
observations_2019_2021 <- esca %>% left_join(obs_an, by = c("ID_MAILLE" = "maille", "cepage" = "cepage",
                                                            "annee" = "an"), relationship = "many-to-many")
save(observations_2019_2021, file = "data/modelisation/observations_2019_2021.RData")



#############################################################################################################
# Garder l'identifiant de parcelle :
load("data/climat/simulations_pheno_bh_obs.RData")
esca <- read.csv("data/ESCA/esca_filtre.csv") # Données ESCA
esca <- esca %>% select(annee, code_commune_INSEE, cepage, age_parcelle_estime, pourcentage_esca, 
                        identifiant_parcelle_analyse)
maille_commune <- read.csv("data/geographie/correspondance_maille_commune.csv", sep=",")
maille_commune <- maille_commune %>% select(ID_MAILLE, INSEE_COM, region_viticole, RU) %>% unique()
esca <- esca %>% left_join(maille_commune, by = c("code_commune_INSEE" = "INSEE_COM"))

observations <- esca %>% left_join(indicateurs_climatiques, by = c("ID_MAILLE" = "maille", "cepage" = "cepage",
                                                                   "annee" = "an"), relationship = "many-to-one")
save(observations, file = "data/modelisation/observations_parcelles.RData")
################################################################################################################


#############################################################################################################
# Garder l'identifiant de parcelle et le climat n-1 :
load("data/climat/simulations_pheno_bh_obs.RData")
esca <- read.csv("data/ESCA/esca_filtre.csv") # Données ESCA
esca <- esca %>% select(annee, code_commune_INSEE, cepage, age_parcelle_estime, pourcentage_esca, 
                        identifiant_parcelle_analyse)
maille_commune <- read.csv("data/geographie/correspondance_maille_commune.csv", sep=",")
maille_commune <- maille_commune %>% select(ID_MAILLE, INSEE_COM, region_viticole, RU) %>% unique()
esca <- esca %>% left_join(maille_commune, by = c("code_commune_INSEE" = "INSEE_COM"))

observations <- esca %>% left_join(indicateurs_climatiques, by = c("ID_MAILLE" = "maille", "cepage" = "cepage",
                                                                   "annee" = "an"), relationship = "many-to-one")
save(observations, file = "data/modelisation/observations_parcelles.RData")
################################################################################################################

