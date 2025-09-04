
# Fonctions nécessaires pour la simulation des incidences d'esca

# Simulation bilan hydrique : ----
bilan_hydrique <- function(meteo_maille, description_parc, cep, id_maille, region, ru=150, altitude = 45, daily = FALSE){
  climat <- meteo_maille %>% # Récupération du climat de la maille
    filter(ID_MAILLE == id_maille, cepage == cep) %>% 
    select(tasmin = tm, tasmax = tx, pr, Time = date, hu = hu) 
  to_use <- description_parc[description_parc$region == region,] # Pratiques agricoles régionales

  res <- AgroEcoclim_index(df = climat, lon = 45, lat = 45, alt = altitude, cepage = cep,
                           mat_sugar_content = to_use$mat_sugar_content, canopy.height = to_use$canopy_height,
                           canopy.width = to_use$canopy_width, distance.between.rows = to_use$rows_distance,
                           row.porosity = to_use$row_porosity, row.azimut = to_use$row_azimut, ru = ru,
                           daily.wb.outputs = daily)
  return(res)
}

# Les noms des variables à récupérer pour chaque période :
param_periodes <- c("an", "longueur_periode", "tn","tx","tm", "rr", "VPD", "hu", "rain.days", 
                    "ftsw","isv", "bh0", "bhv", "bh", "et0", "tv", "auc_isv",
                    "sum.days.isv.faible", "sum.days.isv.fai_mod", "sum.days.isv.mod_sev", "sum.days.isv.sev",
                    "sum.frost.days.0", "sum.heat.days.25", "sum.heat.days.30", "sum.heat.days.35",
                    "isv.faible.seq.5", "isv.faible.seq.10","isv.faible.seq.15",
                    "isv.fai_mod.seq.5", "isv.fai_mod.seq.10","isv.fai_mod.seq.15",
                    "isv.mod_sev.seq.5", "isv.mod_sev.seq.10","isv.mod_sev.seq.15",
                    "isv.sev.seq.5", "isv.sev.seq.10","isv.sev.seq.15")



sim_esca <- function(meteo_maille_jour, description_parcelles, info_mailles){

  cl <- makeCluster(8) # Pour paralléliser les simulations des variables écolcimatiques
  clusterExport(cl, varlist = ls(globalenv()))
  clusterEvalQ(cl, {c(library(dplyr), library(tis))})

  sim.result.list <- parLapply(cl, unique(meteo_maille_jour$ID_MAILLE), function(maille){
    info_maille <- info_mailles %>% filter(ID_MAILLE == maille)
    res_list <- list()
    for(cepage in unique(meteo_maille_jour[meteo_maille_jour$ID_MAILLE == maille,]$cepage)){ 
      res <- bilan_hydrique(meteo_maille_jour, description_parcelles, cepage, maille, 
                            info_maille$region_viticole, info_maille$RU, info_maille$altitude)
      bh <- data.frame(maille = maille, cepage = cepage, an = res$ph$an, debourrement = res$ph$C, 
                       floraison = res$ph$I)
      bh <- bh %>% filter(!an == min(bh$an)) # Enlever la première année car pas de débourrement bien simulé
      bh <- bh %>% left_join(res[["st.an.pheno"]][,param_periodes], by = "an")
      bh <- bh %>% left_join(res[["st.dormance"]][,param_periodes], by = "an", suffix = c("",".dormance"))
      bh <- bh %>% left_join(res[["st.deb.to.flo"]][,param_periodes], by = "an", suffix = c("",".deb_flo"))
      bh <- bh %>% left_join(res[["st.symptomes"]][,param_periodes], by = "an", suffix = c("",".symptomes"))
      bh <- bh %>% left_join(res[["st.deb.to.end"]][,param_periodes], by = "an", suffix = c("",".deb_end"))
      res_list[[length(res_list) + 1]] <- bh
    }
    do.call(rbind, res_list)
  })
  stopCluster(cl)
  indicateurs_climatiques <- do.call(rbind, sim.result.list)
  
  indicateurs_climatiques <- indicateurs_climatiques %>%
    left_join(info_mailles[, c("ID_MAILLE", "region_viticole", "RU")], by = c("maille" = "ID_MAILLE"))


# Simulations de l'incidence d'esca ----

  to_add <- matrix(ncol = 2)
  ages <- c(10,20,30)
  cates_esca <- c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100[")
  for(age in ages){
    for(cat in cates_esca) to_add <- rbind(to_add, matrix(c(age, cat), nrow = 1))
  }
  to_add <- as.data.frame(to_add[-1,])
  colnames(to_add) <- c("age_parcelle_estime", "esca_categoriel")
  
  to_project <- cross_join(to_add, indicateurs_climatiques)
  
  # Mettre les catégories d'esca, cépages et les régions en facteurs
  to_project$cepage <- as.factor(to_project$cepage)
  to_project$region_viticole <- as.factor(to_project$region_viticole)
  to_project$age_parcelle_estime <- as.numeric(to_project$age_parcelle_estime)
  to_project$age_10_20 <- 0
  to_project$age_10_20[to_project$age_parcelle_estime <= 20] <- 1
  to_project$age_10_20 <- as.factor(to_project$age_10_20)
  
  to_project$esca_categoriel <- factor(to_project$esca_categoriel, 
                                       levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))
  
  
  # Simulations des incidences d'esca ---- :
  
  # GLM :
  
  best_lm <- readRDS("modeles/glm_final.rds")
  pred_lm <- predict(best_lm, to_project)
  pred_lm[pred_lm < 0] <- 0
  projections <- cbind(to_project[,1:5], pred_lm)
  
  # Processus Gaussien
  library(reticulate)
  
  # transférer dans l'espace Python
  py$to_project <- r_to_py(to_project)

  
  py_run_string('
import pickle
import pandas as pd
import numpy as np
  
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
  
to_use = ["age_parcelle_estime", "et0.symptomes", "tv.deb_flo", "rr.deb_end", "VPD.symptomes", "rain.days", "tm.dormance", "rr.deb_flo", "auc_isv.symptomes", "tm.symptomes", "isv.deb_flo", "hu.dormance", "bh0.symptomes", "VPD.dormance", "VPD.deb_flo", "sum.heat.days.25.deb_flo", "rr", "sum.days.isv.faible", "sum.days.isv.fai_mod.dormance", "rr.symptomes", "sum.days.isv.sev.dormance", "auc_isv.dormance", "RU", "debourrement", "floraison"]
to_keep = list(["cepage", "region_viticole"]) + to_use
  
# Récuperer le modèle enregistré :
gp = pickle.load(open("modeles/gaussian_process_final.sav", "rb"))
pca = pickle.load(open("modeles/acp.sav", "rb"))
scaler = pickle.load(open("modeles/scaler.sav", "rb"))

to_project = pd.DataFrame(to_project)

to_project["cepage"] = to_project["cepage"].astype("category")
to_project["region_viticole"] = to_project["region_viticole"].astype("category")

to_project = to_project.assign(esca_categoriel_cont = 0)
to_project.loc[to_project.esca_categoriel == "[0;1[","esca_categoriel_cont"] = 0
to_project.loc[to_project.esca_categoriel == "[1;2[","esca_categoriel_cont"] = 1 
to_project.loc[to_project.esca_categoriel == "[2;5[","esca_categoriel_cont"] = 2
to_project.loc[to_project.esca_categoriel == "[5;10[","esca_categoriel_cont"] = 3
to_project.loc[to_project.esca_categoriel == "[10;20[","esca_categoriel_cont"] = 4
to_project.loc[to_project.esca_categoriel == "[20;100[","esca_categoriel_cont"] = 5
  
X = pd.get_dummies(to_project[to_keep])
X = X[scaler.feature_names_in_]
X_scaled = scaler.transform(X) 
X_pca = pd.DataFrame(pca.transform(X_scaled))
X_pca = X_pca.rename(columns={0: "Dim_1", 1: "Dim_2", 2: "Dim_3", 3: "Dim_4", 4: "Dim_5", 5: "Dim_6", 6: "Dim_7", 7: "Dim_8", 8: "Dim_9", 9: "Dim_10", 10: "Dim_11"})

X_pca["esca_categoriel"] = to_project["esca_categoriel_cont"]
pred = gp.predict(X_pca, return_std=True)
pred_gp = pred[0]
pred_gp[pred_gp < 0] = 0
# gp_sd = pred[1]
')
  
  pred_gp <- py_to_r(py$pred_gp)  # récupérer les projections
  # gp_sd <- py_to_r(py$gp_sd) # Si on veut récuperer l'écart type des projections pour le gp
  
  projections <- cbind(projections, pred_gp) # , gp_sd
  
  return(projections)
}



####################################### Simulation climat uniquement :


sim_climat <- function(meteo_maille_jour, description_parcelles, info_mailles){
    
    cl <- makeCluster(8) # Pour paralléliser les simulations des variables écolcimatiques
    clusterExport(cl, varlist = ls(globalenv()))
    clusterEvalQ(cl, {c(library(dplyr), library(tis))})
    
    sim.result.list <- parLapply(cl, unique(meteo_maille_jour$ID_MAILLE), function(maille){
      info_maille <- info_mailles %>% filter(ID_MAILLE == maille)
      res_list <- list()
      for(cepage in unique(meteo_maille_jour[meteo_maille_jour$ID_MAILLE == maille,]$cepage)){ 
        res <- bilan_hydrique(meteo_maille_jour, description_parcelles, cepage, maille, 
                              info_maille$region_viticole, info_maille$RU, info_maille$altitude)
        bh <- data.frame(maille = maille, cepage = cepage, an = res$ph$an, debourrement = res$ph$C, 
                         floraison = res$ph$I)
        bh <- bh %>% filter(!an == min(bh$an)) # Enlever la première année car pas de débourrement bien simulé
        bh <- bh %>% left_join(res[["st.an.pheno"]][,param_periodes], by = "an")
        bh <- bh %>% left_join(res[["st.dormance"]][,param_periodes], by = "an", suffix = c("",".dormance"))
        bh <- bh %>% left_join(res[["st.deb.to.flo"]][,param_periodes], by = "an", suffix = c("",".deb_flo"))
        bh <- bh %>% left_join(res[["st.symptomes"]][,param_periodes], by = "an", suffix = c("",".symptomes"))
        bh <- bh %>% left_join(res[["st.deb.to.end"]][,param_periodes], by = "an", suffix = c("",".deb_end"))
        res_list[[length(res_list) + 1]] <- bh
      }
      do.call(rbind, res_list)
    })
    stopCluster(cl)
    indicateurs_climatiques <- do.call(rbind, sim.result.list)
    
    indicateurs_climatiques <- indicateurs_climatiques %>%
      left_join(info_mailles[, c("ID_MAILLE", "region_viticole", "RU")], by = c("maille" = "ID_MAILLE"))

    return(indicateurs_climatiques)
}