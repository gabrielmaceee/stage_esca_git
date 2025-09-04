library(lme4)
library(dplyr)
library(openxlsx)

load("data/modelisation/observations.RData")

# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[12:47] # ne pas garder la longueur de la période = 365 ou 366
var_an <- colnames(observations)[49:84] # ne pas garder la longueur de la période
var_dormance <- colnames(observations)[85:121]
var_deb_to_flo <- setdiff(colnames(observations)[122:158] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[159:195], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[196:232]
var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)

# Joindre les Cotes-du-Rhone nord et sud en "Cotes-du-Rhone", car que 6 observations pour le sud
observations$region_viticole[observations$region_viticole %in% c("Cotes-du-Rhone nord", "Cotes-du-Rhone sud")] <- "Cotes-du-Rhone"

observations$cepage <- as.factor(observations$cepage) # Transformation des cépages en catégorielle
observations$region_viticole <- as.factor(observations$region_viticole) # Transformation des régions en catégorielle

observations$age_10_20 <- 0
observations$age_10_20[observations$age_parcelle_estime <= 20] <- 1

observations$age_20_30 <- 0
observations$age_20_30[observations$age_parcelle_estime > 20] <- 1

for(region in unique(observations$region_viticole)){ # Pour chaque région
  print(region)
  # recuperer les observations d'une seule région
  data <- observations[observations$region_viticole == region, -c(1,2,6,7) ]
  var2drop = c()
  for(var in setdiff(colnames(data), "cepage")){ # Enlever les variables avec une seule valeur (ça arrive)
    if(length(unique(data[,var]))==1) var2drop = c(var2drop, var)
  }
  data = data[, setdiff(colnames(data), var2drop)]
  
  var2keep <- c("cepage", "age_parcelle_estime", "age_10_20", "age_20_30") # variables initiales
  # Formule initiale
  form <- as.formula(paste("pourcentage_esca ~(age_10_20 + age_parcelle_estime:age_10_20 + age_20_30 + age_parcelle_estime:age_20_30 | cepage)"))
  lm_select <- lmer(form, data = data[,c(var2keep, "pourcentage_esca")])
  r2s = c()
  r2s <- c(cor(predict(lm_select, data[,var2keep]), data$pourcentage_esca)^2) # 0.1687542
  epsilone <- 0.001
  
  effects_0 <- "age_10_20 + age_parcelle_estime:age_10_20 + age_20_30 + age_parcelle_estime:age_20_30" # Effet d'intéraction avec effet fixe

  for(i in 1:10){ # garder 10 variables en plus de cépages, régions viticoles, âges
    r2_var <- c()
    for(var in setdiff(colnames(data), c(var2keep, "pourcentage_esca"))){
      effects <- paste0(var, " + ", effects_0)
      effects <- paste0("(", effects, "|cepage)") # Ajouter la variable à tester à la formule
      form <- as.formula(paste("pourcentage_esca ~", effects)) 
      lm_select <- lmer(form, data = data[,c(var2keep, var, "pourcentage_esca")]) # Modéliser
      r2_var <- c(r2_var, cor(predict(lm_select, data[,c(var2keep,var)]), data$pourcentage_esca)^2) # Calcul du r²
    }
    if(max(r2_var) > r2s[i] + epsilone){ # Si une variable permet d'obtenir un meilleur r² que l'itération précédente
      r2s <- c(r2s, max(r2_var)[1]) # Récuperer le meilleur r²
      var <- setdiff(var_tt, var2keep)[which(r2_var == max(r2_var))[1]] # Récupérer la variable maximisant le r²
      var2keep <- c(var2keep, var)
      effects_0 <- paste0(var, " + ", effects_0) # Ajouter la variable améliorant le plus le modèle 
    }
    else{break} # Si on ne peut pas améliorer significativement le r² : arrêter de chercher
    print(paste(i, ":", var2keep[i+4])) # Print la variable ajoutée
    print(paste(i, ":", r2s[i+1])) # Print le r² associé
  }
  res <- data.frame(variable = var2keep, r2 = c(rep(r2s[1],3), r2s)) # verif les 4 premières lignes = les mêmes
  write.xlsx(res, paste0("data/resultats/res_select_", region, ".xlsx"), rowNames=FALSE)
}
