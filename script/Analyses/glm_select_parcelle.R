# Le but de ce script est de sélectionner les variables et les effets d'interactions qui permettent de
# maximiser le r² du glm à effet aléatoire sur la médiane de la parcelle 
# et d'interpréter le rôle de chaque variables.

library(lme4)
library(lmerTest) # Pour interpréter l'effet de chaque variable
library(dplyr)
library(openxlsx)
library(parallel) # Pour executer la sélection de variables sur plusieurs processeurs en parallèles
library(Metrics) # RMSE


# Chargement du jeu de données avec les identifiants des parcelles
load("data/modelisation/observations_parcelles.RData")

# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[13:48] # ne pas garder la longueur de la période = 365 ou 366
var_an <- colnames(observations)[50:85] # ne pas garder la longueur de la période
var_dormance <- colnames(observations)[86:122]
var_deb_to_flo <- setdiff(colnames(observations)[123:159] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[160:196], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[197:233]
var_tt <- c("RU", "debourrement","age_parcelle_estime", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)

# Mettre les cépages et les régions en facteurs
observations$cepage <- as.factor(observations$cepage)
observations$region_viticole <- as.factor(observations$region_viticole)


# Ajouter les intéractions intéressantes :
climat = c("tm", "VPD", "sum.heat.days.30", "rr")
transpi = c("swi", "isv", "ftsw", "et0", "bhv", "tv", "auc_isv", "isv.sev.seq.10") 
climat <- c(paste0(climat, ".dormance"), paste0(climat, ".deb_flo"), paste0(climat, ".symptomes"))
transpi <- c(paste0(transpi, ".dormance"), paste0(transpi, ".deb_flo"), paste0(transpi, ".symptomes"))
interactions <- c()
for(var1 in climat){
  for(var2 in transpi){
    interactions <- c(interactions, paste0(var1, ":", var2))
  }
}


var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end, interactions)
# Sans années calendaires, avec interactions

# Division de l'âge en deux car sont lien avec l'esca et non linéaire
obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# obs[,var_tt] <- scale(obs[,var_tt]) # centrer réduire

# Récupération de la médiane d'esca par parcelle
med_parcelle <- obs %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
obs <- obs %>% left_join(med_parcelle, by = "identifiant_parcelle_analyse")
obs$med_parcelle <- as.factor(round(obs$med_parcelle, 0))


# Prendre initialement que le cépage et la région viticole et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("med_parcelle","cepage", "region_viticole", "age_parcelle_estime", "age_10_20", "age_20_30")

# Variance expliquée juste par la médiane / parcelle
form <- as.formula(paste("pourcentage_esca ~(1|med_parcelle)"))
lm_select <- lmer(form, data = obs[,c("med_parcelle", "pourcentage_esca")])
cor(predict(lm_select, obs[,var2keep]), obs$pourcentage_esca)^2 # 0.5864637

# Créer le modèle initial
form <- as.formula(paste("pourcentage_esca ~(age_parcelle_estime:age_10_20 + age_parcelle_estime:age_20_30) + (1|age_10_20) + (1|med_parcelle) + (1|region_viticole) + (1|cepage)"))
lm_select <- lmer(form, data = obs[,c(var2keep, "pourcentage_esca")])
r2s <- c(cor(predict(lm_select, obs[,var2keep]), obs$pourcentage_esca)^2) # 0.59

epsilone <- 0.0001 # Amélioration minimum du r² pour ajouter une variable

effects_0 <- "age_parcelle_estime:age_10_20 + age_parcelle_estime:age_20_30" # Effet d'intéraction avec effet fixe

for(i in 1:100){ # garder 100 variables en plus
  cl <- makeCluster(8) # Création des cluters de calcul (un par cpu)
  clusterExport(cl, varlist = ls(globalenv())) # Variable de l'environnement à exporter dans les clusters
  clusterEvalQ(cl, {c(library(dplyr), library(lme4))}) # Librairies à exporter dans les clusters
  r2.list <- parLapply(cl, setdiff(var_tt, var2keep), function(var){
    effects <- paste0(var, " + ", effects_0) # Ajouter la variable à tester
    form <- as.formula(paste("pourcentage_esca ~", effects, "+ (1|age_10_20) + (1|region_viticole) + (1|med_parcelle) + (1|cepage)")) 
    lm_select <- lmer(form, data = obs)
    pred <- predict(lm_select, obs)
    pred[pred < 0] <- 0
    cor(pred, obs$pourcentage_esca)^2 # Calcul du r² avec la variable ajoutée
  })
  stopCluster(cl)
  r2_var <- do.call(rbind, r2.list) # Tous les nouveaux r²
  if(max(r2_var) > r2s[i] + epsilone){
    r2s <- c(r2s, max(r2_var)[1]) # Récuperer le meilleur r²
    var <- setdiff(var_tt, var2keep)[which(r2_var == max(r2_var))[1]] # Récupérer la variable maximisant le r²
    var2keep <- c(var2keep, var)
    effects_0 <- paste0(var, " + ", effects_0) # Ajouter la variable améliorant le plus le modèle 
  }
  else{break}
  print(paste(i, ":", var2keep[i+6])) # Print la variable ajoutée
  print(paste(i, ":", r2s[i+1])) # Print le r² associé
}

plot(r2s, type = "l", ylab = "R²", main = "Evolution du R² en ajoutant des variables", ylim = c(0,1))

res <- data.frame(variable = var2keep, r2 = c(rep(r2s[1], 5), r2s))
write.xlsx(res, "data/resultats/res_select_mediane_interactions.xlsx", rowNames=FALSE)




####### Récupération des coeffs pour intéprétations :

# Centrer et réduire les variables continues pour pouvoir comparer leurs coeffs :
obs_cr <- obs
obs_cr[, c(9:233)] <- scale(obs_cr[,c(9:233)])

form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ (1|age_10_20) + (1|region_viticole) + (1|med_parcelle) + (1|cepage)"))
best_lm <- lmer(form, data = obs_cr)
rl <- lmerTest::lmer(form, data = obs_cr)
summary_cr <- summary(rl) 
coef(rl) 

coeffs_cr <- summary_cr$coefficients[,c(1,5)]
colnames(coeffs_cr) <- c("coeff_cr", "p_value_cr")




# Sans centrer et réduire les variables continues pour pouvoir récuperer les vrais coeffs :
form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ (1|age_10_20) + (1|region_viticole) + (1|med_parcelle) + (1|cepage)"))
best_lm <- lmer(form, data = obs)

VarCorr(best_lm)

library(performance)


# Détails des composants de variance
variance_decomposition(best_lm)


rl <- lmerTest::lmer(form, data = obs)
summary <- summary(rl)

coeffs <- summary$coefficients[,c(1,5)]
colnames(coeffs) <- c("coeff", "p_value")

intercept_med_parcelle <- coeffs$med_parcelle[,1]
names(intercept_med_parcelle) <- row.names(coeffs$med_parcelle)

intercept_cepage <- coeffs$cepage[,1]
names(intercept_cepage) <- row.names(coeffs$cepage)

intercept_region <- coeffs$region_viticole[,1]
names(intercept_region) <- row.names(coeffs$region_viticole)

###
coeffs_cr <- cbind(coeffs_cr, coeffs[,1], rownames(coeffs_cr))
colnames(coeffs_cr) <- c("coeff_cr", "p_value_cr", "coeff_reel", "variable")
coeffs_cr <- as.data.frame(coeffs_cr)

coeffs_cr <- left_join(coeffs_cr, res, by = "variable")

# Le gentil modèle linéaire c'est permis de changer l'ordre de certains effet d'intéraction :
# Ex : "rr.symptomes:et0.symptomes" -> "et0.symptomes:rr.symptomes"

for(i in 4:dim(coeffs_cr)[1]){
  if(is.na(coeffs_cr$r2[i])){
    parties <- strsplit(coeffs_cr$variable[i], ":")[[1]]
    coeffs_cr$variable[i] <- paste(parties[2], parties[1], sep = ":")
  }
}

coeffs_cr <- coeffs_cr[,-c(5)]
coeffs_cr <- left_join(coeffs_cr, res, by = "variable")
coeffs_cr <- coeffs_cr[,c("variable", "r2", "coeff_reel", "coeff_cr", "p_value_cr")]
coeffs_cr <- coeffs_cr[order(coeffs_cr$r2),]

# Character to numeric :
coeffs_cr$coeff_reel <- as.numeric(coeffs_cr$coeff_reel)
coeffs_cr$coeff_cr <- as.numeric(coeffs_cr$coeff_cr)
coeffs_cr$p_value_cr <- as.numeric(coeffs_cr$p_value_cr)

# Finalement, le r² de la saléection de variables n'est pas si important pour l'interprétation des coeffs
# On peut l'enlever
coeffs_cr <- coeffs_cr[,-2]
coeffs_cr <- coeffs_cr[order(abs(coeffs_cr$coeff_cr), decreasing = TRUE),]

# Résumé des coefficients par variable, ordonnée selon le coefficient via normalisation (cr)
write.xlsx(coeffs_cr, "data/resultats/glm_mediane_interactions.xlsx", rowNames=FALSE)


print(coeffs_cr[coeffs_cr$p_value_cr < 0.05,]$variable)


##############################################################################
# Calcul des métriques avec entrainement / test :

# Refaire le data frame juste pour être sur :

obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# On définit la proportion de données qu'on veut dans le jeu d'entraînement
train_ratio <- 0.8

train_indices <- obs %>%
  group_by(identifiant_parcelle_analyse) %>%
  group_split() %>%
  lapply(function(group) {
    n <- nrow(group)
    if(n <= 5) sample(seq_len(n), size = ceiling(n/2)) # ceiling : arrondi au supérieur, floor : arrondi à l'inférieur
    else sample(seq_len(n), size = ceiling(train_ratio * n))
  })

# Construction du jeu d'entraînement et de test
# Associer les indices à chaque groupe
train <- bind_rows(
  Map(function(group, idx) group[idx, ], group_split(obs, obs$identifiant_parcelle_analyse), train_indices)
)

test <- anti_join(obs, train, by = colnames(obs))

# Récupération de la médiane train d'esca par parcelle
med_parcelle <- train %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
train <- train %>% left_join(med_parcelle, by = "identifiant_parcelle_analyse")
train$med_parcelle <- as.factor(round(train$med_parcelle, 0))

test <- test %>% left_join(med_parcelle, by = "identifiant_parcelle_analyse")
test$med_parcelle <- as.factor(round(test$med_parcelle, 0))


form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ (1|age_10_20) + (1|region_viticole) + (1|med_parcelle) + (1|cepage)"))

lm_train <- glmer(form, data = train)

# Train :
pred_train <- predict(lm_train, train)
pred_train[pred_train < 0] <- 0
print("RMSE train :")
rmse(train$pourcentage_esca, pred_train)
print("Erreur moyenne absolue train :")
mean(abs(train$pourcentage_esca - pred_train))
print("R² train :")
cor(pred_train, train$pourcentage_esca)^2

# Test :
pred_test <- predict(lm_train, test)
pred_test[pred_test < 0] <- 0
print("RMSE test :")
rmse(test$pourcentage_esca, pred_test)
print("Erreur moyenne absolue test :")
mean(abs(test$pourcentage_esca - pred_test))
print("R² test :")
cor(pred_test, test$pourcentage_esca)^2

plot(test$pourcentage_esca, pred_test, xlab = "Incidence réelle", ylab = "Incidence prédite",
     main = "Incidence prédite en fonction de l'incidence réelle, test")

plot(test$pourcentage_esca, abs(pred_test - test$pourcentage_esca), xlab = "Incidence réelle", ylab = "Erreur de prédiction",
     main = "Erreur en fonction de l'incidence réelle, test")


######################################################################################################
### Tester avec des catégories d'incidence à la place de la médiane :
# Refaire le data frame juste pour être sur :
obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# Récup les variables et intéractions à utiliser :
var2use <- read.xlsx("data/resultats/res_select_mediane_interactions.xlsx", rowNames=FALSE)

effects_0 <- paste(var2use$variable[-c(1:6)], collapse = " + ")

form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ age_parcelle_estime:age_10_20 + age_parcelle_estime:age_20_30 + (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))

# On définit la proportion de données qu'on veut dans le jeu d'entraînement
train_ratio <- 0.8

train_indices <- obs %>%
  group_by(identifiant_parcelle_analyse) %>%
  group_split() %>%
  lapply(function(group) {
    n <- nrow(group)
    if(n <= 5) sample(seq_len(n), size = ceiling(n/2)) # ceiling : arrondi au supérieur, floor : arrondi à l'inférieur
    else sample(seq_len(n), size = ceiling(train_ratio * n))
  })

# Construction du jeu d'entraînement et de test
# Associer les indices à chaque groupe
train <- bind_rows(
  Map(function(group, idx) group[idx, ], group_split(obs, obs$identifiant_parcelle_analyse), train_indices)
)

test <- anti_join(obs, train, by = colnames(obs))

# Récupération de la médiane train d'esca par parcelle
med_parcelle <- train %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
med_parcelle$esca_categoriel <- "[0;1["
med_parcelle[med_parcelle$med_parcelle >= 1,]$esca_categoriel <- "[1;2["
med_parcelle[med_parcelle$med_parcelle >= 2,]$esca_categoriel <- "[2;5["
med_parcelle[med_parcelle$med_parcelle >= 5,]$esca_categoriel <- "[5;10["
med_parcelle[med_parcelle$med_parcelle >= 10,]$esca_categoriel <- "[10;20["
med_parcelle[med_parcelle$med_parcelle >= 20,]$esca_categoriel <- "[20;100["
med_parcelle$esca_categoriel <- factor(med_parcelle$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

train <- train %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")

test <- test %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")

lm_train <- glmer(form, data = train)

# Train :
pred_train <- predict(lm_train, train)
pred_train[pred_train < 0] <- 0
print("RMSE train :")
rmse(train$pourcentage_esca, pred_train)
print("Erreur moyenne absolue train :")
mean(abs(train$pourcentage_esca - pred_train))
print("R² train :")
cor(pred_train, train$pourcentage_esca)^2

# Test :
pred_test <- predict(lm_train, test)
pred_test[pred_test < 0] <- 0
print("RMSE test :")
rmse(test$pourcentage_esca, pred_test)
print("Erreur moyenne absolue test :")
mean(abs(test$pourcentage_esca - pred_test))
print("R² test :")
cor(pred_test, test$pourcentage_esca)^2


######################################################################################################
### Tester l'hypothèse que le climat ne permet pas d'estimer l'espérance de l'incidence d'esca de chaque 
# parcelle mais seulement la variance par années :
# Centrer les incidences en soustrayant la médiane de la parcelle (estimateur non biaisé de l'espérance et
# moins sensible aux outliers que la moyenne (surtout avec des petits échantillons))


# Refaire le data frame juste pour être sur :

obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# On définit la proportion de données qu'on veut dans le jeu d'entraînement
train_ratio <- 0.8

train_indices <- obs %>%
  group_by(identifiant_parcelle_analyse) %>%
  group_split() %>%
  lapply(function(group) {
    n <- nrow(group)
    if(n <= 5) sample(seq_len(n), size = ceiling(n/2)) # ceiling : arrondi au supérieur, floor : arrondi à l'inférieur
    else sample(seq_len(n), size = ceiling(train_ratio * n))
  })

# Construction du jeu d'entraînement et de test
# Associer les indices à chaque groupe
train <- bind_rows(
  Map(function(group, idx) group[idx, ], group_split(obs, obs$identifiant_parcelle_analyse), train_indices)
)

test <- anti_join(obs, train, by = colnames(obs))

# Récupération de la médiane d'esca par parcelle (sur l'ensemble du jeu de données)
med_parcelle <- obs %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
train <- train %>% left_join(med_parcelle, by = "identifiant_parcelle_analyse")
train$pourcentage_esca <- train$pourcentage_esca - train$med_parcelle

test <- test %>% left_join(med_parcelle, by = "identifiant_parcelle_analyse")
test$pourcentage_esca <- test$pourcentage_esca - test$med_parcelle


form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ (1|age_10_20) + (1|region_viticole) + (1|cepage)"))
lm_train <- lmer(form, data = train)

# Train :
pred_train <- predict(lm_train, train)
print("RMSE train :")
rmse(train$pourcentage_esca, pred_train)
print("Erreur moyenne absolue train :")
mean(abs(train$pourcentage_esca - pred_train))
print("R² train :")
cor(pred_train, train$pourcentage_esca)^2

# Test :
pred_test <- predict(lm_train, test)
print("RMSE test :")
rmse(test$pourcentage_esca, pred_test)
print("Erreur moyenne absolue test :")
mean(abs(test$pourcentage_esca - pred_test))
print("R² test :")
cor(pred_test, test$pourcentage_esca)^2

plot(test$pourcentage_esca, pred_test, xlab = "Incidence réelle", ylab = "Incidence prédite",
     main = "Incidence prédite en fonction de l'incidence réelle, test")

plot(test$pourcentage_esca, abs(pred_test - test$pourcentage_esca), xlab = "Incidence réelle", ylab = "Erreur de prédiction",
     main = "Erreur en fonction de l'incidence réelle, test")




###############################################################################################

# Refaire la sélection de variables mais sans swi, en enlevant les variables trop corrélées,
# et en remplaçant la médiane de l'incidence d'esca par des catégories d'incidence :

###############################################################################################

var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
var_tt <- setdiff(var_tt, c("swi", "swi.dormance", "swi.deb_flo", "swi.symptomes", "swi.deb_end"))
to_test <- observations[ ,var_tt]

# Ajouter les intéractions intéressantes :
climat = c("tm", "VPD", "sum.heat.days.30", "rr")
transpi = c("isv", "ftsw", "et0", "bhv", "tv", "auc_isv", "isv.sev.seq.10") 
climat <- c(paste0(climat, ".dormance"), paste0(climat, ".deb_flo"), paste0(climat, ".symptomes"))
transpi <- c(paste0(transpi, ".dormance"), paste0(transpi, ".deb_flo"), paste0(transpi, ".symptomes"))
interactions <- c()
for(var1 in climat){
  for(var2 in transpi){
    interactions <- c(interactions, paste0(var1, ":", var2))
  }
}

var_tt <- c(var_tt, interactions)

for(inter in interactions){
  parties <- strsplit(inter, ":")[[1]]
  to_test[,inter] <- to_test[,parties[1]] * to_test[,parties[2]]
}


# Enlever les variables trop corrélées :
# Calculer les corrélations de chaque variable à l'incidence d'esca
corr_esca <- c()
for(var in var_tt){
  corr_esca <- c(corr_esca, cor(observations$pourcentage_esca, to_test[, var]))
}

score_var <- data.frame(variable = var_tt, corr_esca = abs(corr_esca))

to_rm <- c()
for(var1 in var_tt){
  if(!var1 %in% to_rm){
    for(var2 in setdiff(var_tt, c(var1, to_rm))){
      if(cor(to_test[,var1], to_test[,var2]) > 0.7){
        if(score_var[score_var$variable == var1, "corr_esca"] > score_var[score_var$variable == var1, "corr_esca"]){
          to_rm <- c(to_rm, var1)
          break
        }
        to_rm <- c(to_rm, var2)
      }
    }
  }
}

var_to_test <- setdiff(var_tt, to_rm)
var_to_test



# Division de l'âge en deux car sont lien avec l'esca et non linéaire
obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# obs[,var_tt] <- scale(obs[,var_tt]) # centrer réduire

# Récupération de la médiane d'esca par parcelle
med_parcelle <- obs %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
med_parcelle$esca_categoriel <- "[0;1["
med_parcelle[med_parcelle$med_parcelle >= 1,]$esca_categoriel <- "[1;2["
med_parcelle[med_parcelle$med_parcelle >= 2,]$esca_categoriel <- "[2;5["
med_parcelle[med_parcelle$med_parcelle >= 5,]$esca_categoriel <- "[5;10["
med_parcelle[med_parcelle$med_parcelle >= 10,]$esca_categoriel <- "[10;20["
med_parcelle[med_parcelle$med_parcelle >= 20,]$esca_categoriel <- "[20;100["
med_parcelle$esca_categoriel <- factor(med_parcelle$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

obs <- obs %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")


# Prendre initialement que le cépage et la région viticole et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("esca_categoriel", "cepage", "region_viticole", "age_parcelle_estime", "age_10_20")


# Variance expliquée juste par la médiane / parcelle
form <- as.formula(paste("pourcentage_esca ~(1|esca_categoriel)"))
lm_select <- lmer(form, data = obs[,c("esca_categoriel", "pourcentage_esca")])
cor(predict(lm_select, obs[,var2keep]), obs$pourcentage_esca)^2 # 0.55

# Créer le modèle initial
form <- as.formula(paste("pourcentage_esca ~age_parcelle_estime:age_10_20  + (1|age_10_20) + (1|esca_categoriel) + (1|region_viticole) + (1|cepage)"))
lm_select <- lmer(form, data = obs[,c(var2keep, "pourcentage_esca")])
r2s <- c(cor(predict(lm_select, obs[,var2keep]), obs$pourcentage_esca)^2) # 0.56

epsilone <- 0.0001 # Amélioration minimum du r² pour ajouter une variable

effects_0 <- "age_parcelle_estime:age_10_20" # Effet d'intéraction avec effet fixe

for(i in 1:length(var_to_test)){
  cl <- makeCluster(8) # Création des cluters de calcul (un par cpu)
  clusterExport(cl, varlist = ls(globalenv())) # Variable de l'environnement à exporter dans les clusters
  clusterEvalQ(cl, {c(library(dplyr), library(lme4))}) # Librairies à exporter dans les clusters
  r2.list <- parLapply(cl, setdiff(var_to_test, var2keep), function(var){
    effects <- paste0(var, " + ", effects_0) # Ajouter la variable à tester
    form <- as.formula(paste("pourcentage_esca ~", effects, "+ (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)")) 
    lm_select <- lmer(form, data = obs)
    pred <- predict(lm_select, obs)
    pred[pred < 0] <- 0
    cor(pred, obs$pourcentage_esca)^2 # Calcul du r² avec la variable ajoutée
  })
  stopCluster(cl)
  r2_var <- do.call(rbind, r2.list) # Tous les nouveaux r²
  if(max(r2_var) > r2s[i] + epsilone){
    r2s <- c(r2s, max(r2_var)[1]) # Récuperer le meilleur r²
    var <- setdiff(var_to_test, var2keep)[which(r2_var == max(r2_var))[1]] # Récupérer la variable maximisant le r²
    var2keep <- c(var2keep, var)
    effects_0 <- paste0(var, " + ", effects_0) # Ajouter la variable améliorant le plus le modèle 
  }
  else{break}
  print(paste(i, ":", var2keep[i+5])) # Print la variable ajoutée
  print(paste(i, ":", r2s[i+1])) # Print le r² associé
}

plot(r2s, type = "l", ylab = "R²", main = "Evolution du R² en ajoutant des variables")

res <- data.frame(variable = var2keep, r2 = c(rep(r2s[1], 4), r2s))
# write.xlsx(res, "data/resultats/res_select_groupe_interactions.xlsx", rowNames=FALSE)
write.xlsx(res, "data/resultats/res_select_groupe_interactions_complet.xlsx", rowNames=FALSE)

# Calcul du point d'inflexion du r² en focntion de la sélection de variable : 
var2use <- read.xlsx("data/resultats/res_select_groupe_interactions.xlsx", rowNames=FALSE)
r2s <- var2use$r2[-c(1:5)]

d1 <- diff(r2s)        # dérivée 1 (pente)
d2 <- diff(d1)               # dérivée 2 (changement de pente)
plot(c(0,d1), type = "b", main = "Amélioration du R²", 
     xlab = "Position", ylab = "ΔR²")
plot(c(0,0,d2), type = "b", main = "Inflexion de R² (changement de pente)", 
     xlab = "Position", ylab = "Δ(ΔR²)")

library(inflection)
res <- findiplist(1:length(r2s), r2s, index=1)
plot(r2s, type = "b")
abline(v = res$iplast, col = "red", lty = 2)



#### 
# Distribution des variables sélectionnées :
# Âge : 
png(filename="graphs/distribution_var_selectionnees/age_parcelle.png")
plot(density(observations$age_parcelle_estime), xlab = "Âge de la parcelle", 
     main = "Densité de l'âge de la parcelle")
dev.off()

# Les 16 premièes var-écoclimatiques :
for(var in var2use$variable[6:21]){
  png(filename=paste0("graphs/distribution_var_selectionnees/",var,".png"))
  plot(density(to_test[,var]), xlab = var, 
       main = paste0("Densité de ", var))
  dev.off()
}

# Lien des variables sélectionnées avec l'esca :
# Âge :
png(filename="graphs/lien_var_selectionnees_esca/age_parcelle.png")
plot(observations$age_parcelle_estime, observations$pourcentage_esca,
     xlab = "Âge de la parcelle", ylab = "Incidence d'esca",
     main = "Incidence d'esca en fonction de l'âge de la parcelle")
dev.off()

# Les 20 premièes var-écoclimatiques :
for(var in var2use$variable[6:21]){
  png(filename=paste0("graphs/lien_var_selectionnees_esca/",var,".png"))
  plot(to_test[,var], observations$pourcentage_esca,
       xlab = var, ylab = "Incidence d'esca",
       main = paste0("Incidence d'esca en fonction ", var))
  dev.off()
}


#### Train / test :

obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# Récup les variables et intéractions à utiliser :
var2use <- read.xlsx("data/resultats/res_select_groupe_interactions.xlsx", rowNames=FALSE)
effects_0 <- paste(var2use$variable[6:21], collapse = " + ")
form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ age_parcelle_estime:age_10_20 + (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))

# On définit la proportion de données qu'on veut dans le jeu d'entraînement
train_ratio <- 0.8

train_indices <- obs %>%
  group_by(identifiant_parcelle_analyse) %>%
  group_split() %>%
  lapply(function(group) {
    n <- nrow(group)
    if(n <= 5) sample(seq_len(n), size = ceiling(n/2)) # ceiling : arrondi au supérieur, floor : arrondi à l'inférieur
    else sample(seq_len(n), size = ceiling(train_ratio * n))
  })

# Construction du jeu d'entraînement et de test
# Associer les indices à chaque groupe
train <- bind_rows(
  Map(function(group, idx) group[idx, ], group_split(obs, obs$identifiant_parcelle_analyse), train_indices)
)

test <- anti_join(obs, train, by = colnames(obs))

# Récupération de la médiane train d'esca par parcelle
med_parcelle <- train %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
med_parcelle$esca_categoriel <- "[0;1["
med_parcelle[med_parcelle$med_parcelle >= 1,]$esca_categoriel <- "[1;2["
med_parcelle[med_parcelle$med_parcelle >= 2,]$esca_categoriel <- "[2;5["
med_parcelle[med_parcelle$med_parcelle >= 5,]$esca_categoriel <- "[5;10["
med_parcelle[med_parcelle$med_parcelle >= 10,]$esca_categoriel <- "[10;20["
med_parcelle[med_parcelle$med_parcelle >= 20,]$esca_categoriel <- "[20;100["
med_parcelle$esca_categoriel <- factor(med_parcelle$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

train <- train %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")

test <- test %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")

lm_train <- glmer(form, data = train)

# Train :
pred_train <- predict(lm_train, train)
pred_train[pred_train < 0] <- 0
print("RMSE train :")
rmse(train$pourcentage_esca, pred_train)
print("Erreur moyenne absolue train :")
mean(abs(train$pourcentage_esca - pred_train))
print("R² train :")
cor(pred_train, train$pourcentage_esca)^2

# Test :
pred_test <- predict(lm_train, test)
pred_test[pred_test < 0] <- 0
print("RMSE test :")
rmse(test$pourcentage_esca, pred_test)
print("Erreur moyenne absolue test :")
mean(abs(test$pourcentage_esca - pred_test))
print("R² test :")
cor(pred_test, test$pourcentage_esca)^2

plot(test$pourcentage_esca, pred_test, xlab = "Incidence réelle", ylab = "Incidence prédite",
     main = "Incidence prédite en fonction de l'incidence réelle, test", 
     asp=1)
abline(a = 0, b = 1, col = "red", lty = 2)

plot(test$pourcentage_esca, abs(pred_test - test$pourcentage_esca), xlab = "Incidence réelle", ylab = "Erreur de prédiction",
     main = "Erreur en fonction de l'incidence réelle, test")


####### Récupération des coeffs pour intéprétations :


obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)

# Centrer et réduire les variables continues pour pouvoir comparer leurs coeffs :
obs_cr <- obs
obs_cr[, c(9:233)] <- scale(obs_cr[,c(9:233)])


# Récupération de la médiane train d'esca par parcelle
med_parcelle <- obs_cr %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
med_parcelle$esca_categoriel <- "[0;1["
med_parcelle[med_parcelle$med_parcelle >= 1,]$esca_categoriel <- "[1;2["
med_parcelle[med_parcelle$med_parcelle >= 2,]$esca_categoriel <- "[2;5["
med_parcelle[med_parcelle$med_parcelle >= 5,]$esca_categoriel <- "[5;10["
med_parcelle[med_parcelle$med_parcelle >= 10,]$esca_categoriel <- "[10;20["
med_parcelle[med_parcelle$med_parcelle >= 20,]$esca_categoriel <- "[20;100["
med_parcelle$esca_categoriel <- factor(med_parcelle$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

obs_cr <- obs_cr %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")


# Récup les variables et intéractions à utiliser :
var2use <- read.xlsx("data/resultats/res_select_groupe_interactions.xlsx", rowNames=FALSE)
effects_0 <- paste(var2use$variable[6:21], collapse = " + ")
form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ age_parcelle_estime:age_10_20 + (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
best_lm <- lmer(form, data = obs_cr)
rl <- lmerTest::lmer(form, data = obs_cr)
summary_cr <- summary(rl) 
coef(rl) 

coeffs_cr <- summary_cr$coefficients[,c(1,5)]
colnames(coeffs_cr) <- c("coeff_cr", "p_value_cr")


# Sans centrer et réduire les variables continues pour pouvoir récuperer les vrais coeffs :

# Récupération de la médiane d'esca par parcelle
med_parcelle <- obs %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
med_parcelle$esca_categoriel <- "[0;1["
med_parcelle[med_parcelle$med_parcelle >= 1,]$esca_categoriel <- "[1;2["
med_parcelle[med_parcelle$med_parcelle >= 2,]$esca_categoriel <- "[2;5["
med_parcelle[med_parcelle$med_parcelle >= 5,]$esca_categoriel <- "[5;10["
med_parcelle[med_parcelle$med_parcelle >= 10,]$esca_categoriel <- "[10;20["
med_parcelle[med_parcelle$med_parcelle >= 20,]$esca_categoriel <- "[20;100["
med_parcelle$esca_categoriel <- factor(med_parcelle$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

obs <- obs %>% left_join(med_parcelle[, c("identifiant_parcelle_analyse", "esca_categoriel")], by = "identifiant_parcelle_analyse")


form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
best_lm <- lmer(form, data = obs)

rl <- lmerTest::lmer(form, data = obs)
summary <- summary(rl)

coeffs <- summary$coefficients[,c(1,5)]
colnames(coeffs) <- c("coeff", "p_value")

coeffs_cate <- coef(best_lm)
intercept_esca_categoriel <- coeffs_cate$esca_categoriel[,1]
names(intercept_esca_categoriel) <- row.names(coeffs_cate$esca_categoriel)

intercept_cepage <- coeffs_cate$cepage[,1]
names(intercept_cepage) <- row.names(coeffs_cate$cepage)

intercept_region <- coeffs_cate$region_viticole[,1]
names(intercept_region) <- row.names(coeffs_cate$region_viticole)

###
coeffs_cr <- cbind(coeffs_cr, coeffs[,1], rownames(coeffs_cr))
colnames(coeffs_cr) <- c("coeff_cr", "p_value_cr", "coeff_reel", "variable")
coeffs_cr <- as.data.frame(coeffs_cr)


coeffs_cr <- coeffs_cr[,c("variable", "coeff_reel", "coeff_cr", "p_value_cr")]


# Character to numeric :
coeffs_cr$coeff_reel <- as.numeric(coeffs_cr$coeff_reel)
coeffs_cr$coeff_cr <- as.numeric(coeffs_cr$coeff_cr)
coeffs_cr$p_value_cr <- as.numeric(coeffs_cr$p_value_cr)

coeffs_cr <- coeffs_cr[order(abs(coeffs_cr$coeff_cr), decreasing = TRUE),]

# Résumé des coefficients par variable, ordonnée selon le coefficient via normalisation (cr)
write.xlsx(coeffs_cr, "data/resultats/glm_groupe_interactions.xlsx", rowNames=FALSE)


print(coeffs_cr[coeffs_cr$p_value_cr < 0.05,]$variable)





### Répartition de la variance :

get_variance_partition_lmm <- function(model) {
  # Variance des effets aléatoires + résidu
  vc <- as.data.frame(VarCorr(model))
  re_vars <- setNames(vc$vcov, vc$grp)  # nom
  
  # Variance expliquée par les effets fixes
  X <- model.matrix(model)
  beta <- fixef(model)
  fixed_var <- var(as.vector(X %*% beta))
  
  total_var <- sum(re_vars) + fixed_var
  
  prop <- c(re_vars, fixed = fixed_var)
  prop_ratio <- round(prop / total_var, 3)
  
  return(data.frame(
    Component = names(prop_ratio),
    Variance = round(prop, 3),
    Proportion = prop_ratio
  ))
}


get_variance_partition_lmm(best_lm)

VarCorr(best_lm)

### Détails des effets fixes :

library(partR2)

res <- partR2(
  best_lm,
  partvars = var2use$variable[6:21],
  data = obs,
  R2_type = "conditional",
  nboot = 100
)
print(res)
