# Le but de ce script est d'ajuster le modèle linéaire final, en fonction de variables choisis préalablement,
# de quantifier sa robustesse, et d'interpréter le rôle de chaque variable.


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


#### Train / test :

obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1

obs$age_10_20 <- as.factor(obs$age_10_20)
obs$age_20_30 <- as.factor(obs$age_20_30)


form <- as.formula(paste0("pourcentage_esca ~", "et0.symptomes + rr.symptomes +  rr.deb_flo +
                          VPD.symptomes:tv.deb_flo + sum.heat.days.30.dormance + sum.days.isv.mod_sev.symptomes +
                          sum.heat.days.30.dormance:isv.sev.seq.10.symptomes + tv.deb_flo +
                          age_parcelle_estime:age_10_20 + 
                          (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))

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

form <- as.formula(paste0("pourcentage_esca ~", "et0.symptomes + rr.symptomes +  rr.deb_flo +
                          VPD.symptomes:tv.deb_flo + sum.heat.days.30.dormance + sum.days.isv.mod_sev.symptomes +
                          sum.heat.days.30.dormance:isv.sev.seq.10.symptomes + tv.deb_flo +
                          age_parcelle_estime:age_10_20 + 
                          (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))

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


form <- as.formula(paste0("pourcentage_esca ~", "et0.symptomes + rr.symptomes +  rr.deb_flo +
                          VPD.symptomes:tv.deb_flo + sum.heat.days.30.dormance + sum.days.isv.mod_sev.symptomes +
                          sum.heat.days.30.dormance:isv.sev.seq.10.symptomes + tv.deb_flo +
                          age_parcelle_estime:age_10_20 + 
                          (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
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
write.xlsx(coeffs_cr, "data/resultats/glm_final.xlsx", rowNames=FALSE)

#####
# Enregistrer le modèle :
saveRDS(best_lm, "modeles/glm_final.rds") 
# Pour l'utiliser : best_lm <- readRDS("modeles/glm_final.rds")

##### Créer et enregistrer le modèle sans variables écolimatiques :

form_no_climat <- as.formula(paste0("pourcentage_esca ~", "age_parcelle_estime:age_10_20 + 
                          (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
lm_no_climat <- lmer(form_no_climat, data = obs)
saveRDS(lm_no_climat, "modeles/glm_no_climat.rds") 
# Pour l'utiliser : lm_no_climat <- readRDS("modeles/glm_no_climat.rds")



#######################################
var2plot <- c("et0.symptomes", "rr.symptomes", "rr.deb_flo", "sum.heat.days.30.dormance", 
                  "sum.days.isv.mod_sev.symptomes", "tv.deb_flo")

interactions <- c("VPD.symptomes__tv.deb_flo", "sum.heat.days.30.dormance__isv.sev.seq.10.symptomes")


a_plot <- obs[, c(var2plot, "pourcentage_esca")]

for(inter in interactions){
  parties <- strsplit(inter, "__")[[1]]
  a_plot[,inter] <- obs[,parties[1]] * obs[,parties[2]]
}

var2plot <- c(var2plot, interactions)

#### 
# Distribution des variables sélectionnées :
# Âge : 
png(filename="graphs/distribution_var_selectionnees/age_parcelle.png")
plot(density(observations$age_parcelle_estime), xlab = "Âge de la parcelle", 
     main = "Densité de l'âge de la parcelle")
dev.off()

# Les 16 premièes var-écoclimatiques :
for(var in var2plot){
  png(filename=paste0("graphs/distribution_var_selectionnees/",var,".png"))
  plot(density(a_plot[,var]), xlab = var, 
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
for(var in var2plot){
  png(filename=paste0("graphs/lien_var_selectionnees_esca/",var,".png"))
  plot(a_plot[,var], a_plot$pourcentage_esca,
       xlab = var, ylab = "Incidence d'esca",
       main = paste0("Incidence d'esca en fonction ", var))
  dev.off()
}




# 
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

best_lm <- readRDS("modeles/glm_final.rds")
get_variance_partition_lmm(best_lm)

VarCorr(best_lm)


library(partR2)

res <- partR2(
  best_lm,
  partvars = c("et0.symptomes", "rr.symptomes", "rr.deb_flo", "sum.heat.days.30.dormance", 
               "sum.days.isv.mod_sev.symptomes", "tv.deb_flo","VPD.symptomes:tv.deb_flo", 
               "sum.heat.days.30.dormance:isv.sev.seq.10.symptomes",
               "age_parcelle_estime:age_10_20"),
  data = obs,
  R2_type = "conditional",  # <--- ici on prends en compte les effets aléatoires dans la valeur du r²
  nboot = 100
)

print(res)

res_marg <- partR2(
  best_lm,
  partvars = c("et0.symptomes", "rr.symptomes", "rr.deb_flo", "sum.heat.days.30.dormance", 
               "sum.days.isv.mod_sev.symptomes", "tv.deb_flo","VPD.symptomes:tv.deb_flo", 
               "sum.heat.days.30.dormance:isv.sev.seq.10.symptomes",
               "age_parcelle_estime:age_10_20"),
  data = obs,
  R2_type = "marginal",  # <--- ici on ne prends pas en compte les effets aléatoires dans la valeur du r²
  nboot = 10
)

print(res_marg)



form <- as.formula(paste0("pourcentage_esca ~", "et0.symptomes + rr.symptomes +  rr.deb_flo +
                          VPD.symptomes:tv.deb_flo + sum.heat.days.30.dormance + sum.days.isv.mod_sev.symptomes +
                          sum.heat.days.30.dormance:isv.sev.seq.10.symptomes + tv.deb_flo +
                          age_parcelle_estime:age_10_20 + 
                          (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
lm_cr <- lmer(form, data = obs_cr)
get_variance_partition_lmm(lm_cr)


var2use <- read.xlsx("data/resultats/res_select_groupe_interactions_complet.xlsx", rowNames=FALSE)
effects_0 <- paste(var2use$variable[6:105], collapse = " + ")
form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ age_parcelle_estime:age_10_20 + (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
lm_cr_complet <- lmer(form, data = obs_cr)

get_variance_partition_lmm(lm_cr_complet)


var2use <- read.xlsx("data/resultats/res_select_groupe_interactions.xlsx", rowNames=FALSE)
effects_0 <- paste(var2use$variable[6:21], collapse = " + ")
form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ age_parcelle_estime:age_10_20 + (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
lm_cr_top_15 <- lmer(form, data = obs_cr)

get_variance_partition_lmm(lm_cr_top_15)


var2use <- read.xlsx("data/resultats/res_select_groupe_interactions.xlsx", rowNames=FALSE)
effects_0 <- paste(var2use$variable[6:54], collapse = " + ")
form <- as.formula(paste0("pourcentage_esca ~", effects_0, "+ age_parcelle_estime:age_10_20 + (1|age_10_20) + (1|region_viticole) + (1|esca_categoriel) + (1|cepage)"))
lm_cr_top_50 <- lmer(form, data = obs_cr)

get_variance_partition_lmm(lm_cr_top_50)


get_variance_partition_lmm_2 <- function(model, data, response_var) {
  y <- data[[response_var]]
  var_y <- var(y)
  
  # Prédictions
  y_hat_full  <- predict(model)                      # Fixes + aléatoires
  y_hat_fixed <- predict(model, re.form = NA)        # Fixes only
  residuals   <- residuals(model)
  
  # Variance totale expliquée par le modèle complet
  var_full  <- var(y_hat_full)
  var_fixed <- var(y_hat_fixed)
  var_rand  <- var(y_hat_full - y_hat_fixed)
  var_resid <- var(residuals)
  
  # Empêcher dépassement : renormalisation
  total_model_var <- var_fixed + var_rand
  unexplained <- var_y - total_model_var
  
  results <- data.frame(
    Component = c("Fixed effects", "Random effects", "Residual"),
    Variance = round(c(var_fixed, var_rand, unexplained), 4),
    Proportion = round(c(var_fixed, var_rand, unexplained) / var_y, 4)
  )
  
  total_row <- data.frame(
    Component = "TOTAL",
    Variance = sum(results$Variance),
    Proportion = sum(results$Proportion)
  )
  
  return(rbind(results, total_row))
}

  
res <- get_variance_partition_lmm_2(best_lm, data = obs, response_var = "pourcentage_esca")
print(res)

