library(lme4)
# library(MuMIn) # Pour réuperer les r² avec les glm lme4
library(dplyr)
library(ggplot2)
library(Metrics) # RMSE
library(openxlsx)

load("data/modelisation/observations.RData")

# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[13:48] # ne pas garder la longueur de la période = 365 ou 366
var_an <- colnames(observations)[50:85] # ne pas garder la longueur de la période
var_dormance <- colnames(observations)[86:122]
var_deb_to_flo <- setdiff(colnames(observations)[123:159] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[160:196], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[197:233]
var_tt <- c("RU", "debourrement","age_parcelle_estime", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)

observations$cepage <- as.factor(observations$cepage)
observations$region_viticole <- as.factor(observations$region_viticole)
# features <- observations[ ,-c(1, 2, 5, 6)]
# features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")
i_train = sample(1:dim(observations)[1], replace = FALSE, size = 4400)
# train = features[i_train ,]



######################################################################

var_tt <- c("cepage", "region_viticole", "age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
r2s <- c() # adj.r2 <- r2 - (1 - r2) * (k / (n-k-1))
vars <- c()
rmses <- c()
for(var in var_tt){
  vars <- c(vars, var)
  lm <- lm(formula = paste0("pourcentage_esca~",var), data = observations[, c("pourcentage_esca", var)])
  r2s <- c(r2s, summary(lm)$r.squared)
  rmses <- c(rmses, rmse(observations$pourcentage_esca, predict(lm, var=observations[,var])))
}
res <- data.frame(variable = vars, r2 = r2s, rmse = rmses)
write.xlsx(res[order(res$rmse, decreasing = FALSE),] , "data/resultats/r2_rmse_par_variable.xlsx", rowNames=FALSE)

#####################################################################

var_tt <- c("age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
observations[,var_tt] <- scale(observations[,var_tt]) # centrer réduire


# Prendre initialement que le cépage et la région viticole et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("cepage", "region_viticole")
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(var2keep, "pourcentage_esca")])

r2s <- c(cor(predict(lm_select, observations[,var2keep]), observations$pourcentage_esca)^2) # 0.24
epsilone <- 0.001

for(i in 1:50){ # garder 50 variables en plus de cépages et régions viticoles
  r2_var <- c()
  for(var in setdiff(var_tt, var2keep)){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(var2keep, var, "pourcentage_esca")])
    r2_var <- c(r2_var, cor(predict(lm_select, observations[,c(var2keep,var)]), observations$pourcentage_esca)^2)
}
  if(max(r2_var) > r2s[i] + epsilone){
    r2s <- c(r2s, max(r2_var))
    var2keep <- c(var2keep, setdiff(var_tt, var2keep)[which(r2_var == max(r2_var))])
  }
  else{break}
  print(paste(i, ":", var2keep[i+2]))
}
plot(r2s, type = "l")

res <- data.frame(variable = var2keep, r2 = c(r2s[1], r2s))
write.xlsx(res, "data/resultats/res_select_glm_r2.xlsx", rowNames=FALSE)

########################################################################################################
# Diff avec juste avant : séparation âge = 10-15, 15-25, 25-30, et ajout des effets fixes cépages pour toutes les variables
var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
obs <- observations
obs$age_10_15 <- 0
obs$age_10_15[obs$age_parcelle_estime <= 15] <- 1

obs$age_15_25 <- 0
obs$age_15_25[obs$age_parcelle_estime > 15 & obs$age_parcelle_estime <= 25] <- 1

obs$age_25_30 <- 0
obs$age_25_30[obs$age_parcelle_estime > 25] <- 1
obs[,var_tt] <- scale(obs[,var_tt]) # centrer réduire


# Prendre initialement que le cépage, région viticole, et les catégories d'âges et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("cepage", "region_viticole", "age_parcelle_estime", "age_10_15", "age_15_25", "age_25_30")
form <- as.formula(paste("pourcentage_esca ~(age_parcelle_estime:age_10_15 + age_parcelle_estime:age_15_25 + age_parcelle_estime:age_25_30|cepage) + (1|region_viticole)")) # verif que + (1|cepage)
i_train = sample(1:dim(observations)[1], replace = FALSE, size = 4400)
lm_select <- lmer(form, data = obs[i_train,c(var2keep, "pourcentage_esca")])

r2s <- c(cor(predict(lm_select, obs[-i_train,var2keep]), obs$pourcentage_esca[-i_train])^2) # 0.2986534
epsilone <- 0.001

effects_0 <- "age_parcelle_estime:age_10_15 + age_parcelle_estime:age_15_25 + age_parcelle_estime:age_25_30" # Effet d'intéraction sans effet fixe


for(i in 1:50){ # garder 50 variables en plus de cépages, régions viticoles, âges
  r2_var <- c()
  for(var in setdiff(var_tt, var2keep)){
    effects <- paste0(var, " + ", effects_0)
    effects <- paste0("(", effects, "|cepage)")
    form <- as.formula(paste("pourcentage_esca ~", effects, " + (1|region_viticole)")) # verif que + (1|cepage)
    lm_select <- lmer(form, data = obs[i_train,c(var2keep, var, "pourcentage_esca")])
    r2_var <- c(r2_var, cor(predict(lm_select, obs[-i_train,c(var2keep,var)]), obs$pourcentage_esca[-i_train])^2)
  }
  if(max(r2_var) > r2s[i] + epsilone){
    r2s <- c(r2s, max(r2_var)[1])
    var <- setdiff(var_tt, var2keep)[which(r2_var == max(r2_var))]
    var2keep <- c(var2keep, var)
    effects_0 <- paste0(var, " + ", effects_0)
  }
  else{break}
  print(paste(i, ":", var2keep[i+6]))
}
plot(r2s, type = "l")

res <- data.frame(variable = var2keep, r2 = c(rep(r2s[1],5), r2s))# verif les 6 premières lignes = les mêmes
write.xlsx(res, "data/resultats/res_select_age3_cepage.xlsx", rowNames=FALSE)


########################################################################################################

# Diff avec juste avant : pas de train / test, séparation âge = 10-20, 20-30, et ajout des effets fixes catégories d'âge

var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1

obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1
obs[,var_tt] <- scale(obs[,var_tt]) # centrer réduire


# Prendre initialement que le cépage et la région viticole et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("cepage", "region_viticole", "age_parcelle_estime", "age_10_20", "age_20_30")
form <- as.formula(paste("pourcentage_esca ~(age_10_20 + age_parcelle_estime:age_10_20 + age_20_30 + age_parcelle_estime:age_20_30) + (1|region_viticole)"))
# j'ai oublié le "|cepage" ici, mais pas important pour la sélection de variables
lm_select <- lmer(form, data = obs[,c(var2keep, "pourcentage_esca")])

r2s <- c(cor(predict(lm_select, obs[,var2keep]), obs$pourcentage_esca)^2) # 0.2986534
epsilone <- 0.001

effects_0 <- "age_10_20 + age_parcelle_estime:age_10_20 + age_20_30 + age_parcelle_estime:age_20_30" # Effet d'intéraction avec effet fixe


for(i in 1:12){ # garder 20 variables en plus de cépages et régions viticoles
  r2_var <- c()
  for(var in setdiff(var_tt, var2keep)){
    effects <- paste0(var, " + ", effects_0)
    effects <- paste0("(", effects, "|cepage)")
    form <- as.formula(paste("pourcentage_esca ~", effects, " + (1|region_viticole)")) # verif que + (1|cepage)
    lm_select <- lmer(form, data = obs[,c(var2keep, var, "pourcentage_esca")])
    r2_var <- c(r2_var, cor(predict(lm_select, obs[,c(var2keep,var)]), obs$pourcentage_esca)^2)
  }
  if(max(r2_var) > r2s[i] + epsilone){
    r2s <- c(r2s, max(r2_var)[1])
    var <- setdiff(var_tt, var2keep)[which(r2_var == max(r2_var))]
    var2keep <- c(var2keep, var)
    effects_0 <- paste0(var, " + ", effects_0)
  }
  else{break}
  print(paste(i, ":", var2keep[i+5]))
}
plot(r2s, type = "l")

res <- data.frame(variable = var2keep, r2 = c(rep(r2s[1],4), r2s)) # verif les 5 premières lignes = les mêmes


########################################################################################################

# Pareil avec rmse :
var2keep <- c("cepage", "region_viticole")
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(var2keep, "pourcentage_esca")])

rmses <- c(Metrics::rmse(predict(lm_select, observations[,var2keep]), observations$pourcentage_esca))
epsilone <- 0.005

for(i in 1:50){ # garder 50 variables en plus de cépages et régions viticoles
  rmse_var <- c()
  for(var in setdiff(var_tt, var2keep)){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(var2keep, var, "pourcentage_esca")])
    rmse_var <- c(rmse_var, Metrics::rmse(predict(lm_select, observations[,c(var2keep,var)]), observations$pourcentage_esca))
    
  }
  if(min(rmse_var) < rmses[i] - epsilone){
    rmses <- c(rmses, min(rmse_var))
    var2keep <- c(var2keep, setdiff(var_tt, var2keep)[which(rmse_var == min(rmse_var))])
  }
  else{break}
  print(paste(i, ":", var2keep[i+2]))
}
plot(rmses, type = "l")
res <- data.frame(variable = var2keep, rmse = c(rmses[1], rmses))
write.xlsx(res, "data/resultats/res_select_glm_rmse.xlsx", rowNames=FALSE)

# "et0.symptomes", "ftsw.dormance", "sum.days.isv.sev.symptomes"
# sont sélectionnés par r2 et rmse mais beaucoup moins importants que région viticole et cépage

## Les variables qui handicapent le modèle :

var2rm <- c()
var_tt <- c("cepage", "region_viticole", "age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(var_tt, "pourcentage_esca")])
rmse1 <- c(Metrics::rmse(predict(lm_select, observations[,var_tt]), observations$pourcentage_esca))
rmses <- c()

for(var in setdiff(var_tt, c("cepage", "region_viticole"))){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(setdiff(var_tt, var), "pourcentage_esca")])
    rmse_var <- Metrics::rmse(predict(lm_select, observations[,setdiff(var_tt, var)]), observations$pourcentage_esca)
    if(rmse_var < rmse1){
    rmses <- c(rmses, rmse_var)
    var2rm <- c(var2rm, var)
    }
}
res <- data.frame(variable = var2rm, rmse = rmses)


## Ou faire l'inverse : partir du lm totale et enlevé ce qui ne change pas grand chose au modèle

var2rm <- c()
var_tt <- c("cepage", "region_viticole", "age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(var_tt, "pourcentage_esca")])
rmses <- c(Metrics::rmse(predict(lm_select, observations[,var_tt]), observations$pourcentage_esca))

for(i in 1:(length(var_tt)-2)){
  rmse_vars <- c()
  var2test = setdiff(var_tt, c(var2rm, "cepage", "region_viticole"))
  for(var in var2test){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations[,c(setdiff(var_tt, c(var, var2rm)), "pourcentage_esca")])
    rmse_vars <- c(rmse_vars,Metrics::rmse(predict(lm_select, observations[,setdiff(var_tt, c(var, var2rm))]), observations$pourcentage_esca))
  }
  rmses <- c(rmses, max(rmse_vars))
  var2rm <- c(var2rm, var2tests[which(rmse_vars == max(rmse_vars))])
}
res <- data.frame(variable = c("cepage", "region_viticole", var2rm), rmse = c(rmses[1], rmses[1], rmses[2:219]))
write.xlsx(res, "data/resultats/res_select_inv_glm_rmse.xlsx", rowNames=FALSE)





###########################################################################################################
## Avec effet aléatoire moyenne parcelle

load("data/modelisation/observations_parcelles.RData")
library(parallel)
var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
obs <- observations
obs$age_10_20 <- 0
obs$age_10_20[obs$age_parcelle_estime <= 20] <- 1
obs$age_20_30 <- 0
obs$age_20_30[obs$age_parcelle_estime > 20] <- 1
obs[,var_tt] <- scale(obs[,var_tt]) # centrer réduire
observations$cepage <- as.factor(observations$cepage)
observations$region_viticole <- as.factor(observations$region_viticole)

med_parcelle <- obs %>% group_by(identifiant_parcelle_analyse) %>% summarise(med_parcelle = median(pourcentage_esca))
obs <- obs %>% left_join(med_parcelle, by = "identifiant_parcelle_analyse")
obs$med_parcelle <- as.factor(round(obs$med_parcelle, 0))


# Prendre initialement que le cépage et la région viticole et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("med_parcelle","cepage", "region_viticole", "age_parcelle_estime", "age_10_20", "age_20_30")
form <- as.formula(paste("pourcentage_esca ~(age_10_20 + age_parcelle_estime:age_10_20 + age_20_30 + age_parcelle_estime:age_20_30) + (1|region_viticole) + (1|med_parcelle)"))
# j'ai oublié le "|cepage" ici, mais pas important pour la sélection de variables
lm_select <- lmer(form, data = obs[,c(var2keep, "pourcentage_esca")])

r2s <- c(cor(predict(lm_select, obs[,var2keep]), obs$pourcentage_esca)^2) # 0.59
epsilone <- 0.001

effects_0 <- "age_10_20 + age_parcelle_estime:age_10_20 + age_20_30 + age_parcelle_estime:age_20_30" # Effet d'intéraction avec effet fixe

for(i in 1:50){ # garder 50 variables en plus
  cl <- makeCluster(8)
  clusterExport(cl, varlist = ls(globalenv()))
  clusterEvalQ(cl, {c(library(dplyr), library(lme4))})
  # for(var in setdiff(var_tt, var2keep)){
  r2.list <- parLapply(cl, setdiff(var_tt, var2keep), function(var){
    effects <- paste0(var, " + ", effects_0)
    effects <- paste0("(", effects, "|cepage)")
    form <- as.formula(paste("pourcentage_esca ~", effects, " + (1|region_viticole) + (1|med_parcelle)")) # verif que + (1|cepage)
    lm_select <- lmer(form, data = obs[,c(var2keep, var, "pourcentage_esca")])
    pred <- predict(lm_select, obs[,c(var2keep,var)])
    pred[pred < 0] <- 0
    cor(pred, obs$pourcentage_esca)^2
  })
  stopCluster(cl)
  r2_var <- do.call(rbind, r2.list)
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
plot(r2s, type = "l")

res <- data.frame(variable = var2keep, r2 = c(rep(r2s[1], 5), r2s))
write.xlsx(res, "data/resultats/res_select_glm_mediane.xlsx", rowNames=FALSE)



##########################################################################################################
##### Modèle pénalisés :

library(glmnet)

features <- observations[ ,var_tt]
features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")

ridge <- glmnet(x = features, y = observations$pourcentage_esca, alpha = 0)
Metrics::rmse(predict(ridge, features), observations$pourcentage_esca)

cv_ridge <- cv.glmnet(features, observations$pourcentage_esca, alpha = 0)
ridge <- glmnet(features, observations$pourcentage_esca, alpha = 0, lambda = cv_ridge$lambda.min)
Metrics::rmse(predict(ridge, features), observations$pourcentage_esca)

lasso <- glmnet(x = as.matrix(observations[ ,-c(1, 2, 5, 6)]), y = observations$pourcentage_esca, alpha = 1 )
Metrics::rmse(predict(lasso, as.matrix(observations[ ,-c(1, 2, 5, 6)])), observations$pourcentage_esca)

lasso <- glmnet(x = features, y = observations$pourcentage_esca, alpha = 1) 
Metrics::rmse(predict(lasso, features), observations$pourcentage_esca)

# , family = poisson(link = "log") = pas top

cv_lasso <- cv.glmnet(features, observations$pourcentage_esca, alpha = 1)
lasso <- glmnet(features, observations$pourcentage_esca, alpha = 1, lambda = cv_lasso$lambda.min)
Metrics::rmse(predict(lasso, features), observations$pourcentage_esca)
# 5.017208
res_lasso = data.frame(var = var_tt[-c(1,2)], beta = abs(lasso$beta[36:253]))


score_var <- read.xlsx("data/resultats/r2_rmse_par_variable.xlsx")
vars <- var_tt[-c(1,2)]
to_rm <- c()
for(var1 in vars){
  if(!var1 %in% to_rm){
    for(var2 in setdiff(vars, c(var1, to_rm))){
      if(cor(observations[,var1], observations[,var2]) > 0.8){
        if(score_var[score_var$variable == var1, "rmse"] > score_var[score_var$variable == var1, "rmse"]){
          to_rm <- c(to_rm, var1)
          break
        }
        to_rm <- c(to_rm, var2)
      }
    }
  }
}

to_keep <- c("cepage", "region_viticole",setdiff(vars, to_rm))

library(glmmLasso)
fixed_formula <- reformulate(to_keep[-c(1,2)], response = "pourcentage_esca")
glmLasso <- glmmLasso(fix = fixed_formula, rnd=list(cepage=~1, region_viticole = ~1), 
                      data = observations[ ,c(to_keep, "pourcentage_esca")], lambda = 400,
                      family = poisson(link = "log")) #family = gaussian(link = "identity"))
Metrics::rmse(predict(glmLasso, observations[ ,c(to_keep, "pourcentage_esca")]), observations$pourcentage_esca)
# 6.153954


print("Accuracy train : ")
pred_train <- round(predict(glmLasso, observations[,c(to_keep, "esca_categoriel")]),0)
table(pred_train == observations$esca_categoriel)[1]/ length(pred_train)
print("Train : Table des correspondances (prédictions / réelles)")
table(pred_train, observations$esca_categoriel)



# Grille de lambdas à tester
lambda_grid <- c(200, 300, 400, 500) # 1, 5, 10, 20, 50, 100, 

# Nombre de folds
k <- 5
folds <- sample(rep(1:k, length.out = nrow(observations)))


# Stocker les erreurs moyennes
cv_errors <- numeric(length(lambda_grid))

for (i in seq_along(lambda_grid)) {
  lambda <- lambda_grid[i]
  errors <- numeric(k)
  
  for (j in 1:k) {
    train_data <- observations[folds != j, ]
    test_data <- observations[folds == j, ]
    
    # Ajuster le modèle
    tryCatch({
      fit <- glmmLasso(
        fix = fixed_formula,
        rnd = list(cepage = ~1, region_viticole = ~1),
        data = train_data[,c(to_keep, "pourcentage_esca")],
        lambda = lambda,
        family = poisson(link = "log")
      )
      
      # Prédictions sur le test set
      X_test <- model.matrix(fixed_formula, test_data)
      eta_test <- X_test %*% fit$coefficients
      mu_test <- exp(eta_test)
      
      # Calcul de la deviance du Poisson
      y_test <- test_data[["pourcentage_esca"]]
      mu_test[mu_test == 0] <- 1e-6  # éviter log(0)
      dev <- 2 * sum(y_test * log(ifelse(y_test == 0, 1, y_test / mu_test)) - (y_test - mu_test))
      errors[j] <- dev
      
    }, error = function(e) {
      errors[j] <- NA
      warning(paste("Erreur pour lambda =", lambda, "dans fold", j))
    })
  }
  
  # Moyenne des erreurs (en excluant les NA)
  cv_errors[i] <- mean(errors, na.rm = TRUE)
}

# Résultat : lambda optimal
best_lambda <- lambda_grid[which.min(cv_errors)]
print(data.frame(lambda = lambda_grid, cv_deviance = cv_errors))
cat("Meilleur lambda :", best_lambda, "\n")


observations$esca_categoriel <- "[0;1["
observations[observations$pourcentage_esca >= 1,]$esca_categoriel <- "[1;2["
observations[observations$pourcentage_esca >= 2,]$esca_categoriel <- "[2;5["
observations[observations$pourcentage_esca >= 5,]$esca_categoriel <- "[5;10["
observations[observations$pourcentage_esca >= 10,]$esca_categoriel <- "[10;20["
observations[observations$pourcentage_esca >= 20,]$esca_categoriel <- "[20;100["
observations$esca_categoriel <- factor(observations$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

fixed_formula <- reformulate(to_keep[-c(1,2)], response = "esca_categoriel")
glmLasso <- glmmLasso(fix = fixed_formula, rnd=list(cepage=~1, region_viticole = ~1), 
                      data = observations[ ,c(to_keep, "esca_categoriel")], lambda = 400,
                      family=acat())

print("Accuracy train : ")
pred_train <- round(predict(glmLasso, observations[,c(to_keep, "esca_categoriel")]),0)
table(pred_train == observations$esca_categoriel)[1]/ length(pred_train)
print("Train : Table des correspondances (prédictions / réelles)")
table(pred_train, observations$esca_categoriel)

# library(glasso)
# s = var(observations[,var_tt[-c(1,2)]])
# a<-glasso(s, rho=.01)
