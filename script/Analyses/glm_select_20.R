library(lme4)
# library(MuMIn) # Pour réuperer les r² avec les glm lme4
library(dplyr)
library(ggplot2)
library(Metrics) # RMSE
library(openxlsx)

load("data/modelisation/observations_20.RData")

# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[12:47] # ne pas garder la longueur de la période = 365 ou 366
var_an <- colnames(observations)[49:84] # ne pas garder la longueur de la période
var_dormance <- colnames(observations)[85:121]
var_deb_to_flo <- setdiff(colnames(observations)[122:158] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[159:195], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[196:232]
var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)

observations$cepage <- as.factor(observations$cepage)
observations$region_viticole <- as.factor(observations$region_viticole)
i_train = sample(1:dim(observations)[1], replace = FALSE, size = 4400)

observations_20 <- observations[observations$pourcentage_esca < 20, ]


######################################################################

var_tt <- c("cepage", "region_viticole", "age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
r2s <- c() # adj.r2 <- r2 - (1 - r2) * (k / (n-k-1))
vars <- c()
rmses <- c()
for(var in var_tt){
  vars <- c(vars, var)
  lm <- lm(formula = paste0("pourcentage_esca~",var), data = observations_20[, c("pourcentage_esca", var)])
  r2s <- c(r2s, summary(lm)$r.squared)
  rmses <- c(rmses, rmse(observations_20$pourcentage_esca, predict(lm, var=observations_20[,var])))
}
res <- data.frame(variable = vars, r2 = r2s, rmse = rmses)
write.xlsx(res[order(res$rmse, decreasing = FALSE),] , "data/resultats/r2_rmse_par_variable_20.xlsx", rowNames=FALSE)

#####################################################################

var_tt <- c("age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
observations_20[,var_tt] <- scale(observations_20[,var_tt]) # centrer réduire


# Prendre initialement que le cépage et la région viticole et ajouter une à une de manière itérative, 
# la variable améliorant le plus le r2 du modèle.
var2keep <- c("cepage", "region_viticole")
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(var2keep, "pourcentage_esca")])

r2s <- c(cor(predict(lm_select, observations_20[,var2keep]), observations_20$pourcentage_esca)^2) # 0.24
epsilone <- 0.001

for(i in 1:50){ # garder 50 variables max en plus de cépages et régions viticoles
  r2_var <- c()
  for(var in setdiff(var_tt, var2keep)){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(var2keep, var, "pourcentage_esca")])
    r2_var <- c(r2_var, cor(predict(lm_select, observations_20[,c(var2keep,var)]), observations_20$pourcentage_esca)^2)
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
write.xlsx(res, "data/resultats/res_select_glm_r2_20.xlsx", rowNames=FALSE)



# Pareil avec rmse :
var2keep <- c("cepage", "region_viticole")
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(var2keep, "pourcentage_esca")])

rmses <- c(Metrics::rmse(predict(lm_select, observations_20[,var2keep]), observations_20$pourcentage_esca))
epsilone <- 0.001

for(i in 1:50){ # garder 50 variables max en plus de cépages et régions viticoles
  rmse_var <- c()
  for(var in setdiff(var_tt, var2keep)){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(var2keep, var, "pourcentage_esca")])
    rmse_var <- c(rmse_var, Metrics::rmse(predict(lm_select, observations_20[,c(var2keep,var)]), observations_20$pourcentage_esca))
    
  }
  if(min(rmse_var)[1] < rmses[i] - epsilone){
    rmses <- c(rmses, min(rmse_var)[1])
    var2keep <- c(var2keep, setdiff(var_tt, var2keep)[which(rmse_var == min(rmse_var)[1])[1]])
  }
  else{break}
  print(paste(i, ":", var2keep[i+2]))
}
plot(rmses, type = "l")
res <- data.frame(variable = var2keep, rmse = c(rmses[1], rmses))
write.xlsx(res, "data/resultats/res_select_glm_rmse_20.xlsx", rowNames=FALSE)

# "et0.symptomes", "ftsw.dormance", "sum.days.isv.sev.symptomes"
# sont sélectionnés par r2 et rmse mais beaucoup moins importants que région viticole et cépage

## Les variables qui handicapent le modèle :

var2rm <- c()
var_tt <- c("cepage", "region_viticole", "age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(var_tt, "pourcentage_esca")])
rmse1 <- c(Metrics::rmse(predict(lm_select, observations_20[,var_tt]), observations_20$pourcentage_esca))
rmses <- c()

for(var in setdiff(var_tt, c("cepage", "region_viticole"))){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(setdiff(var_tt, var), "pourcentage_esca")])
    rmse_var <- Metrics::rmse(predict(lm_select, observations_20[,setdiff(var_tt, var)]), observations_20$pourcentage_esca)
    if(rmse_var < rmse1){
    rmses <- c(rmses, rmse_var)
    var2rm <- c(var2rm, var)
    }
}
res <- data.frame(variable = var2rm, rmse = rmses)


## Ou faire l'inverse : partir du lm totale et enlevé ce qui ne change pas grand chose au modèle

var2rm <- c()
var_tt <- c("cepage", "region_viticole", "age_parcelle_estime", "RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)
lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(var_tt, "pourcentage_esca")])
rmses <- c(Metrics::rmse(predict(lm_select, observations_20[,var_tt]), observations_20$pourcentage_esca))

for(i in 1:(length(var_tt)-2)){
  rmse_vars <- c()
  var2test = setdiff(var_tt, c(var2rm, "cepage", "region_viticole"))
  for(var in var2test){
    lm_select <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = observations_20[,c(setdiff(var_tt, c(var, var2rm)), "pourcentage_esca")])
    rmse_vars <- c(rmse_vars,Metrics::rmse(predict(lm_select, observations_20[,setdiff(var_tt, c(var, var2rm))]), observations_20$pourcentage_esca))
  }
  rmses <- c(rmses, max(rmse_vars)[1])
  var2rm <- c(var2rm, var2test[which(rmse_vars == max(rmse_vars)[1])[1]])
}
res <- data.frame(variable = c("cepage", "region_viticole", var2rm), rmse = c(rmses[1], rmses[1], rmses[2:219]))
write.xlsx(res, "data/resultats/res_select_inv_glm_rmse_20.xlsx", rowNames=FALSE)




########################################################################################################
##### Modèles pénalisés :

library(glmnet)

features <- observations_20[ ,var_tt]
features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")

ridge <- glmnet(x = features, y = observations_20$pourcentage_esca, alpha = 0)
Metrics::rmse(predict(ridge, features), observations_20$pourcentage_esca)

cv_ridge <- cv.glmnet(features, observations_20$pourcentage_esca, alpha = 0)
ridge <- glmnet(features, observations_20$pourcentage_esca, alpha = 0, lambda = cv_ridge$lambda.min)
Metrics::rmse(predict(ridge, features), observations_20$pourcentage_esca)

lasso <- glmnet(x = as.matrix(observations_20[ ,-c(1, 2, 5, 6)]), y = observations_20$pourcentage_esca, alpha = 1 )
Metrics::rmse(predict(lasso, as.matrix(observations_20[ ,-c(1, 2, 5, 6)])), observations_20$pourcentage_esca)

lasso <- glmnet(x = features, y = observations_20$pourcentage_esca, alpha = 1) 
Metrics::rmse(predict(lasso, features), observations_20$pourcentage_esca)

# , family = poisson(link = "log") = pas top

cv_lasso <- cv.glmnet(features, observations_20$pourcentage_esca, alpha = 1)
lasso <- glmnet(features, observations_20$pourcentage_esca, alpha = 1, lambda = cv_lasso$lambda.min)
Metrics::rmse(predict(lasso, features), observations_20$pourcentage_esca)
# 5.017208
res_lasso = data.frame(var = var_tt[-c(1,2)], beta = abs(lasso$beta[36:253]))


score_var <- read.xlsx("data/resultats/r2_rmse_par_variable.xlsx")
vars <- var_tt[-c(1,2)]
to_rm <- c()
for(var1 in vars){
  if(!var1 %in% to_rm){
    for(var2 in setdiff(vars, c(var1, to_rm))){
      if(cor(observations_20[,var1], observations_20[,var2]) > 0.8){
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
                      data = observations_20[ ,c(to_keep, "pourcentage_esca")], lambda = 400,
                      family = poisson(link = "log")) #family = gaussian(link = "identity"))
Metrics::rmse(predict(glmLasso, observations_20[ ,c(to_keep, "pourcentage_esca")]), observations_20$pourcentage_esca)
# 6.153954


print("Accuracy train : ")
pred_train <- round(predict(glmLasso, observations_20[,c(to_keep, "esca_categoriel")]),0)
table(pred_train == observations_20$esca_categoriel)[1]/ length(pred_train)
print("Train : Table des correspondances (prédictions / réelles)")
table(pred_train, observations_20$esca_categoriel)



# Grille de lambdas à tester
lambda_grid <- c(200, 300, 400, 500) # 1, 5, 10, 20, 50, 100, 

# Nombre de folds
k <- 5
folds <- sample(rep(1:k, length.out = nrow(observations_20)))


# Stocker les erreurs moyennes
cv_errors <- numeric(length(lambda_grid))

for (i in seq_along(lambda_grid)) {
  lambda <- lambda_grid[i]
  errors <- numeric(k)
  
  for (j in 1:k) {
    train_data <- observations_20[folds != j, ]
    test_data <- observations_20[folds == j, ]
    
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


observations_20$esca_categoriel <- "[0;1["
observations_20[observations_20$pourcentage_esca >= 1,]$esca_categoriel <- "[1;2["
observations_20[observations_20$pourcentage_esca >= 2,]$esca_categoriel <- "[2;5["
observations_20[observations_20$pourcentage_esca >= 5,]$esca_categoriel <- "[5;10["
observations_20[observations_20$pourcentage_esca >= 10,]$esca_categoriel <- "[10;20["
observations_20[observations_20$pourcentage_esca >= 20,]$esca_categoriel <- "[20;100["
observations_20$esca_categoriel <- factor(observations_20$esca_categoriel, levels = c("[0;1[", "[1;2[", "[2;5[", "[5;10[", "[10;20[", "[20;100["))

fixed_formula <- reformulate(to_keep[-c(1,2)], response = "esca_categoriel")
glmLasso <- glmmLasso(fix = fixed_formula, rnd=list(cepage=~1, region_viticole = ~1), 
                      data = observations_20[ ,c(to_keep, "esca_categoriel")], lambda = 400,
                      family=acat())

print("Accuracy train : ")
pred_train <- round(predict(glmLasso, observations_20[,c(to_keep, "esca_categoriel")]),0)
table(pred_train == observations_20$esca_categoriel)[1]/ length(pred_train)
print("Train : Table des correspondances (prédictions / réelles)")
table(pred_train, observations_20$esca_categoriel)

# library(glasso)
# s = var(observations_20[,var_tt[-c(1,2)]])
# a<-glasso(s, rho=.01)
