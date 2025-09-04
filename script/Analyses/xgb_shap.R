library(xgboost)
library(caret) # cross validation
library(dplyr)
library(ggplot2)
library(Metrics)
source("R_functions/shap.R")

load("data/modelisation/observations.RData")


# RU to catégorielle ordinale
observations$RU <- as.integer(factor(observations$RU))
observations$cepage <- as.factor(observations$cepage)
observations$region_viticole <- as.factor(observations$region_viticole)
features <- observations[ ,-c(1, 2, 5, 6)]
features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")
i_train = sample(1:dim(observations)[1], replace = FALSE, size = 3500)
#train = observations[i_train ,c("cepage", "region_viticole", "age_parcelle_estime","RU", "debourrement", "floraison", var_an_pheno)]
train = features[i_train ,]

var_an_pheno <- colnames(observations)[12:47] # ne pas garder la longueur de la période = 365 ou 366
var_an <- colnames(observations)[49:84] # ne pas garder la longueur de la période
var_dormance <- colnames(observations)[85:121]
var_deb_to_flo <- setdiff(colnames(observations)[122:158] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[159:195], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[196:232]
var_tt <- c(var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)


# XGboost : Forêt aléatoire optimisé par descente de gradient :  


# Uniquement sur cépage et région viticole :
model_xgb = xgboost(data = as.matrix(features[, 1:35]), 
                    nround = 1000, 
                    objective= "count:poisson", # "reg:squarederror" 
                    label= observations[, "pourcentage_esca"],
                    eta = 0.1,
                    max_depth = 100)

plot(model_xgb$evaluation_log$train_poisson_nloglik)
rmse(observations$pourcentage_esca, predict(model_xgb, as.matrix(features[,1:35])))
# 5.329678
mean(abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features[,1:35]))))
# 3.47226
plot(observations$pourcentage_esca, abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features[,1:35]))))
plot(density(abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features[,1:35])))))




# Sur toutes les variables :
# Le max du max :
i_train = sample(1:dim(observations)[1], replace = FALSE, size = 3500)
model_xgb = xgboost(data = as.matrix(features[i_train,]), nround = 1000, 
                     objective= "count:poisson", # "reg:squarederror" 
                     label= observations[i_train, "pourcentage_esca"],
                    eta = 0.1,
                    max_depth = 100)  

plot(model_xgb$evaluation_log$train_poisson_nloglik)
rmse(observations$pourcentage_esca[i_train], predict(model_xgb, as.matrix(features[i_train,])))
mean(abs(observations$pourcentage_esca[i_train] - predict(model_xgb, as.matrix(features[i_train,]))))

rmse(observations$pourcentage_esca[-i_train], predict(model_xgb, as.matrix(features[-i_train,])))
mean(abs(observations$pourcentage_esca[-i_train] - predict(model_xgb, as.matrix(features[-i_train,]))))

plot(observations$pourcentage_esca, abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features))),
     xlab = "Incidence d'esca", ylab = "Erreur absolue", 
     main = "Erreur selon l'incidence d'esca, XGBoost max du max")
plot(density(abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features)))))


# plot(model_xgb$evaluation_log$train_poisson_nloglik)
# rmse(observations$pourcentage_esca[-i_train], predict(model_xgb, as.matrix(features[-i_train ,])))
# mean(abs(observations$pourcentage_esca[-i_train] - predict(model_xgb, as.matrix(features[-i_train ,]))))
# plot(observations$pourcentage_esca[-i_train], abs(observations$pourcentage_esca[-i_train] - predict(model_xgb, as.matrix(features[-i_train ,]))))
# plot(density(abs(observations$pourcentage_esca[-i_train] - predict(model_xgb, as.matrix(features[-i_train ,])))))
# 
# rmse(observations$pourcentage_esca[i_train], predict(model_xgb, as.matrix(features[i_train ,])))
# mean(abs(observations$pourcentage_esca[i_train] - predict(model_xgb, as.matrix(features[i_train ,]))))

# Ces paramètres engendre un fort sur-apprentissage ! (surtout le max_depth=100)
# -> il ne faudra pas les utiliser pour de la prédiction

# XGBoost cross-validé :
# Définition de la grille d’hyperparamètres
grid <- expand.grid(
  nrounds = c(100, 200, 500),
  max_depth = c(3, 6, 9),
  eta = c(0.001, 0.01, 0.05, 0.1, 0.3),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

# Cross-validation
control <- trainControl(
  method = "cv",
  number = 10,
  verboseIter = TRUE
)

# Entraînement
model_xgb_cv <- train(
  x = features,
  y = observations[, "pourcentage_esca"],
  method = "xgbTree",
  trControl = control,
  tuneGrid = grid
)

# Meilleur modèle (sans sur-apprentissage !) == 100 rounds, profondeur max = 6, eta = 0.1

# model_xgb_cv = xgb.cv(data = as.matrix(train), label= observations[i_train, "pourcentage_esca"],
#                       nfold = 10, # 10 sous jeu de données
#                    nround = 200, 
#                    objective = "count:poisson",
#                    eta = 0.05,
#                    max_depth = 6,
#                    metrics = "rmse",
#                    callbacks =  list(cb.cv.predict(save_models = TRUE)) ) 
# plot(model_xgb_cv$evaluation_log$test_rmse_mean, type = "l")
# rmse(observations$pourcentage_esca[i_train], model_xgb_cv$pred)
# mean(abs(observations$pourcentage_esca[i_train] - model_xgb_cv$pred))
# 
# rmse(observations$pourcentage_esca[i_train], predict(model_xgb_cv, as.matrix(features[i_train ,])))
# mean(abs(observations$pourcentage_esca[i_train] - predict(model_xgb_cv$models[[1]], as.matrix(features[i_train ,]))))


nrs <- seq(20, 100, 10)
rmses <- c()
for(nr in nrs){
  model_xgb = xgboost(data = as.matrix(features[i_train,]), nround = nr, 
                      objective= "count:poisson", # "reg:squarederror" 
                      label= observations[i_train, "pourcentage_esca"],
                      eta = 0.1,
                      max_depth = 6)  
  rmses <- c(rmses, rmse(observations$pourcentage_esca[-i_train], predict(model_xgb, as.matrix(features[-i_train,]))))
}

plot(nrs, rmses)

nrs[which(rmses == min(rmses))] # 4.724608

model_xgb = xgboost(data = as.matrix(features[i_train,]), nround = 80, 
                    objective= "count:poisson", # "reg:squarederror" 
                    label= observations[i_train, "pourcentage_esca"],
                    eta = 0.1,
                    max_depth = 6)  

plot(model_xgb$evaluation_log$train_poisson_nloglik, type = "l")

rmse(observations$pourcentage_esca[i_train], predict(model_xgb, as.matrix(features[i_train,])))
mean(abs(observations$pourcentage_esca[i_train] - predict(model_xgb, as.matrix(features[i_train,]))))

rmse(observations$pourcentage_esca[-i_train], predict(model_xgb, as.matrix(features[-i_train,])))
mean(abs(observations$pourcentage_esca[-i_train] - predict(model_xgb, as.matrix(features[-i_train,]))))



# plot(observations$pourcentage_esca, abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features))))
# plot(density(abs(observations$pourcentage_esca - predict(model_xgb, as.matrix(features)))))
# plot(observations$pourcentage_esca, predict(model_xgb, as.matrix(features)))

shap = shap.score.rank(xgb_model = model_xgb, 
                                   X_train = as.matrix(features),
                                   shap_approx = F)
var_importance(shap, top_n=10)

## Prepare data for top N variables
shap_long = shap.prep(shap = shap,
                           X_train = as.matrix(features) , 
                           top_n = 10
                           )

## Plot shap overall metrics
plot.shap.summary(data_long = shap_long)


xgb.plot.shap(data = features, # input data
              model = model_xgb, # xgboost model
              features = names(shap$mean_shap_score[1:10]), # only top 10 var
              n_col = 3, # layout option
              plot_loess = T # add red line to plot
)

shap_ind <- cbind(shap_long[1:5489,]$value)

for(i in 2:10){
  shap_ind <- cbind(shap_ind, shap_long[(5489*(i-1)+1) : (5489*i),]$value)
}
shap_ind <- as.data.frame(shap_ind)
colnames(shap_ind) <- names(shap$mean_shap_score[1:10])
# ou direct shap$score




##############################################################"
# Avec sélection de variable :
# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[11:47]
var_an <- colnames(observations)[48:84]
var_dormance <- colnames(observations)[85:121]
var_deb_to_flo <- setdiff(colnames(observations)[122:158] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[159:195], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[196:232]
var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)


enlever <- c("bh0", "bh", "tn", "tx", "rr", "et0", "sum.heat.days.25", "sum.heat.days.30", "hu", "isv", "longueur_periode",
             "isv.faible.seq.5", "isv.faible.seq.10", "isv.faible.seq.15",
             "isv.fai_mod.seq.5", "isv.fai_mod.seq.10", "isv.fai_mod.seq.15",
             "isv.mod_sev.seq.5", "isv.mod_sev.seq.10", "isv.mod_sev.seq.15",
             "isv.sev.seq.5", "isv.sev.seq.10", "isv.sev.seq.15", "auc_isv") # sum.days.isv.sev"
# Enlever les phéno et garder que tm ?
to_use <- setdiff(c(var_an_pheno,"debourrement", "floraison", "annee", "cepage", "age_parcelle_estime", "region_viticole"), enlever) 

features <- observations[ ,to_use]
features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")


model_xgb = xgboost(data = as.matrix(features[i_train,]), nround = 80, 
                    objective= "count:poisson", # "reg:squarederror" 
                    label= observations[i_train, "pourcentage_esca"],
                    eta = 0.1,
                    max_depth = 6)  

plot(model_xgb$evaluation_log$train_poisson_nloglik, type = "l")

rmse(observations$pourcentage_esca[i_train], predict(model_xgb, as.matrix(features[i_train, ])))
mean(abs(observations$pourcentage_esca[i_train] - predict(model_xgb, as.matrix(features[i_train,]))))

rmse(observations$pourcentage_esca[-i_train], predict(model_xgb, as.matrix(features[-i_train,])))
mean(abs(observations$pourcentage_esca[-i_train] - predict(model_xgb, as.matrix(features[-i_train,]))))
