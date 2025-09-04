library(caret)
library(dplyr)
library(ggplot2)
library(Metrics)

load("data/modelisation/observations.RData")

# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[11:47]
var_an <- colnames(observations)[48:84]
var_dormance <- colnames(observations)[85:121]
var_deb_to_flo <- setdiff(colnames(observations)[122:158] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[159:195], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[196:232]
var_tt <- c("RU", "debourrement", "floraison", var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)


# RU to catégorielle ordinale
observations$RU <- as.integer(factor(observations$RU))
features <- observations[ ,-c(1, 2, 5, 6)]
features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")
to_train <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = observations[ ,-c(1, 2, 6)], sep = "_")
i_train = sample(1:dim(observations)[1], replace = FALSE, size = 3500)
#train = observations[i_train ,c("cepage", "region_viticole", "age_parcelle_estime","RU", "debourrement", "floraison", var_an_pheno)]
train = features[i_train ,]


#########################################
library(sgd)
observations_enc <- model.matrix(~ cepage + region_viticole + . - 1, data = observations[, -c(1, 2, 5, 6)])
y <- observations[, "pourcentage_esca"]
sgd_lm <- sgd(x = observations_enc[i_train,], y = y[i_train], model = "lm", 
              model.control = list(family =gaussian()),
              sgd.control = list(lr = 0.01))#poisson()


y_pred <- predict(sgd_lm, newdata = observations_enc[-i_train,], type = "response")
rmse(y[-i_train], y_pred)
y_pred <- predict(sgd_lm, newdata = observations_enc[i_train,], type = "response")
rmse(y[i_train], y_pred)

sgd_lm <- sgd(x = observations_enc, y = y, model = "lm",
              model.control = list(family = gaussian()))
y_pred <- predict(sgd_lm, newdata = observations_enc, type = "response")
rmse(y, y_pred)



############################################################
# GLM optimisé par descente de gradient via pytorch

library(torch)

# Effets aléatoires : 
X <- observations[, -c(1,2,5,6)] %>%
  mutate(cepage = as.factor(cepage),
    region_viticole = as.factor(region_viticole))

# Donne à chaque modalité une valeur entière
cepage_id <- as.integer(X$cepage)
region_id <- as.integer(X$region_viticole)

# Pour les effets aléatoires :
X <- model.matrix(~.,X)[,-1]

# Définition de y :
y <- observations$pourcentage_esca


# Conversion en tenseurs
y_tensor <- torch_tensor(as.numeric(y), dtype = torch_float())
X_tensor <- torch_tensor(as.matrix(X), dtype = torch_float())
cepage_tensor <- torch_tensor(cepage_id, dtype = torch_long())
region_tensor <- torch_tensor(region_id, dtype = torch_long())

# Définition du modèle :

model_mixed <- nn_module(
  initialize = function(input_dim, n_cepages, n_regions) {
    self$fixed = nn_linear(input_dim, 1)
    self$cepage_effect = nn_embedding(n_cepages + 1, 1)
    self$region_effect = nn_embedding(n_regions + 1, 1)
  },
  forward = function(x, cepage, region) {
    fixed_part = self$fixed(x)
    cepage_part = self$cepage_effect(cepage)$squeeze(2)
    region_part = self$region_effect(region)$squeeze(2)
    output = fixed_part$squeeze(2) + cepage_part + region_part
    return(output)
  }
)


# Entrainement du modèle :
net <- model_mixed(input_dim = ncol(X), n_cepages = max(cepage_id), n_regions = max(region_id))
optimizer <- optim_adam(net$parameters, lr = 0.0001)

num_epochs <- 100000
for (epoch in 1:num_epochs) {
  optimizer$zero_grad()
  output <- net(X_tensor, cepage_tensor, region_tensor)
  loss <- nnf_mse_loss(output, y_tensor)
  loss$backward()
  optimizer$step()
  if (epoch %% 10 == 0) {
    cat("Epoch:", epoch, "Loss:", loss$item(), "\n")
  }
}

y_pred <- net(X_tensor, cepage_tensor, region_tensor)
rmse(observations$pourcentage_esca, y_pred) # 5.05241
mean(abs(observations$pourcentage_esca - y_pred)) # 3.39331
