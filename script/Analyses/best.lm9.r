# A Function to retrieve the "best linear model", using all the possible
# combinations between all variables
# Created by B Bois, 2008-09-18

################ Modifications ##########################
# B. BOIS : 2013-07-21 --> redone totally the function to optimize the calculation
# time removing the "for" loops (used lapply and apply instead)
# Much faster now!
# B. BOIS : 2013-07-21 --> changed the outputs, and left the possibility
# to give as an output or not a lm object(if not, only coefficients are provided)
# B. BOIS : 2014-05-26 --> the models can now be selected with an additionnal criterion : Leave One Out Cross Validation (LOOCV) Effciency. The function loocv.lm was added
# B. BOIS : 2017-07-04 --> it is now possible to test also 2 order polynomial  in the linear model with the option "polynom2=T"
# B. BOIS : 2017-07-12 --> now best.lm does always LOOCV, no matter the criterion, and provides the LOOCV obs and predicted data. Note that the model coefficients change for each step of the LOOCV process, so this validation is totally independent.

########################################################


################# Functions : #########################
####### Main function
# best.lm : calculates all the possible combinations within the depending variables in a matrix or data.frame "indata". the depending (i.e. explanatory) variables are in the columns of "indata" which colnames are identified by the argument "varnames". Then all the possible models derived from these combinations are tested against the variable of "indata" to explain which colname is identified by the argment "X".
# best.lm provide the possibility to select the best model according to 3 criterions : AIC, adjusted R squared or LOOCV (leave one out crosse validation) efficiency

# Other sub functions (on which best.lm rely)
# do.lm --> calculates in a more faster manner than the lm function the essential information required to estimate an linear model (created to fasten calculations)
# loocv.lm --> calculates a leave one out cross validation performance of a given linear model
#######################################################

# Warning : the time of
# execution of this function might be very long
# it mainly depends on the number of covariates and the maximum number of
# co-variates used within the model


#### best.lm -----------------------------
best.lm <- function(indata=NULL, X=NULL, random_var = NULL, varnames=NULL, max.var=3, 
                    criterion=c("AIC","adj.r2", "LOOCV"),
                    give.model=F, kill.collinearity=F, cor.t.test.thres=0.05,
                    return.combn = F, polynom2=F,
                    comb = NULL)
{
  ### Arguments : 
  # indata : the input data.frame : y, X
  # X : the colname of indata containing the variable to explain
  # varnames : the names of the depending variables, i.e. "covariates" -> variables à tester
  # max.var : the maximum number of "covariates" that are to be used in the model
  # criterion : statistical criterion to be used to selected the best model. For LOOCV, a leave one out cross validation is performed for each model, and the model providing the highest efficiency is returned
  # give.model : in the returned list, a lm object is returned if give.model = TRUE. If FALSE, a much lighter (in memory) list is return (in order to save calculation time)
  # kill.collinearity : as using correlated "covariates" in a linear model (i.e. multicollinearity) lead to (1) (1) inaccurate model parameterization, (2) decreased statistical power, and (3) exclusion of significant predictor variables during model creation, the kill.colinearity argument set to true allows to eliminate multicollinarity by (1) comparing each pair of covariate and testing their correlation using a T test of Preason correlation coefficient, (2) remove all possible combination of models in which pairs of correlated covariates are found
    ### See the followin g paper ofr example : Graham, Michael H. 2003. ?Confronting multicollinearity in ecological multiple regression?. Ecology 84 (11): 2809-15. doi:10.1890/02-3114.
  # cor.t.test.thres : a numeric value giving the p.value equal or under whichcovariates are considered to be correlated.
  # return.combn : if TRUE no model is test, only a list provinding all the possible combination of covariates identified by "varnames"
  # polynom2 : if TRUE additionnal covariates are provided : the second order polynoms of the covariates, i.e. indata[,varnames]^2
  # comb : a list of numerical vectors identifying the combinations of covariates identified by "varnames". Typically a list return by best.lm when return.comb is set to TRUE
  
  
  if(length(criterion)>1) criterion <- criterion[1]
  icrit <- which(c("AIC","adj.r2", "LOOCV") == criterion)
  
  if(polynom2 == T)
  {
    if(return.combn == T) 
    {
      print("Combination provided are only between input 'varnames' (order 2 polynoms are excluded)")
    } else {
      varnames2 <- paste0(varnames,".2")
      var2 <- indata[, varnames]^2
      colnames(var2) <- varnames2
      varnames <- c(varnames, varnames2)
      indata <- cbind(indata,var2)
    }
  }
  dframe <- indata[,c(X,random_var,varnames)]
  dframe <- na.omit(dframe)
  m <- as.matrix(dframe)
  max.var <- min(max.var,length(varnames))
  vars <- as.list(0:max.var)
  if(is.null(comb))
    comb <- unlist(lapply(vars, function(x) combn(length(varnames), x, simplify=F)), recursive=F)[-1]
  if(kill.collinearity==T)
  {
    var.pairs <- combn(length(varnames),2)
    test.cor <- apply(var.pairs, 2, function(x) 
      cor.test(as.numeric(m[,varnames[x[1]]]), as.numeric(m[,varnames[x[2]]])))
    # avoid.comb <- var.pairs[, test.cor$p.value <= cor.t.test.thres]
    avoid.comb <- var.pairs[, test.cor$estimate >= 0.7]
    if(!is.matrix(avoid.comb)) avoid.comb <- t(t(avoid.comb))
    in.comb <- lapply(comb, function(x) apply(avoid.comb,2, function(y) sum(y %in% x) == 2))                   
    in.comb <- unlist(lapply(in.comb, function(x) sum(x) >= 1))
    comb <- comb[!in.comb]
  }
  if(return.combn == T)
  {
    return(comb)
  }
  if(criterion == "LOOCV")
  {
    lm.list <- lapply(comb, function(x) {
      loocv.lm(vnames=varnames[x], random_var = random_var, X=X, m=dframe)})    
    criterion.value <- unlist(lm.list)
    which.best <- which.max(criterion.value)
    best.loocv.efficiency <- criterion.value[which.best]
  } else {
    lm.list <- lapply(comb, function(x) {
      do.lm(vnames=varnames[x], random_var = random_var, X=X, m=dframe)})
    criterion.value <- unlist(lapply(lm.list, function(x) x[[icrit]]))
    which.best <- which.min(criterion.value*c(1,-1)[icrit])
    best.loocv.efficiency <- NA
  }
  best.vars <- varnames[comb[[which.best]]]
  Xvar <- m[,c(random_var, best.vars)]
  colnames(Xvar) <- c(random_var, best.vars) 
  best.form <- paste(X, "~", paste("(1|", random_var,")", collapse = "+"), "+", paste( best.vars, collapse="+"))
  best.mod <- lmer(formula(best.form), data = dframe)
  # lm.fit( Xvar, m[,X])
  best.coefs <- coef(best.mod)
  k <- length(best.vars) + 2
  n <- nrow(m)
  RSS <- sum(summary(best.mod)$residuals^2) # sum(best.mod$residuals^2)
  best.aic <- 2*(k+1) + n*log(RSS/n)
  best.r2 <- cor(fitted.values(best.mod), dframe[,X])^2 # cor(best.mod$fitted.values , m[,X])^2
  best.adj.r2 <- best.r2 - (1 - best.r2) * (k / (n-k-1))
  #names(best.coefs)[1] <- "(Intercept)"
  
    #  LOOCV data if LOOCV criterion was selected
  if(criterion == "LOOCV"){
  loocv.data <- loocv.lm(vnames = best.vars, random_var = random_var, X = X, m = m, return.eff = F)
  if(is.na(best.loocv.efficiency)) 
  {  best.loocv.efficiency <- 1 - (sum( (loocv.data$pred-loocv.data$obs)^2 ) / sum( (loocv.data$obs-mean(loocv.data$obs))^2 ))
  }  }
  else{loocv.data=NA}
  
  if(give.model==T)
  {
    # best.mod <- lm(formula(best.form), indata)
    return(list(model=best.mod, AIC=best.aic, Rsquared = best.r2, 
                adjRsquared=best.adj.r2, formula=best.form, coefficients=best.coefs, 
                best.loocv.efficiency=best.loocv.efficiency, loocv.data=loocv.data))
  } else {
    return(list(AIC=best.aic, Rsquared = best.r2, 
                adjRsquared=best.adj.r2, formula=best.form, coefficients = best.coefs,
                best.loocv.efficiency=best.loocv.efficiency, loocv.data=loocv.data))
  }
}

#### do.lm  -----------------------------
do.lm <- function(vnames, random_var, X, m)
{
  Xvar <- m[,c("cepage", "region_viticole", vnames)]
  colnames(Xvar) <- c("cepage", "region_viticole", vnames)
  
  form <- paste(X, "~", paste("(1|", random_var,")", collapse = "+"), "+", paste( vnames, collapse="+"))
  mod <- lmer(formula(form), data = m)
  
  # mod <- glm.fit( Xvar, m[,X])
  # mod <- lmer(pourcentage_esca~. + (1|cepage) + (1|region_viticole), data = m[,c(X, colnames(Xvar))])
  k <-length(vnames) + 2
  n <- nrow(m)
  RSS <- sum(summary(mod)$residuals^2) # sum(mod$residuals^2)
  aic <- 2*(k+1) + n*log(RSS/n)
  r2 <- cor(fitted.values(mod), m[,X])^2 # cor(mod$fitted.values , m[,X])^2
  adj.r2 <- r2 - (1 - r2) * (k / (n-k-1))
  return(c(aic, adj.r2)) 
}

#### loocv.lm -----------------------------
loocv.lm <- function(vnames, random_var, X, m, return.eff=T)
      {
        n <- as.list(1:nrow(m))
        preds <- lapply(n, function(y) 
                {
                  form <- paste(X, "~", paste("(1|", random_var,")", collapse = "+"), "+", paste( vnames, collapse="+"))
                  mod <- lmer(formula(form), data = m[-y,])
                  pred <- predict(mod, m[, c(random_var, vnames)])
                  return(pred)
                }
              )
          preds <- unlist(preds)
          Efficiency <- 1 - (sum( (preds-m[,X])^2 ) / sum( (m[,X]-mean(m[,X]))^2 ))
          if(return.eff)
          {
            return(Efficiency)
          } else {
            return(data.frame(cbind(obs=m[,X],pred=preds)))
          }
        }



######################################## Test ############################################

library(lme4)

load("data/modelisation/observations.RData")

#
# récupération des variables par type de période :
var_an_pheno <- colnames(observations)[12:47] # ne pas garder la longueur de la période = 365 ou 366
var_an <- colnames(observations)[49:84] # ne pas garder la longueur de la période
var_dormance <- colnames(observations)[85:121]
var_deb_to_flo <- setdiff(colnames(observations)[122:158] , c("sum.heat.days.35.deb_flo", "isv.faible.seq.15.deb_flo", "isv.fai_mod.seq.15.deb_flo", "isv.mod_sev.seq.10.deb_flo", "isv.mod_sev.seq.15.deb_flo"))
# Enlever sum.heat.days.35.deb_flo car que 0
var_symptomes <- setdiff(colnames(observations)[159:195], "sum.frost.days.0.symptomes") # Enlever sum.frost.days.0.symptomes car que 0
var_deb_to_end <- colnames(observations)[196:232]
var_tt <- c(var_an_pheno, var_an, var_dormance, var_deb_to_flo, var_symptomes, var_deb_to_end)

# Pas besoin de centrer réduire si on ne compare pas les coeffs de régression entre eux
# observations[, c("age_parcelle_estime","RU", "debourrement", "floraison", var_tt)] <- scale(observations[, c("age_parcelle_estime","RU", "debourrement", "floraison", var_tt)])

# RU to catégorielle ordinale
#observations$RU <- as.integer(factor(observations$RU))
observations$cepage <- as.factor(observations$cepage)
observations$region_viticole <- as.factor(observations$region_viticole)
X_names <- c("cepage", "region_viticole", "age_parcelle_estime","RU", "debourrement", "floraison", var_tt)
features <- observations[ ,c("pourcentage_esca", X_names)]
# features <- Matrix::sparse.model.matrix(~ cepage + region_viticole + . - 1, data = features, sep = "_")
# features <- as.data.frame(as.matrix(features))
b.lm = best.lm(indata = features, X = "pourcentage_esca", 
        #varnames = c(var_an_pheno,"debourrement", "floraison", "age_parcelle_estime", colnames(features)[1:35]), 
        random_var = X_names[c(1,2)],
        varnames = X_names[-c(1,2)],
        criterion = "AIC", # LOOCV, AIC, adj.r2
        max.var=3,
        kill.collinearity = TRUE
        )

b.lm.r2 = best.lm(indata = features, X = "pourcentage_esca", 
               #varnames = c(var_an_pheno,"debourrement", "floraison", "age_parcelle_estime", colnames(features)[1:35]), 
               random_var = X_names[c(1,2)],
               varnames = X_names[-c(1,2)],
               criterion = "adj.r2", # LOOCV, AIC, adj.r2
               max.var=3,
               kill.collinearity = TRUE
)

# Année phéno :
# max.var = 3 :
# AIC -> sum.heat.days.30
# adj.r2 -> sum.heat.days.30 + isv.fai_mod.seq.15"
# LOOCV -> 
     


# Tt :
# Année phéno :
# max.var = 3 :
# AIC -> 
# adj.r2 -> 
# LOOCV -> 