# 
# Script Sébastien ZITO ; sebastien.zito@inrae.fr
#        Gabriel MACE ; gabriel.mace@inrae.fr
#        Chloe DELMAS ; chloe.delmas@inrae.fr
# 
#  avril 2025
# 
# Ce script contient nos fonctions permettant de simuler des stades phénologiques :
library(stringr) # manipulation de caractères
library(readxl)
library(tidyr)
library(dplyr)
library(data.table) # large data manipulation


### DEGREE DAYS PHENOLOGICAL MODELLING -----
dd.pheno <- function(doy, year=NULL, TN, TX, TM=NULL, t0=60, Tb=0, Fstar=NULL, give.F=F){
  if(is.null(year)) year=rep(0,length(TN))
  # Output table
  #
  years <- sort(unique(year))
  # 
  DDoutput <- setNames(vector("list", length(years)), years)
  doypredoutput <- c()
  #
  for(i in 1:length(years)){
    # i = 1
    index <- which(year == years[i])
    
    TNi <- TN[index]
    TXi <- TX[index]
    doyi <- doy[index]
    
    if(!is.null(TM)) TMi <- TM[index]
    
    if(is.null(Fstar))
    {
      cat("Fstar value not provided : using default parameters for Pinot noir blooming\n from Parker et al 2011\n")
      Fstar <-  gfv.param.values[gfv.param.values$varieties == "pinot noir", "Fstar.flo"]    
    }
    if(is.null(TM)){
      # controle TN and TX
      TNi[TNi > TXi] <- TXi[TNi > TXi]-0.1
      TMi <- (TNi+TXi)/2
    }
    if(sum(is.na(TMi)) > 0)
    {
      DD <- NA
      doypred <- NA
    } else {
      TMi[doyi < t0] <- NA
      DD <- cumsum(pmax(TMi-Tb, 0, na.rm=T))
      doypred <- which(DD < Fstar)
      
      if(length(doypred) == 0)
      {
        doypred <- NA
      } else {
        doypred <- max(doypred) + 1
      }
    }
    doypredoutput <- c(doypredoutput, doypred)
    if(give.F==T) DDoutput[[i]] <- DD
  }
  names(doypredoutput) <- years
  if(give.F==T)
  {
    return(DDoutput)
  } else {
    return(doypredoutput)
  }
}

# Wang and Engel / Garcia de Cortazar model ---------
we <- function(TM, Topt=27.4, Tmin=0, Tmax=40)
{
  alpha <- log(2) / log( (Tmax - Tmin)/(Topt - Tmin) )
  
  CTt <- ifelse(TM < Tmin | TM > Tmax, 0,
                (2 * (TM - Tmin)^alpha * (Topt-Tmin)^alpha - (TM-Tmin)^(2*alpha)) /
                  ( (Topt - Tmin)^(2*alpha) ) )
  CTt[is.na(CTt)] <- 0
  CTt <- pmax(CTt, 0)
  return(CTt)
}
# 
SUWE <- function(tn, tx, doy, year, SouthernHemisphere=F,doy.dorm.start = 213, Cc=c(chardonnay=199.01, pinot.noir=198.53)[1], Fc.bb=c(chardonnay=12.50, pinot.noir=12.98)[1],Fc.flo=c(chardonnay=22.26, pinot.noir=21.85)[1],Fc.ver=c(chardonnay=63.28, pinot.noir=59.50)[1], tm1 = -6.70, topt =  7.94, minT = -0.17, tn2 = 40.58, topt.we.bb=26.1, topt.we.flo=29.34, topt.we.ver=21.84,tmin.we=0, tmax.we=40, provide.daily.output = F, date_debourrement = NULL) {
  
  td <- (tn+tx)/2
  cu <- td
  cu[] <- NA
  years <- sort(unique(year))
  y_init <- 2
  cu[td < tm1] <-   1/(1+exp(-4*((td[td < tm1]-tm1)/(topt-tm1))))
  cu[tm1 <= td & td < topt]  <-  1 + (-0.5*(td[tm1 <= td  & td < topt]-topt)^2)/(tm1-topt)^2
  cu[topt <= td & td  < tn2]  <-  1 - ((1-minT)*(td[topt <= td  & td < tn2]-topt)^2/(2*(tn2-topt)^2))
  cu[tn2 <= td] <-  minT + ((1-minT)/(1+exp(-4*((tn2-td[tn2 <= td])/(tn2-topt)))))
  #webb <- we(td, Topt=topt.we, Tmin = tmin.we, Tmax = tmax.we)
  
  if(SouthernHemisphere == T) y_init <- 1
  output <- lapply(y_init:length(years), function(i) {
    index <- c(which(doy >= doy.dorm.start & year == years[i-1]),
               which(doy < doy.dorm.start & year == years[i]))
    if(SouthernHemisphere == T) index <- c(which(doy > doy.dorm.start & year == years[i]))
    doyi <- doy[index]
    cui <- cu[index]
    yeari <- year[index]
    # Sum of cold accumulation :
    sum.cui <- cumsum(cui)
    if(is.na(sum.cui[length(sum.cui)]))
    {
      warning(paste0("Calculation was not possible due to missing values for year :\n", years[i]))
    }
    if(sum.cui[length(sum.cui)] < Cc){
      output <- c(years[i], NA,NA,NA,NA,NA,NA,NA,NA)
      if(provide.daily.output == T) output.daily <- list(FU.BBFLO = NA,  FU.FLOVER = NA)
    } else {
      # Dormancy break index :
      dbi.index <- which(sum.cui >= Cc)[1]
      # Dormancy break day
      dbi.day <- doyi[dbi.index]
      dbi.year <- yeari[dbi.index]
      # Forcing
      sum.index <- index[dbi.index:length(index)]
      sum.thi <- we(td[sum.index], Topt=topt.we.bb, Tmin = tmin.we, Tmax = tmax.we)
      sum.thi <- cumsum(sum.thi)
      # sum.thi <- webb[sum.index]
      if(max(sum.thi) < Fc.bb){
        output <- c(years[i], dbi.day, dbi.year, NA, NA, NA,NA,NA,NA)
        if(provide.daily.output == T) output.daily <- list(FU.BBFLO = NA,  FU.FLOVER = NA)
      } else {
        sum.ph.index <- which(sum.thi >= Fc.bb)
        
        if(is.null(date_debourrement)){
          ph.index <- sum.index[sum.ph.index][1]
          # print(paste("index 1 : ", ph.index))
          bb.day <- doy[ph.index]
          bb.year <- year[ph.index] 
        }
        else{
            bb.day <- date_debourrement[date_debourrement$year==years[i],]$debourrement
            bb.year <- years[i] 
            # print(paste(bb.year, bb.day))
            if(is.na(bb.day) || identical(bb.day, numeric(0))){ 
              ph.index <- sum.index[sum.ph.index][1]}
            else ph.index <- which(doy==round(bb.day,0))[year[i] - year[1] + 2]
            # print(paste("index 2 : ", ph.index))
        }       

        sum.index <- ph.index:index[length(index)]
        we.debflo <- we(td[sum.index], Topt=topt.we.flo, Tmin = tmin.we, Tmax = tmax.we)
        sum.flo <- cumsum(we.debflo)
        sum.ph.index <- which(sum.flo >= Fc.flo)
        ph.index <- sum.index[sum.ph.index][1]
        flo.day <- doy[ph.index]
        flo.year <- year[ph.index]
        
        # print(ph.index)
        if(!is.na(ph.index))
        {
          sum.index <- ph.index:(ph.index+200) # Flowering to veraison can't last longer than 200 days (which is already quite long!)
          we.flover <- we(td[sum.index], Topt=topt.we.ver, Tmin = tmin.we, Tmax = tmax.we)
          sum.ver <- cumsum(we.flover)
          sum.ph.index <- which(sum.ver >= Fc.ver)
          #print(paste(sum.ph.index[1] , max(sum.ver), range(doy[sum.index])))
          ph.index <- sum.index[sum.ph.index][1]
          ver.day <- doy[ph.index]
          ver.year <- year[ph.index]
        } else {
          ver.day <- NA
          ver.year <- NA
          we.flover <- NA
        }  
        # bb.day <- doy[sum.index][sum.bb.index][1]
        # bb.year <- year[sum.index][sum.bb.index][1]
        #output <- rbind(output,c(years[i], dbi.day, dbi.year, bb.day, bb.year))
        if(provide.daily.output == T) output.daily <- list(FU.BBFLO = we.debflo[which(sum.flo <= Fc.flo)],  FU.FLOVER = we.flover[which(sum.ver <= Fc.ver)])
        
        output <- c(years[i], dbi.day, dbi.year, bb.day, bb.year, flo.day, flo.year, ver.day, ver.year)
      } 
    }
    if(provide.daily.output == T){
      return(list(Yearly.Data = output, Daily.Data = output.daily))
    } else {
      return(output)
    }
    
  })
  if(provide.daily.output == T)
  {
    Yearly <- lapply(output, function(x) return(x$Yearly.Data))
    Daily <- lapply(output,function(x) return(x$Daily.Data))
    Yearly <- as.data.frame(do.call("rbind",Yearly))
    names(Yearly) <- c("year","dormancy.break.doy", "dormancy.break.year","budbreak.doy", "budbreak.year", "flowering.doy", "flowering.year","veraison.doy","veraison.year" )
    names(Daily) <- Yearly$year
    return(list(Yearly=Yearly, Daily=Daily))
  } else {
    output <- as.data.frame(do.call("rbind",output))
    names(output) <- c("year","dormancy.break.doy", "dormancy.break.year","budbreak.doy", "budbreak.year", "flowering.doy", "flowering.year","veraison.doy","veraison.year" )
    return(output)
  }
}

SUWE.param.values <- data.frame(varieties=c("cabernet sauvignon",
                                            "chardonnay",
                                            "chasselas",
                                            "grenache",
                                            "merlot",
                                            "monastrell",
                                            "pinot noir",
                                            "riesling",
                                            "sauvignon blanc",
                                            "syrah",
                                            "ugni blanc",
                                            "late ripening"),
                                Cc=c(166.02,
                                     199.01,
                                     196.68,
                                     197.68,
                                     152.14,
                                     216.56,
                                     198.53,
                                     185.38,
                                     203.48,
                                     198.63,
                                     199.74,
                                     232.67),
                                Fc.bb=c(20.24,
                                        12.5,
                                        12.4,
                                        14.02,
                                        19.91,
                                        7.64,
                                        12.98,
                                        14.08,
                                        13.06,
                                        13.55,
                                        17.8,
                                        23.39
                                ),
                                Fc.flo=c(23.61,
                                         22.26,
                                         24.46,
                                         24.71,
                                         24.96,
                                         28.04,
                                         21.85,
                                         22.13,
                                         25.34,
                                         23.91,
                                         24.41,
                                         29.59
                                ),
                                Fc.ver=c(65.3,
                                         63.28,
                                         53.4,
                                         66.9,
                                         64.45,
                                         62.07,
                                         59.5,
                                         64.63,
                                         59.27,
                                         59.99,
                                         67.36,
                                         70.85
                                ))
rownames(SUWE.param.values) <- SUWE.param.values$varieties


# Les relation entre les dates de débourrement des différents cépages ont été trouvées via le site "plantgrape"

debourrement.plant.grape <- data.frame(varieties=c("meunier", "chenin", "cabernet franc", "melon",
                                                   "gewurztraminer", "trousseau", "savagnin", "semillon",
                                                   "pinot auxerrois", "gamay", "poulsard", "muscat de hambourg",
                                                   "cinsault", "carignan", "mourvedre", "muscat petits grains",
                                                   "niellucciu", "vermentinu", "fer servadou"),
                                ref = c("pinot noir", "chardonnay", "grenache", "chasselas",
                                        "chasselas", "pinot noir", "chasselas", "sauvignon blanc",
                                        "ugni blanc", "chasselas", "chasselas", "chasselas",
                                        "cabernet sauvignon", "syrah", "cabernet sauvignon", "chasselas",
                                        "chasselas", "sauvignon blanc", "grenache"),
                                ecart=c(1, 0, -1, 2,
                                        0, 0, 0, -2,
                                        0, 0, -1, 1,
                                        0, 2, 0, 0,
                                        -1, 0, 0))
rownames(debourrement.plant.grape) <- debourrement.plant.grape$varieties


# THE GFV Model -------------
##Parker, et al, 2013. Classification of varieties for their timing of flowering and veraison using a modelling approach: A case study for the grapevine species Vitis vinifera L. Agricultural and Forest Meteorology 180, 249-264. https://doi.org/10.1016/j.agrformet.2013.06.005
# Sugar : 200 g/L in Parker, A.K., et al, 2020. Temperature-based grapevine sugar ripeness modelling for a wide range of Vitis vinifera L. cultivars. Agricultural and Forest Meteorology 285-286, 107902. https://doi.org/10.1016/j.agrformet.2020.107902

# GFV and GSR parameters -----
gfv.param.values <- data.frame(
  varieties=c("cabernet franc","cabernet sauvignon", "chardonnay",
              "chasselas", "grenache", "merlot","pinot noir",
              "riesling", "sauvignon blanc", "syrah", "ugni blanc", "meunier", "aligote",
              "chenin","gewurztraminer", "trousseau", "savagnin",
              "semillon", "gamay", "poulsard", "muscat de hambourg",
              "cinsault", "carignan", "mourvedre", "muscat petits grains"
              ) ,
  Fstar.flo = c(1225,	1270,	1217,		
                1274,	1269, 1266,	1219,
                1242,	1238,	1277, 1376, 1120, 1227, 
                1280, 1230, 1179, 1155,
                1317, 1219, 1122, 1207,
                1265, 1288, 1354, 1200
                ) ,
  Fstar.ver = c(2655,	2641,	2547,		
                2342,	2750, 2627,	2511,
                2584,	2517,	2598, 2777, 2379, 2502, 
                NA, NA, NA, NA,
                NA, NA, NA, NA,
                NA, NA, NA, NA),
  Fstar.mat170 = c(2683, 2797, 2723,
                   NA,2786,2696, 2695, 
                   2693,2602,2817,3108,2398, NA, 
                   NA, NA, NA, NA,
                   NA, NA, NA, NA,
                   NA, NA, NA, NA),
  Fstar.mat180 = c(2759, 2865, 2772,
                   NA,2836,2732, 2734, 
                   3002,2671,2853,3181,2474, NA, 
                   NA, NA, NA, NA,
                   NA, NA, NA, NA,
                   NA, NA, NA, NA),
  Fstar.mat190 = c(2837, 2926, 2813,
                   NA,2890,2794, 2788, 
                   3069,2719,2934,NA,2561, NA, 
                   NA, NA, NA, NA,
                   NA, NA, NA, NA,
                   NA, NA, NA, NA),
  Fstar.mat200 = c(2909, 3031, 2892,
                   NA,2967,2856, 2838, 
                   3225,2820,2965,NA,2634, NA, 
                   NA, NA, NA, NA,
                   NA, NA, NA, NA,
                   NA, NA, NA, NA),
  Fstar.mat210 = c(2964, 3055, NA,
                   NA,3062,2904, 2899, 
                   NA,2854,2987,NA,2660, NA, 
                   NA, NA, NA, NA,
                   NA, NA, NA, NA,
                   NA, NA, NA, NA),
  Fstar.mat220 = c(3065, 3136, NA,
                   NA,3108,2962, 2933, 
                   NA,2895,3034,NA,2727, NA, 
                   NA, NA, NA, NA,
                   NA, NA, NA, NA,
                   NA, NA, NA, NA))
rownames(gfv.param.values) <- gfv.param.values$varieties



floraison.ref <- data.frame(varieties=c("melon", "pinot auxerrois"),
                                       ref = c("chardonnay", "gewurztraminer"),
                                       ecart=c(1, 0))
rownames(floraison.ref) <- floraison.ref$varieties

### Fonctions de simulations des phénologies ----

# Prédiction de la date de débourrement en fonction de :
# cepage = le cépage : String : cf : gfv et suwe params.
# tn = température minimale par jour : Numeric
# tx = température maximale par jour : Numeric
# doy = les jours juliens correspondats : Numeric
# year = les années

# Retourne : Un dataframe contenant l'année (year) 
# et le jour julien de débourrement prédit (debourrement_SUWE).

debourrement <- function(cepage, tn, tx, doy, year){
  if(cepage %in% SUWE.param.values$varieties){ # Si la date de débourrement est directement simulable 
    params = SUWE.param.values[SUWE.param.values$varieties == cepage,]
    res_SUWE = SUWE(tn = tn, tx = tx, doy = doy, year = year, Cc = params$Cc, 
                    Fc.bb = params$Fc.bb, Fc.flo = params$Fc.flo, Fc.ver = params$Fc.ver)
    return(data.frame(year = res_SUWE$budbreak.year, fin_dormance = res_SUWE$dormancy.break.doy,
                      debourrement_SUWE = res_SUWE$budbreak.doy))
  }
  else{ # Si la date de débourrement du cépage dépend de celle d'un autre
    # recuperer le cepage de reference
    ref = debourrement.plant.grape[debourrement.plant.grape$varieties == cepage, "ref"] 
    # Recuperer l'écart en jour des deux cépages
    ecart = debourrement.plant.grape[debourrement.plant.grape$varieties == cepage, "ecart"]
    res = debourrement(ref, tn, tx, doy, year) # Simuler la date de debourrement du cépage de référence
    return(data.frame(year = res$year, fin_dormance = NA, debourrement_SUWE = res$debourrement_SUWE + ecart))
  }
}


# ! Années non complètes traitées différement par gfv et SUWE !

# Prédiction de la date de floraison en fonction de :
# cepage = le cépage : String : cf : gfv et suwe params.
# tn = température minimale par jour : Numeric
# tx = température maximale par jour : Numeric
# doy = les jours juliens correspondats : Numeric
# year = les années

# Retourne : Un dataframe contenant l'année (year)
# Le jour julien de floraison prédit par gfv "flo_gfv" et SUWE "flo_SUWE" si dispo (sinon NA).
# Ajout de pouvoir donner le jour de débourrement au modèle SUWE :
floraison <- function(cepage, tn, tx, doy, year, date_debourrement = NULL){
  if(cepage %in% gfv.param.values$varieties){ # Si la date de floraison est directement simulable
    if(cepage %in% SUWE.param.values$varieties){ # Si la floraison est simulable via SUWE
      params = SUWE.param.values[SUWE.param.values$varieties == cepage,]
      res_SUWE = SUWE(tn = tn, tx = tx, doy = doy, year = year, Cc = params$Cc, 
                      Fc.bb = params$Fc.bb, Fc.flo = params$Fc.flo, Fc.ver = params$Fc.ver, 
                      date_debourrement = date_debourrement)
      params = gfv.param.values[gfv.param.values$varieties == cepage,]
      res_flo_gfv = dd.pheno(doy = doy, year = year, TN = tn, TX = tx, Fstar = params$Fstar.flo)
      res = data.frame(year = unique(year), flo_gfv = res_flo_gfv)
      res <- res %>% left_join(res_SUWE[, c("year", "flowering.doy")], by = "year") # jointure par année
      colnames(res)[3] <- "flo_SUWE"
      return(res)
    }
    else{ # Si la floraison est simulable qui via gfv
      params = gfv.param.values[gfv.param.values$varieties == cepage,]
      res_flo_gfv = dd.pheno(doy = doy, year = year, TN = tn, TX = tx, Fstar = params$Fstar.flo)
      return(data.frame(year = unique(year), flo_gfv = res_flo_gfv, flo_SUWE = NA))
    }
  }
  else{ # Si la date de floraison du cépage dépend de celle d'un autre
    # recuperer le cepage de reference
    ref = floraison.ref[floraison.ref$varieties == cepage, "ref"] 
    # Recuperer l'écart en jour des deux cépages
    ecart = floraison.ref[floraison.ref$varieties == cepage, "ecart"]
    res = floraison(ref, tn, tx, doy, year) # Simuler la date de floraison du cépage de référence
    return(data.frame(year = res$year, flo_gfv = res$flo_gfv + ecart, flo_SUWE = res$flo_SUWE + ecart))
  }
}


# Version adaptée pour le bilan hydrique :
# SUWE : booléan indiquant si on utilise le modèle SUWE (sinon : gfv)
# Retourne : Un dataframe contenant  :
#l'année et le jour julien de floraison prédit par gfv ou SUWE.
floraison_bh <- function(cepage, tn, tx, doy, year, SUWE = FALSE, date_debourrement = NULL){
  if(SUWE == TRUE){
    if(cepage %in% SUWE.param.values$varieties){ # Si la floraison est simulable via SUWE
      params = SUWE.param.values[SUWE.param.values$varieties == cepage,]
      res_SUWE = SUWE(tn = tn, tx = tx, doy = doy, year = year, Cc = params$Cc, 
                      Fc.bb = params$Fc.bb, Fc.flo = params$Fc.flo, Fc.ver = params$Fc.ver,
                      date_debourrement = date_debourrement)
      return(data.frame(year = res_SUWE$year, floraison = res_SUWE$flowering.doy))
    }
    else warning(paste0("Suwe pas dispo pour : ", cepage, " : floraison selon gfv"))
  }
  # Si SUWE n'est pas demandé ou non disponible pour ce cépage : simulation gfv :
  if(cepage %in% gfv.param.values$varieties){ # Si la date de floraison est directement simulable
    params = gfv.param.values[gfv.param.values$varieties == cepage,]
    res_flo_gfv = dd.pheno(doy = doy, year = year, TN = tn, TX = tx, Fstar = params$Fstar.flo)
    return(data.frame(year = unique(year), floraison = res_flo_gfv))
  }
  else{ # Si la date de floraison du cépage dépend de celle d'un autre
    # recuperer le cepage de reference
    ref = floraison.ref[floraison.ref$varieties == cepage, "ref"] 
    # Recuperer l'écart en jour des deux cépages
    ecart = floraison.ref[floraison.ref$varieties == cepage, "ecart"]
    res = floraison_bh(ref, tn, tx, doy, year, SUWE = FALSE) # Simuler la date de floraison du cépage de référence
    return(data.frame(year = res$year, floraison = res$floraison + ecart))
  }
}


### Proximité inter cépages des dates de floraison ----

# Les variètes que l'on ne peut pas simuler : 

# esca <- fread("data/ESCA/esca_filtre.csv") # jeu de données déjà trié
# cepage_esca <- str_to_lower(unique(esca$cepage))
# 
# # SUWE
# cepage_esca[!cepage_esca %in% SUWE.param.values$varieties]
# 
# 
# #gfv
# cepage_esca[!cepage_esca %in% gfv.param.values$varieties]
# 
# 
# # Floraison : 
# data_flo <- read_xlsx("data/DATA FLO.xlsx")[,c("année", "Id", "cepage", "DATE 50 ST", "DATE 50 JJ")]
# 
# data_flo$`DATE 50 JJ` <- as.numeric(data_flo$`DATE 50 JJ`)
# data_flo$`DATE 50 ST` <- as.numeric(data_flo$`DATE 50 ST`)
# data_flo <- drop_na(data_flo)
# 
# mean_pheno <- data_flo %>% group_by(cepage) %>% summarise(date_st = mean(`DATE 50 ST`), date_jj = mean(`DATE 50 JJ`))
# 
# # mettre en minuscule et supprimer les n et b
# data_flo$cepage <- substr(str_to_lower(data_flo$cepage),1,nchar(data_flo$cepage)-2) 
# data_flo$cepage <- str_replace(data_flo$cepage, "-", " ")
# data_flo[data_flo$cepage == "sauvignon",]$cepage <- "sauvignon blanc"
# 
# # Vérif : 
# sort(unique(data_flo$cepage)[!unique(data_flo$cepage) %in% cepage_esca])
# 
# # On ne prend que les cépages qui nous intéressent :
# 
# data_flo <- data_flo[data_flo$cepage %in% c(cepage_esca, gfv.param.values$varieties, SUWE.param.values$varieties) ,]
# unique(data_flo$cepage)
# # pas de gewurztraminer
# 
# mean_pheno <- data_flo %>% group_by(cepage) %>% summarise(date_st = mean(`DATE 50 ST`), date_jj = mean(`DATE 50 JJ`))
# 
# table(data_flo$cepage)
# 
