# V0.8.4 : Gabriel Macé
# Version modifiée de 0.8.3
# Avec moins d'indicateurs : sans le swi, et sans les années calendaires


source("R_functions/fonctions_pheno.R")
source("script/Bilan_hydrique/et0.r")
source("script/Bilan_hydrique/gv.grow.r")
source("script/Bilan_hydrique/gv.rad.int.dev.r")
source("script/Bilan_hydrique/wb.vine.1.4.r")
source("script/Bilan_hydrique/gccm.r")

get_indice <- function(per,an,can,x,parmean,parsum,parmax,parmin){
  isv_seq_faible <- ifelse((x[per[[can]],]$isv > 0.66 & x[per[[can]],]$isv <= 0.84), 0, 1)
  isv_seq_faible <- rle(isv_seq_faible)$lengths[rle(isv_seq_faible)$values == 0] # seq de jours de stress hydrique faible
  isv_seq_fai_mod <- ifelse((x[per[[can]],]$isv > 0.32 & x[per[[can]],]$isv <= 0.66), 0, 1)
  isv_seq_fai_mod <- rle(isv_seq_fai_mod)$lengths[rle(isv_seq_fai_mod)$values == 0] # seq de jours de stress hydrique fiable à modéré
  isv_seq_mod_sev <- ifelse((x[per[[can]],]$isv > 0.08 & x[per[[can]],]$isv <= 0.32), 0, 1)
  isv_seq_mod_sev <- rle(isv_seq_mod_sev)$lengths[rle(isv_seq_mod_sev)$values == 0] # seq de jours de stress hydrique modère à sévère
  isv_seq_sev <- ifelse(x[per[[can]],]$isv <=0.08, 0, 1)
  isv_seq_sev <- rle(isv_seq_sev)$lengths[rle(isv_seq_sev)$values == 0] # seq de jours de stress hydrique sévère
  
  c(an=an,
    longueur_periode = dim(x[per[[can]],])[1], # Récuperer la longueur de la période
    apply(x[per[[can]],parmean],2, mean),
    rain.days=sum(x[per[[can]],"rr"] > 1, na.rm=T),
    isv.faible.seq.5=sum(isv_seq_faible>=5), # Nbre de séq d'au moins 5 jours de stress hydrique faible d'affilés
    isv.faible.seq.10=sum(isv_seq_faible>=10),
    isv.faible.seq.15=sum(isv_seq_faible>=15),
    isv.fai_mod.seq.5=sum(isv_seq_fai_mod>=5), # Nbre de séq d'au moins 5 jours de stress hydrique faible à modéré d'affilés
    isv.fai_mod.seq.10=sum(isv_seq_fai_mod>=10),
    isv.fai_mod.seq.15=sum(isv_seq_fai_mod>=15),
    isv.mod_sev.seq.5=sum(isv_seq_mod_sev>=5), # Nbre de séq d'au moins 5 jours de stress hydrique modéré à sévère d'affilés
    isv.mod_sev.seq.10=sum(isv_seq_mod_sev>=10),
    isv.mod_sev.seq.15=sum(isv_seq_mod_sev>=15),
    isv.sev.seq.5=sum(isv_seq_sev>=5), # Nbre de séq d'au moins 5 jours de stress hydrique sévere d'affilés
    isv.sev.seq.10=sum(isv_seq_sev>=10),
    isv.sev.seq.15=sum(isv_seq_sev>=15),
    
    sum.days.isv.faible = sum(x[per[[can]],"isv"] <= 0.84 & x[per[[can]],"isv"] > 0.66, na.rm=T),
    sum.days.isv.fai_mod = sum(x[per[[can]],"isv"] <= 0.66 & x[per[[can]],"isv"] > 0.31, na.rm=T),
    sum.days.isv.mod_sev = sum(x[per[[can]],"isv"] <= 0.31 & x[per[[can]],"isv"] > 0.08, na.rm=T),
    sum.days.isv.sev = sum(x[per[[can]],"isv"] <= 0.08 , na.rm=T),
    sum.frost.days.0=sum(x[per[[can]],"tn"] < 0, na.rm=T),
    sum.heat.days.25=sum(x[per[[can]],"tx"] > 25, na.rm=T),
    sum.heat.days.30=sum(x[per[[can]],"tx"] > 30, na.rm=T),
    sum.heat.days.35=sum(x[per[[can]],"tx"] > 35, na.rm=T),
    # Aire sous la courbe de l'isv : 
    auc_isv = sum(lintegrate(x = x[per[[can]],"doy"], y = x[per[[can]],"isv"], xint = x[per[[can]],"doy"]))
  )
}

get_indice_dormance <- function(per, per_n_1, an, can,x, x_1,parmean,parsum,parmax,parmin){
  if(dim(x_1)[1]<367){ # première année
    return(get_indice(per,an,can,x,parmean,parsum,parmax,parmin))
  }
  else{
    n_1 = max(per_n_1[[can]]) # récuperer la longueur de l'année n-1
    isv_seq_faible <- c(ifelse((x[per[[can]],]$isv > 0.66 & x[per[[can]],]$isv <= 0.84), 0, 1), ifelse((x_1[per_n_1[[can]],]$isv > 0.66 & x_1[per_n_1[[can]],]$isv <= 0.84), 0, 1))
    isv_seq_faible <- rle(isv_seq_faible)$lengths[rle(isv_seq_faible)$values == 0] # seq de jours de stress hydrique faible
    isv_seq_fai_mod <- c(ifelse((x[per[[can]],]$isv > 0.32 & x[per[[can]],]$isv <= 0.66), 0, 1), ifelse((x_1[per_n_1[[can]],]$isv > 0.32 & x_1[per_n_1[[can]],]$isv <= 0.66), 0, 1))
    isv_seq_fai_mod <- rle(isv_seq_fai_mod)$lengths[rle(isv_seq_fai_mod)$values == 0] # seq de jours de stress hydrique fiable à modéré
    isv_seq_mod_sev <- c(ifelse((x[per[[can]],]$isv > 0.08 & x[per[[can]],]$isv <= 0.32), 0, 1), ifelse((x_1[per_n_1[[can]],]$isv > 0.08 & x_1[per_n_1[[can]],]$isv <= 0.32), 0, 1))
    isv_seq_mod_sev <- rle(isv_seq_mod_sev)$lengths[rle(isv_seq_mod_sev)$values == 0] # seq de jours de stress hydrique modère à sévère
    isv_seq_sev <- c(ifelse(x[per[[can]],]$isv <=0.08, 0, 1), ifelse(x_1[per_n_1[[can]],]$isv<=0.08, 0, 1))
    isv_seq_sev <- rle(isv_seq_sev)$lengths[rle(isv_seq_sev)$values == 0] # seq de jours de stress hydrique sévère
    
    return(c(an=an,
             longueur_periode = dim(x[per[[can]],])[1] + dim(x_1[per_n_1[[can]],])[1], # Récuperer la longueur de la période
             apply(x_1[c(per_n_1[[can]], per[[can]]+n_1),parmean],2, mean),
             rain.days=sum(x_1[per_n_1[[can]],"rr"] > 1, na.rm=T) + sum(x[per[[can]],"rr"] > 0, na.rm=T),
             isv.faible.seq.5=sum(isv_seq_faible>=5), # Nbre de séq d'au moins 5 jours de stress hydrique faible d'affilés
             isv.faible.seq.10=sum(isv_seq_faible>=10),
             isv.faible.seq.15=sum(isv_seq_faible>=15),
             isv.fai_mod.seq.5=sum(isv_seq_fai_mod>=5), # Nbre de séq d'au moins 5 jours de stress hydrique faible à modéré d'affilés
             isv.fai_mod.seq.10=sum(isv_seq_fai_mod>=10),
             isv.fai_mod.seq.15=sum(isv_seq_fai_mod>=15),
             isv.mod_sev.seq.5=sum(isv_seq_mod_sev>=5), # Nbre de séq d'au moins 5 jours de stress hydrique modéré à sévère d'affilés
             isv.mod_sev.seq.10=sum(isv_seq_mod_sev>=10),
             isv.mod_sev.seq.15=sum(isv_seq_mod_sev>=15),
             isv.sev.seq.5=sum(isv_seq_sev>=5), # Nbre de séq d'au moins 5 jours de stress hydrique sévere d'affilés
             isv.sev.seq.10=sum(isv_seq_sev>=10),
             isv.sev.seq.15=sum(isv_seq_sev>=15),
             
             sum.days.isv.faible = sum(x[per[[can]],"isv"] <= 0.84 & x[per[[can]],"isv"] > 0.66, na.rm=T) +
               sum(x_1[per_n_1[[can]],"isv"] <= 0.84 & x_1[per_n_1[[can]],"isv"] > 0.66, na.rm=T),
             sum.days.isv.fai_mod = sum(x[per[[can]],"isv"] <= 0.66 & x[per[[can]],"isv"] > 0.31, na.rm=T) +
               sum(x_1[per_n_1[[can]],"isv"] <= 0.66 & x_1[per_n_1[[can]],"isv"] > 0.31, na.rm=T),
             sum.days.isv.mod_sev = sum(x[per[[can]],"isv"] <= 0.31 & x[per[[can]],"isv"] > 0.08, na.rm=T) +
               sum(x_1[per_n_1[[can]],"isv"] <= 0.31 & x_1[per_n_1[[can]],"isv"] > 0.08, na.rm=T),
             sum.days.isv.sev = sum(x[per[[can]],"isv"] <= 0.08 , na.rm=T) + sum(x_1[per_n_1[[can]],"isv"] <= 0.08 , na.rm=T),
             sum.frost.days.0=sum(x[per[[can]],"tn"] < 0, na.rm=T) + sum(x_1[per_n_1[[can]],"tn"] < 0, na.rm=T),
             sum.heat.days.25=sum(x[per[[can]],"tx"] > 25, na.rm=T) + sum(x_1[per_n_1[[can]],"tx"] > 25, na.rm=T),
             sum.heat.days.30=sum(x[per[[can]],"tx"] > 30, na.rm=T) + sum(x_1[per_n_1[[can]],"tx"] > 30, na.rm=T),
             sum.heat.days.35=sum(x[per[[can]],"tx"] > 35, na.rm=T) + sum(x_1[per_n_1[[can]],"tx"] > 35, na.rm=T),
             # Aire sous la courbe de l'isv :
             auc_isv = sum(lintegrate(x = x[per[[can]],"doy"], y = x[per[[can]],"isv"], xint = x[per[[can]],"doy"])) + sum(lintegrate(x = x_1[per_n_1[[can]],"doy"], y = x_1[per_n_1[[can]],"isv"], xint = x_1[per_n_1[[can]],"doy"]))
    )
    )
  }
}


AgroEcoclim_index <- function(df, lon, lat, alt, 
                              VarNames = c(Time="Time",tn="tasmin",tx="tasmax",rr="pr", hu="hu",
                                           un=NULL, ux=NULL, rg=NULL, et0=NULL, v=NULL),
                              cepage, mat_sugar_content,
                              canopy.height,canopy.width,distance.between.rows,row.porosity,row.azimut,tx_enh = 0.35,
                              ru,fr.thres = 0, heat.thres = 35,
                              SUWE_floraison=F, # Do we use SUWE for the flowering prediction 
                              obs.pheno = data.frame(an=NA, C=NA, I=NA, V=NA, mat=NA ),
                              daily.wb.outputs=F,
                              calculate.incomplete.last.year=F,
                              i_model=NA,i_scen = NA,i_reg = NA,
                              phenologie=NULL)
{
  # 
  # 
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  # CHECK/SET CLIMATE DATA ####
  # 
  ###############################################################################################################'
  ###############################################################################################################'
  # Change column names
  colnames(df)[match(VarNames, colnames(df))] <- names(VarNames)
  
  #transform to data.frame
  df <- as.data.frame(df)
  
  # Check time serie
  if(nrow(df[which(substr(df$Time, 6, 10) == "02-30"),])>0)return(paste0("Le modèle utilisé contient des 30 février, indicateurs non calculés"))
  # Set Time and order
  df$Time <- as.Date(df$Time)
  df <- df[order(df$Time),]
  # 
  # Identifying days that are not followed by the next one....
  whichdaymiss <- which(diff(df$Time) == 2)
  if(length(whichdaymiss)>0)
  {
    df[sort(c(whichdaymiss,whichdaymiss+1)), ]
    # copy the data from the day before the missing date
    dfTail <- df[whichdaymiss,]
    # change the date so it corresponds to the missing date
    dfTail$Time <- dfTail$Time+1
    # add these dummy values (from missing dates) to the data frame
    df <- rbind(df, dfTail)
    # sort
    df <- df[order(df$Time),]
  }
  df$an <- as.numeric(format(df$Time,"%Y" ))
  df$doy <- as.numeric(format(df$Time,"%j" ))
  df$mois <- as.numeric(format(df$Time,"%m" ))
  df$date <- df$Time
  # 
  # Check that all years are complete and data is complete
  # and remove years with uncomplete data
  ans <- unique(df$an)
  NaRows <- which(apply(df[,c("date","tn","tx","rr")],1, function(x) sum(is.na(x)) > 0))
  if(length(NaRows) > 0)
  {
    NaAns <- ans[ans %in% df$an[NaRows]]
    df <- df[!(df$an %in% NaAns),]
  }
  AnsLength <- tapply(df$an, df$an, length)
  AnsUncomplete <- as.numeric(names(AnsLength[which(AnsLength < 365)]))
  if(calculate.incomplete.last.year==T) # if option calculate.incomplete.last.year 
    # is selected last year is kept even is uncomplete
    AnsUncomplete <- AnsUncomplete[ AnsUncomplete != max(df$an)]  
  
  df <- df[!(df$an %in% AnsUncomplete), ]
  # Removing potential inconsistencies (e.g. TN > TX!)
  lesquels <- which(df$tn >= df$tx)
  if(length(lesquels > 0))
  {
    txlq <- df$tn[lesquels]
    tnlq <- df$tx[lesquels]
    df$tn[lesquels] <- tnlq
    df$tx[lesquels] <- txlq
  }
  
  
  # Calcul du vpd : ####
  # Average daily temperature / Temperature moyenne quotidienne
  if(is.null(df$tm)) df$tm <- (df$tn+df$tx)/2
  psat <- 0.611 * exp((17.6 *df$tm) / (df$tm + 243) )
  df$VPD <- (psat - ((psat*df$hu)/100) )
  rm(psat)
  gc()
  
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  # PHENOLOGY ####
  # 
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  
  if(!SUWE_floraison)
  {
    # Budburst paramtrization using heat summation over 5 degrees starting on January 1st
    # García de Cortázar-Atauri, I., Brisson, N., & Gaudillere, J. P. (2009). Performance of several models for predicting budburst date of grapevine (Vitis vinifera L.). International Journal of Biometeorology, 53(4), 317‑326. https://doi.org/10.1007/s00484-009-0217-4
    # bpar <- as.list(gdd.param.values[gdd.param.values$varieties==cepage,])
    # Flowering / Veraison and Maturity 
    # Parker, et al, 2013 / Parker, A.K., et al, 2020
    Fstar <- as.list(gfv.param.values[gfv.param.values$varieties==cepage,])
    # 
    ph <- data.frame(an = unique(df$an),
                     C = c(NA, debourrement(cepage = cepage,
                                            tn = df$tn, tx = df$tx,
                                            doy = df$doy,
                                            year = df$an)$debourrement_SUWE), ## BB
                     I = floraison_bh(cepage = cepage,
                                      tn = df$tn, tx = df$tx,
                                      doy = df$doy,
                                      year = df$an,
                                      SUWE = FALSE)$floraison)# FLO) # 
    # Estim 1st year of BB
    ph$C[1] <- round(ph$I[1] - mean(ph$I[-1]-ph$C[-1]),0)
    
    
  } else {
    # Cabernet franc do not exist in SUWE, used Merlot which is the nearest 
    SUWE.param.values <- rbind(SUWE.param.values, SUWE.param.values[SUWE.param.values$varieties == "merlot",])
    SUWE.param.values$varieties[nrow(SUWE.param.values)] <- "cabernet franc"
    # 
    # Budburst paramtrization using Morales Castilla et al (PNAS, 2020) budburst model-based on Smoothed-Utah and Weng & Engel approaches
    bpar <- as.list(SUWE.param.values[SUWE.param.values$varieties==cepage,])
    # Flowering / Veraison and Maturity 
    # Parker, et al, 2013 / Parker, A.K., et al, 2020
    Fstar <- as.list(gfv.param.values[gfv.param.values$varieties==cepage,])
    # 
    ph <- data.frame(an = unique(df$an),
                     C = c(NA, debourrement(cepage = cepage,
                                            tn = df$tn, tx = df$tx,
                                            doy = df$doy,
                                            year = df$an)$debourrement_SUWE), ## BB
                     I = c(NA, floraison_bh(cepage = cepage,
                                            tn = df$tn, tx = df$tx,
                                            doy = df$doy,
                                            year = df$an,
                                            SUWE = TRUE)$floraison)# FLO
    ) # 
    # Estim 1st year of BB
    ph$C[1] <- round(ph$I[1] - mean(ph$I[-1]-ph$C[-1]),0)
    # Estim 1st year of flo
    ph$I[1] <- round(ph$V[1] - mean(ph$V[-1]-ph$I[-1]),0)
    
    
  }
  
  
  
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  # WATER BALANCE ####
  # 
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  # Potential evapotranspiration calculation / Calcul de l'ETP (PM ou Hargreaves selon dispo donnees
  if(is.null(df$et0))
  {  
    if(is.null(df$un) | is.null(df$ux) | is.null(df$rg) | is.null(df$v))
    {
      df$et0 <- etht(df$tn, df$tx, latdeg = lat, doy = df$doy)
      # Global radiation calculation if measurement is missing / Calcul du rayonnement global (Hargreaves --> Temperature)
      if(is.null(df$rg))
      {
        df$rg <- Gh.ht(df$tn, df$tx, latdeg=lat, doy=df$doy, a = 0.16)
        df$et0 <- etht(df$tn, df$tx, latdeg = lat, doy = df$doy)
      } else {
        # If Glob rad is available, then et0 is calculated with Hargreaves radiation
        df$et0 <- ethr(df$tn, df$tx,Gs = df$rg)
      }
      # Calcul de l'ET0 et remplacement de valeur manquante par la formule de Hargreaves
    } else {
      etpp <- et0(df$tn, df$tx, df$rg, df$un, df$ux, df$v,
                  latdeg = lat, alt = alt, doy=df$doy)
      # Pourcentage de valeurs manquantes :
      sum(is.na(etpp))/length(etpp)
      etpht <- etht(df$tn, df$tx, latdeg = lat, doy = df$doy)
      df$et0 <- ifelse(is.na(etpp), etpht, etpp)
      df$etht <- etpht
    }
  }
  if(is.null(df$rg)) df$rg <- Gh.ht(df$tn, df$tx, latdeg=lat, doy=df$doy, a = 0.16)
  
  dlhours <- daylength(latdeg=lat, doy=df$doy, year=df$an)
  df$durins <- sd.angstrom(Gh=df$rg, daylight.hours = dlhours, doy =  df$doy, latdeg = lat, year=df$an)
  df$fracins <- df$durins / dlhours
  
  
  # To simulate grapevine vegetative growth we use heat summation (base 10 C) from january 1st
  # Calcul des degrees jours en base 10 depuis le 1er janvier. Cela sert a simuler la croissance de la canopee
  df$st10 <- NA
  df$st10 <- as.numeric(unlist(tapply(df$tm, df$an, function(x) cumsum(pmax(x-10,0)))))
  
  # GRAPEVINE GROWTH MODELLING & WATER BALANCE
  # DATE FOR TRIMMING (FULL VEGETATION ON THE ROW) 
  # Date de rognage fictive (15 jours apres floraison) : on considere alors que la vegetation est pleinement developpee a cette date
  ph$rogn <- ph$I+15
  ans <- unique(df$an)
  # Variables of canopy development / Variable de developpement de la canopee
  pargrow <- c("poro", "row.height", "row.width", "r0")
  # Water balance related variables / Variables issues du modele de bilan hydrique
  parwb <-c("bh0","ftsw", "isv","tv","es","bh","bhv","esenh","k")
  # Output data.frame for water balance modelling / tableau de sortie du modele de bilan hydrique
  out <- df
  out[,pargrow] <- NA
  out[,parwb] <- NA
  # 
  # Loop of of water balance modelling for each year / Boucle de simulation du bilan hydrique pour chaque annee
  for(i in 1:length(ans))
  {
    x <- out[out$an == ans[i],]
    dstart <- ph[ph$an == ans[i],"C"]
    if(length(dstart) == 0) dstart = NULL
    dend <- ph[ph$an == ans[i],"rogn"]
    if(length(dend) == 0) dend = NULL
    # Grapevine canopy growth modelling / Modelisation de la croissance de la canopee
    croiss <- gv.grow(gdd = x$st10, 
                      doy.start = dstart, doy.end = dend,
                      poro.max = row.porosity, height.max = canopy.height,
                      width.max = canopy.width, k.max = 0.6, leafall.beg = as.numeric(format(as.Date("2001-11-01"),"%j")), leafall.end = as.numeric(format(as.Date("2001-11-30"),"%j")))
    #print(croiss$row.height[1])
    # Model of sunbeam interception by grapevine canopy (Riou et al 1989) / Modele d'interception du rayonnement par la vigne (Riou et al. 1989)
    kvine <- k.vine(r_lat = lat, r_longit = lon, i_year=x$an, i_doy=x$doy,
                    d_r = distance.between.rows, hf_r = croiss$row.height,
                    lf_r = croiss$row.width, az_r = row.azimut,
                    poro_r=croiss$poro, alb_f = 0.22, alb_s = 0.18,
                    rDt = 1800, r_rgi=x$rg)
    
    # Soil water balance from Lebon et al. (2003) / Bilan hydrique de Lebon et la. (2003)
    wbvine <- wb.vine(k=kvine, rr=x$rr, et0=x$et0, ru = ru,
                      ru_enh = 30, ru_sup = 5, tx_enh = tx_enh, b_rit = 0.7)
    out[out$an == ans[i],pargrow] <- as.data.frame(croiss[pargrow])
    out[out$an == ans[i],parwb] <- as.data.frame(wbvine[parwb])
  } # Water balance loop ending / fin de la boucle de calcul du bilan hydrique
  # 
  # 
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  # AGRO - ECO CLIMATIC indices ####
  # 
  ###############################################################################################################'
  ###############################################################################################################'
  # 
  #  ' Define periods ####
  # 
  # 
  flo <- ph$I ; names(flo) <- ph$an
  deb <- ph$C ; names(deb) <- ph$an
  # 
  # PERIODS ON WHICH WE CALCULATE ALL INDICES EAHC YEAR 
  rownames(ph) <- ph$an
  
  # Période de dormance
  dormance_n <- apply(ph[,c("C","I")], 1, function(x) (1:x[1]), simplify = F) # Dormance sur l'année n
  dormance_n_1 <- mapply(function(an) { # Dormance sur l'année n-1
    if (an %% 4 == 1) { # Année n-1 bissextile
      return(245:366)
    } else {
      return(244:365)
    }
  }, ph$an, SIMPLIFY = FALSE)
  names(dormance_n_1) <- ph$an 
  
  jours.an <- mapply(function(an) { # Jours par an
    if (an %% 4 == 0) { # Année n-1 bissextile
      return(1:366)
    } else {
      return(1:365)
    }
  }, ph$an, SIMPLIFY = FALSE)
  names(jours.an) <- ph$an 
  
  annee.pheno <- mapply(function(an) { # Pour récup sur l'année phénologique
    if (an %% 4 == 0) { # Année n-1 bissextile
      return(1:244)
    } else {
      return(1:243)
    }
  }, ph$an, SIMPLIFY = FALSE)
  names(annee.pheno) <- ph$an 
  
  # Débourrement à floraison
  deb.to.flo <- apply(ph[,c("C","I")], 1, function(x)  (x[1] + 1):(x[2]), simplify = F) # Débourrement à floraison 
  # Période d'observation des symptomes d'ESCA
  symptomes <- mapply(function(i, an) {
    if (an %% 4 == 0) { # Année n bissextile
      return((i + 1):244)
    } else {
      return((i + 1):243)
    }
  }, ph[,"I"], ph$an, SIMPLIFY = FALSE)
  names(symptomes) <- ph$an
  
  # Débourrement au 1er septembre
  deb.to.end <- mapply(function(i, an) {
    if (an %% 4 == 0) { # Année n bissextile
      return((i + 1):244)
    } else {
      return((i + 1):243)
    }
  }, ph[,"C"], ph$an, SIMPLIFY = FALSE)
  names(deb.to.end) <- ph$an
  
  annee <- (as.numeric(format(as.Date(paste(2001, 1, 1,sep="-")),"%j"))):
    (as.numeric(format(as.Date(paste(2001, 12, 31,sep="-")),"%j")))
  
  
  # 
  #  ' Create Synthesis tables ####
  # 
  # 
  # Synthesis tables (st)
  st.deb <-  st.flo <- c() # pheno-phase
  st.dormance <- st.deb.to.flo <- st.symptomes <- st.deb.to.end <- c() # define period around pheno-phase
  st.an.pheno <- c() # Growing seasons and whole year
  # 
  # 
  #  ' Calc Synthesis tables ####
  # 
  #  
  # **** BBOIS 2023-02-08 14:00 : Update -int?gration de extr?mes thermiques par p?riode
  #  PARMEAN / PARSUM 
  parmean <- c("tn","tx","tm","rg", "durins","fracins","ftsw","bh0","isv","tv", "bh","bhv",
               "rr",  "VPD", "hu", "et0")
  # 
  ans <- sort(unique(df$an))
  
  for(i in 1:length(ans))
  {
    an <- ans[i]
    can <- as.character(an)
    x <- as.data.frame(out[out$an == an,])
    x_1 <- as.data.frame(out[out$an %in% c(an-1, an),]) # Pour la dormance qui est sur 2 années
    st.deb   <- rbind(st.deb,get_indice(deb,an,can,x,parmean))
    st.flo   <- rbind(st.flo,get_indice(flo,an,can,x,parmean))
    st.dormance <- rbind(st.dormance,get_indice_dormance(dormance_n, dormance_n_1,an,can,x, x_1,parmean))
    st.an.pheno <- rbind(st.an.pheno,get_indice_dormance(annee.pheno, dormance_n_1,an,can,x, x_1,parmean))
    st.deb.to.flo <- rbind(st.deb.to.flo,get_indice(deb.to.flo,an,can,x,parmean))
    st.symptomes <- rbind(st.symptomes,get_indice(symptomes,an,can,x,parmean))
    st.deb.to.end <- rbind(st.deb.to.end, get_indice(deb.to.end,an,can,x,parmean))
  }
  
  #  list of st to transform in df
  st2df  <- c("st.deb", "st.flo","st.dormance", "st.deb.to.flo", "st.symptomes", "st.deb.to.end","st.an.pheno")
  for(i_st in st2df) assign(i_st, as.data.frame(get(i_st)))
  # 
  
  
  
  # 
  #  ' Aggregating all as a list ####
  # 
  #   
  outlist <- Hmisc::llist(ph, # Phenology
                          st.deb, st.flo,  #Indices at a given stage
                          st.dormance, st.deb.to.flo, st.symptomes, st.deb.to.end, # Indices during phenological periods
                          st.an.pheno)
  # If "daily.wb.outputs" is TRUE, daily data from water balance modelling is provided
  if(daily.wb.outputs==T){
    require(tidyr)
    outlist$daily.wb <- out[,c("an", "doy" ,"ftsw", "tv", "isv")]
  }
  
  return(outlist)
}
# 
# 

