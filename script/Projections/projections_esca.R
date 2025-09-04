#
### Gabriel Macé - 07/2025
##
# Le but de ce script est de mettre en place une démarche permettant de simuler l'incidence d'esca en fonction
# de données météorologiques simulées.

library(terra) # données spatiales
library(dplyr)
library(parallel)


# Récupération des mailles ----
mailles <- read.csv("data/Projections/Qmap_coord_used.csv")
# Les lignes de ce tableau correspondantes aux colonnes des tableaux de projections climatiques

### Choisir NOS mailles et ajouter les cépages
maille_cepage <- read.csv2("data/geographie/correspondance_maille_cepage.csv")  

# Transformation en données géographiques pour pouvoir transformer les coordonnées
# et les joindre avec nos mailles d'intérêt
coordinates <- vect(mailles,geom = c("lon", "lat"), crs = "EPSG:4326")
coordinates <- terra::project(coordinates, "EPSG:27572")
coordinates <- as.data.frame(geom(coordinates)[,c(3,4)])
colnames(coordinates) <- c("X_maille", "Y_maille")
coordinates$X_maille <- round(coordinates$X_maille / 100, 0)
coordinates$Y_maille <- round(coordinates$Y_maille / 100, 0)
coordinates$index_maille <- as.numeric(rownames(coordinates))

coordinates <- coordinates %>% 
  left_join(maille_cepage, by = c("X_maille", "Y_maille"), relationship = "one-to-many") %>%
  na.exclude()

mailles2keep <- unique(coordinates$index_maille) 



########## Pour bilan hydrique :

source("script/Bilan_hydrique/AE.clim.ind_V0.8.4_Gabriel.r") # Fonctions bilan hydrique

# load("data/climat/climat_maille_cepage_jour_1995_2014.RData") # Données climatiques par maille / cépage
description_parcelles <- read.csv("data/Reg_viti_water_balance_param.csv") # Pratiques agricoles régionales
description_parcelles[description_parcelles$cepage == "cabernet-sauvignon",]$cepage <- "cabernet sauvignon"

# Récup les RU et les altitues des mailles
info_mailles <- read.csv("data/geographie/correspondance_maille_commune.csv", sep=",")
info_mailles <- info_mailles %>% 
  select(X_maille, Y_maille, ID_MAILLE, RU, altitude, region_viticole) %>% unique()


############

source("R_functions/functions_projections.r") # Fonction projection incidence esca

# CMIP6 models list

# cmip6_models <- c("MPI-ESM1-2-HR") # ACCESS-ESM1-5 # Pour tester avec un seul modèle climatique

### Liste totale des modèles :
# cmip6_models <- c("ACCESS-CM2","ACCESS-ESM1-5","BCC-CSM2-MR","CanESM5","CMCC-ESM2",
#                   "EC-Earth3","EC-Earth3-CC","EC-Earth3-Veg","EC-Earth3-Veg-LR","FGOALS-g3",
#                   "GFDL-ESM4","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G",
#                   "KIOST-ESM","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR", "NESM3",
#                   "NorESM2-LM","NorESM2-MM","TaiESM1")
# "KACE-1-0-G" -> 360 jours

cmip6_models <- c("ACCESS-CM2", "CanESM5", "CMCC-ESM2", "EC-Earth3-CC",
                  "EC-Earth3-Veg-LR", "EC-Earth3-Veg", "EC-Earth3", "FGOALS-g3",
                  "GFDL-ESM4", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
                  "KIOST-ESM", "MIROC6", "MPI-ESM1-2-HR",
                  "MPI-ESM1-2-LR", "NorESM2-LM", "NorESM2-MM", "TaiESM1")

# Pb : "KIOST-ESM" : hu du scénario ssp245 : 365 jours piles par an


repository <- "D:/" # "data/Projections/Climat/"
cmip6_files <- list.files(repository)
# ne pas considérer les fichiers inutiles du disque dur
cmip6_files <- setdiff(cmip6_files, c("System Volume Information", "$RECYCLE.BIN")) 

# Tout simuler et récup que les incidences :
for(i_model in cmip6_models){
  print(paste0("*** J'en suis au modèle : ", i_model, " ***"))
    for(i_scen in c("ssp585", "ssp245")){ # "historical", 
      if(i_scen == "ssp245" & i_model == "KIOST-ESM") break
      print(paste0("Et au scénario : ", i_scen))
  
        my_period <- if(i_scen ==  "historical") "1995-2014" else c("2081-2100", "2046-2065")
        
        #  Check if pr tasmax and tasmin exist
      if(paste0("tasmax_",i_scen,"_",i_model, ".RData") %in% cmip6_files & 
         paste0("tasmin_",i_scen,"_",i_model, ".RData") %in% cmip6_files &
         paste0("pr_",i_scen,"_",i_model, ".RData") %in% cmip6_files){
          # load data
        tm <- get(load(paste0(repository, "tasmin_",i_scen,"_",i_model, ".RData")))[,mailles2keep]
        date <- format(as.Date(rownames(tm)), "%Y%m%d")
        for(period in my_period){
        if(i_scen == "historical") idx_date = which(date >= 19940000 & date < 20150000) 
        else if(period == "2046-2065") idx_date = which(date >= 20450000 & date < 20660000) # de 2046 à 2065
        else idx_date = which(date >= 20800000 & date < 21010000) # de 2081 à 2100
        date <- date[idx_date]
        tm <- tm[idx_date,]
            
        tx <- get(load(paste0(repository, "tasmax_",i_scen,"_",i_model, ".RData")))[idx_date,mailles2keep]
        pr <- get(load(paste0(repository, "pr_",i_scen,"_",i_model, ".RData")))[idx_date,mailles2keep]
        hu <- get(load(paste0(repository, "hurs_",i_scen,"_",i_model, ".RData")))[idx_date,mailles2keep]
        rm(myqmap_data); gc()
        
        ### Joindre toutes les projections climatiques
        meteo <- data.frame(
          date = rep(date, times = length(mailles2keep)),
          idx_maille = rep(mailles2keep, each = length(date)),
          tm = as.vector(tm) - 273.15,
          tx = as.vector(tx) - 273.15,
          pr = as.vector(pr),
          hu = as.vector(hu)
        )
        
        
        # Si la pluie est négative la mettre à 0 :
        meteo$pr <- pmax(meteo$pr, 0)
        # Si l'humidité est sup à 100 la mettre à 100
        meteo$hu <- pmin(meteo$hu, 100)
        # Si la température max est inférieure à la température min : la mettre à la températue min
        meteo$tx <- pmax(meteo$tx, meteo$tm)
        # Si la température min est supérieure à la température min : la mettre à la températue max
        meteo$tm <- pmin(meteo$tx, meteo$tm)
        
        ### Joindre les données climatiques aux mailles
        meteo_maille_jour <- coordinates %>%
          left_join(meteo, by = c("index_maille" = "idx_maille"), relationship = "many-to-many", multiple = "all")
        meteo_maille_jour <- meteo_maille_jour %>% select(!index_maille)
        # Récup des années des projections climatiques
        meteo_maille_jour <- meteo_maille_jour %>%
          mutate(date = as.Date((date), format = "%Y%m%d"), YEAR = format(date, "%Y"))
        projections <- sim_esca(meteo_maille_jour, description_parcelles, info_mailles)
        save(projections, file = paste0("data/Projections/Res/esca_",i_scen,"_",i_model, "_", period, ".RData"))
        }
      }
    }
}
 

# Ne simuler et ne récupérer que les variables écoclimatiques 
for(i_model in cmip6_models){
  print(paste0("*** J'en suis au modèle : ", i_model, " ***"))
  for(i_scen in c("ssp245", "ssp585")){ # "historical", 
    if(i_scen == "ssp245" & i_model == "KIOST-ESM") break
    print(paste0("Et au scénario : ", i_scen))
    
    my_period <- if(i_scen ==  "historical") "1995-2014" else c("2046-2065") # "2081-2100"
    
    #  Check if pr tasmax and tasmin exist
    if(paste0("tasmax_",i_scen,"_",i_model, ".RData") %in% cmip6_files & 
       paste0("tasmin_",i_scen,"_",i_model, ".RData") %in% cmip6_files &
       paste0("pr_",i_scen,"_",i_model, ".RData") %in% cmip6_files){
      # load data
      tm <- get(load(paste0(repository, "tasmin_",i_scen,"_",i_model, ".RData")))[,mailles2keep]
      date <- format(as.Date(rownames(tm)), "%Y%m%d")
      for(period in my_period){
        if(i_scen == "historical"){ idx_date = which(date >= 19940000 & date < 20150000) 
        } else if(period == "2046-2065"){ idx_date = which(date >= 20450000 & date < 20660000) # de 2046 à 2065
        } else idx_date = which(date >= 20800000 & date < 21010000) # de 2081 à 2100
        date <- date[idx_date]
        tm <- tm[idx_date,]
        
        tx <- get(load(paste0(repository, "tasmax_",i_scen,"_",i_model, ".RData")))[idx_date,mailles2keep]
        pr <- get(load(paste0(repository, "pr_",i_scen,"_",i_model, ".RData")))[idx_date,mailles2keep]
        hu <- get(load(paste0(repository, "hurs_",i_scen,"_",i_model, ".RData")))[idx_date,mailles2keep]
        rm(myqmap_data); gc()
        
        ### Joindre toutes les projections climatiques
        meteo <- data.frame(
          date = rep(date, times = length(mailles2keep)),
          idx_maille = rep(mailles2keep, each = length(date)),
          tm = as.vector(tm) - 273.15,
          tx = as.vector(tx) - 273.15,
          pr = as.vector(pr),
          hu = as.vector(hu)
        )
        
        
        # Si la pluie est négative la mettre à 0 :
        meteo$pr <- pmax(meteo$pr, 0)
        # Si l'humidité est sup à 100 la mettre à 100
        meteo$hu <- pmin(meteo$hu, 100)
        # Si la température max est inférieure à la température min : la mettre à la températue min
        meteo$tx <- pmax(meteo$tx, meteo$tm)
        # Si la température min est supérieure à la température min : la mettre à la températue max
        meteo$tm <- pmin(meteo$tx, meteo$tm)
        
        ### Joindre les données climatiques aux mailles
        meteo_maille_jour <- coordinates %>%
          left_join(meteo, by = c("index_maille" = "idx_maille"), relationship = "many-to-many", multiple = "all")
        meteo_maille_jour <- meteo_maille_jour %>% select(!index_maille)
        # Récup des années des projections climatiques
        meteo_maille_jour <- meteo_maille_jour %>%
          mutate(date = as.Date((date), format = "%Y%m%d"), YEAR = format(date, "%Y"))
        projections_climat <- sim_climat(meteo_maille_jour, description_parcelles, info_mailles)
        save(projections_climat, file = paste0("data/Projections/Climat/climat_",i_scen,"_",i_model, "_", period, ".RData"))
      }
    }
  }
}
  


