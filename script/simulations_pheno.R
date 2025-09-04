# 
# Script Sébastien ZITO ; sebastien.zito@inrae.fr
#        Gabriel MACE ; gabriel.mace@inrae.fr
#        Chloe DELMAS ; chloe.delmas@inrae.fr
# 
#  avril 2025
# 
# Ce script à pour but de simuler les dates de débourrement et de floraison pour nos différents cépages
# et nos différentes mailles en fonction des observations climatiques correspondantes (2003-2023)
# Mais aussi en fonction des simulations de climat RCP 4-5 et RCP 8-5

library(readxl)
library(dplyr)
library(Metrics) # Qualité des prédictions
library(ggplot2)

# Chargement des fonctions de phénologies
source("R_functions/fonctions_pheno.R")



# Simulations phénologies sur nos données observées SAFRAN/ ESCA ----
load(file = "data/climat/climat_maille_cepage_jour_complet.RData") # Données SAFRAN : "complet"
# Filtre des cépages et des années entre 2003 et 2023
# Enlever le vermentino car peu d'individus, 
# et le sciccarellu et le niellucciu car pas de floraison simulable
complet <- complet %>% filter(!cepage %in% c("sciaccarellu", "vermentinu", "niellucciu"),
                              DATE >= 20030000 & DATE < 20240000)
complet$doy <- as.numeric(strftime(as.Date(as.character(complet$DATE),"%Y%m%d"), format = "%j"))
complet$year <- as.numeric(strftime(as.Date(as.character(complet$DATE),"%Y%m%d"), format = "%Y"))

# Simulation des dates de débourrement :
date_debourrement <- data.frame(maille = NULL, cepage = NULL, year = NULL, debourrement_SUWE = NULL)
for(maille in unique(complet$ID_MAILLE)){ # Pour chaque maille
  for(cepage in unique(complet[complet$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    data <- complet[complet$ID_MAILLE == maille & complet$cepage == cepage ,] # Données correspondantes
    pheno <- debourrement(cepage = cepage, tn = data$TINF_H_Q, tx = data$TSUP_H_Q, doy = data$doy, 
                          year = data$year)
    date_debourrement <- rbind(date_debourrement, data.frame(maille = maille, cepage = cepage, 
                                                             year = pheno$year, debourrement_SUWE = pheno$debourrement_SUWE))
  }
}

# Simulation des dates de floraison :
date_floraison <- data.frame(maille = NULL, cepage = NULL, year = NULL, flo_gfv = NULL, flo_SUWE = NULL)
for(maille in unique(complet$ID_MAILLE)){ # Pour chaque maille
  for(cepage in unique(complet[complet$ID_MAILLE == maille,]$cepage)){ # Pour cahque cépage dans la maille
    data <- complet[complet$ID_MAILLE == maille & complet$cepage == cepage,] # Données correspondantes
    pheno <- floraison(cepage = cepage, tn = data$TINF_H_Q, tx = data$TSUP_H_Q, doy = data$doy, 
                      year = data$year)
    date_floraison <- rbind(date_floraison, data.frame(maille = maille, cepage = cepage, 
                                                       year = pheno$year, flo_gfv = pheno$flo_gfv,
                                                       flo_SUWE = pheno$flo_SUWE))
  }
}


# Aggrégation des deux dates
date_pheno <- date_floraison %>% left_join(date_debourrement, by = c("maille", "cepage", "year"), 
                                           relationship = "one-to-one")

# Stats descriptives ----
# Récupération des données ESCA
esca <- read.csv("data/ESCA/esca_filtre.csv") # Données ESCA
esca <- esca %>% filter(!cepage %in% c("sciaccarellu", "vermentinu", "niellucciu")) 
# Jointure  maille / commune
maille_commune <- read.csv("data/geographie/correspondance_maille_commune.csv", sep = ",")
maille_cepage_annee <- esca %>% left_join(maille_commune[, 
                                            c("ID_MAILLE", "INSEE_COM","X_maille", "Y_maille")],
                                            by = c("code_commune_INSEE" = "INSEE_COM"),
                                            relationship = "many-to-many", multiple = "any")
# Nombre d'observations si on prends une observations par maille, cépage, année :
dim(maille_cepage_annee %>% distinct(cepage, annee, ID_MAILLE))[1] # ESCA : 5489 obs -> 3443
# Table des maille, cépage, année avec plusieurs observations (parcelles) :
count <- maille_cepage_annee %>% group_by(annee, ID_MAILLE, cepage) %>% summarise(n = n())
table(count$n)
rm(count)
gc()
# Enregistrer en RData:
maille_cepage_annee <- maille_cepage_annee %>% 
  select(annee, code_commune_INSEE, cepage, ID_MAILLE, X_maille, Y_maille)
save(maille_cepage_annee, file = "data/maille_cepage_annee.RData")

# Jointure phénologie / incidence esca par cépage :
incidence_cepage <- esca %>% 
  group_by(cepage) %>%
  summarise(moy_esca = round(mean(pourcentage_esca),2)) 
date_pheno <- date_pheno %>% left_join(incidence_cepage, by = "cepage")

# Ajout des régions et départements pour la visualisation : 
maille_commune <- maille_commune %>% left_join(esca[, c("code_commune_INSEE", "nom_region",
                                                        "nom_departement")],
                                               by = c("INSEE_COM" = "code_commune_INSEE"),
                                               relationship = "one-to-one", multiple = "any")

date_pheno <- date_pheno %>% left_join(maille_commune[, c("ID_MAILLE", "nom_region", "nom_departement")],
                                       by = c("maille" = "ID_MAILLE"),
                                       relationship = "many-to-many", multiple = "any")

# Débourrement par cépage :
data_counts <- date_pheno %>% filter(!is.na(debourrement_SUWE)) %>%
  group_by(cepage) %>%
  summarise(n = n(), .groups = "drop")
ggplot(date_pheno %>% filter(!is.na(debourrement_SUWE)), aes(
  x = reorder(cepage, debourrement_SUWE), y = debourrement_SUWE, fill = moy_esca)) + 
  geom_boxplot() + 
  geom_text(data = data_counts, 
            aes(x = as.factor(cepage), 
            y =  max(date_pheno$debourrement_SUWE, na.rm = TRUE) + 2, label = n, fill = 1), 
            position = position_dodge(width = 0.75), size = 3) +  
  theme(axis.text.x = element_text(angle = 90, size = 11),
        text = element_text(size = 10)) + 
  scale_fill_gradient2(midpoint = median(incidence_cepage$moy_esca), 
                       low = "blue", mid = "white", high = "red", space = "Lab" ) +
  labs(x = "Cépage", y = "Jour de débourrement", fill = "Incidence esca moyenne",
       title = "Jour de débourrement simulé par cépage, des couples mailles / cépages d'intêrets ESCA, 2003-2023")

# Débourrement par année :
data_counts <- date_pheno %>% filter(!is.na(debourrement_SUWE)) %>%
  group_by(year) %>%
  summarise(n = n(), .groups = "drop")
date_pheno %>% group_by(year) %>%
  summarize(debourrement_SUWE = mean(debourrement_SUWE, na.rm = TRUE)) %>%
  ggplot(aes(x = year, y = debourrement_SUWE)) + geom_line() +
  geom_text(data = data_counts, 
            aes(x = year, 
                y =  mean(date_pheno$debourrement_SUWE, na.rm = TRUE) + 10, label = n), 
            position = position_dodge(width = 0.75), size = 3) + 
  labs(x = "Année", y = "Jour de débourrement SUWE moyen", 
       title = "Jour de débourrement simulé moyen par année, des couples mailles / cépages d'intêrets ESCA, 2003-2023")

# Floraison GFV par cépage :
data_counts <- date_pheno %>%
  group_by(cepage) %>%
  summarise(n = n(), .groups = "drop")
ggplot(date_pheno, aes(
  x = reorder(cepage, flo_gfv, mean), y = flo_gfv, fill = moy_esca)) + 
  geom_boxplot() + 
  geom_text(data = data_counts, 
            aes(x = as.factor(cepage), 
                y =  max(date_pheno$flo_gfv, na.rm = TRUE) + 2, label = n, fill = 1), 
            position = position_dodge(width = 0.75), size = 3) + 
  theme(axis.text.x = element_text(angle = 90, size = 11),
        text = element_text(size = 10), ) + 
  scale_fill_gradient2(midpoint = median(incidence_cepage$moy_esca), 
                       low = "blue", mid = "white", high = "red", space = "Lab" ) +
  labs(x = "Cépage", y = "Jour de floraison gfv", fill = "Incidence esca moyenne",
       title = "Jour de floraison simulé (GFV) par cépage, des couples mailles / cépages d'intêrets ESCA, 2003-2023")

# FLoraison GFV par année :
data_counts <- date_pheno %>%
  group_by(year) %>%
  summarise(n = n(), .groups = "drop")
date_pheno %>% group_by(year) %>%
  summarize(flo_gfv = mean(flo_gfv, na.rm = TRUE)) %>%
  ggplot(aes(x = year, y = flo_gfv)) + geom_line() +
  geom_text(data = data_counts, 
            aes(x = year, 
                y =  mean(date_pheno$flo_gfv, na.rm = TRUE) + 12, label = n), 
            position = position_dodge(width = 0.75), size = 3) + 
  labs(x = "Année", y = "Jour de floraison GFV moyen", 
       title = "Jour de floraison simulé (GFV) moyen par année, des couples mailles / cépages d'intêrets ESCA, 2003-2023")


# Floraison SUWE par cépage :
data_counts <- date_pheno %>% filter(!is.na(flo_SUWE)) %>%
  group_by(cepage) %>%
  summarise(n = n(), .groups = "drop")
ggplot(subset(date_pheno, !is.na(flo_SUWE)), aes(
  x = reorder(cepage, flo_SUWE, mean), y = flo_SUWE, fill = moy_esca)) + 
  geom_boxplot() + 
  geom_text(data = data_counts, 
            aes(x = as.factor(cepage), 
                y =  max(date_pheno$flo_SUWE, na.rm = TRUE) + 2, label = n, fill = 1), 
            position = position_dodge(width = 0.75), size = 3) +
  theme(axis.text.x = element_text(angle = 90, size = 11),
        text = element_text(size = 10), ) + 
  scale_fill_gradient2(midpoint = median(incidence_cepage$moy_esca), 
                       low = "blue", mid = "white", high = "red", space = "Lab" ) +
  labs(x = "Cépage", y = "Jour de floraison SUWE", fill = "Incidence esca moyenne",
       title = "Jour de floraison simulé (SUWE) par cépage, des couples mailles / cépages d'intêrets ESCA, 2003-2023")

# Floraison SUWE par année :
data_counts <- date_pheno %>% filter(!is.na(flo_SUWE)) %>%
  group_by(year) %>%
  summarise(n = n(), .groups = "drop")
date_pheno %>% group_by(year) %>%
  summarize(flo_SUWE = mean(flo_SUWE, na.rm = TRUE)) %>%
  ggplot(aes(x = year, y = flo_SUWE)) + geom_line() +
  geom_text(data = data_counts, 
            aes(x = year, 
                y =  mean(date_pheno$flo_SUWE, na.rm = TRUE) + 16, label = n), 
            position = position_dodge(width = 0.75), size = 3) + 
  labs(x = "Année", y = "Jour de floraison SUWE moyen", 
       title = "Jour de floraison simulé (SUWE) moyen par année, des couples mailles / cépages d'intêrets ESCA, 2003-2023")



# Pour chaque cépage : boxplot par an : 
for(cepage in unique(date_pheno$cepage)){
  plot(ggplot(date_pheno %>% filter(cepage == cepage), aes(
    x = as.factor(year), y = debourrement_SUWE)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, size = 11),
          text = element_text(size = 10)) + 
    labs(x = "Année", y = "Jour de débourrement",
         title = paste("Jour de débourrement simulé par an, pour chaque maille où", cepage, "est présent")))
  plot(ggplot(date_pheno %>% filter(cepage == cepage), aes(
    x = as.factor(year), y = flo_gfv)) + 
      geom_boxplot() + 
      theme(axis.text.x = element_text(angle = 90, size = 11),
            text = element_text(size = 10)) + 
      labs(x = "Année", y = "Jour de floraison",
           title = paste("Jour de floraison simulé (gfv) par an, pour chaque maille où", cepage, "est présent")))
}

for(cepage in unique(date_pheno[!is.na(date_pheno$flo_SUWE),]$cepage)){
plot(ggplot(date_pheno %>% filter(cepage == cepage), aes(
    x = as.factor(year), y = flo_SUWE)) + 
      geom_boxplot() + 
      theme(axis.text.x = element_text(angle = 90, size = 11),
            text = element_text(size = 10)) + 
      labs(x = "Année", y = "Jour de floraison",
           title = paste("Jour de floraison simulé (SUWE) par an, pour chaque maille où", cepage, "est présent")))
}

# ggplot(date_pheno, aes(
#   x = reorder(nom_region, debourrement_SUWE, mean), y = debourrement_SUWE)) + 
#   geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90, size = 11),
#         text = element_text(size = 10), ) + 
#   labs(x = "Région", y = "Jour de débourrement", fill = "Incidence esca moyenne",
#        title = "Jour de débourrement par région, mailles / cépages d'intêrets ESCA, 2003-2023")
# 
# ggplot(date_pheno, aes(
#   x = reorder(nom_region, flo_gfv, mean), y = flo_gfv)) + 
#   geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90, size = 11),
#         text = element_text(size = 10), ) + 
#   labs(x = "Région", y = "Jour de floraison gfv", fill = "Incidence esca moyenne",
#        title = "Jour de floraison gfv par région, mailles / cépages d'intêrets ESCA, 2003-2023")

rm(data_counts, data, complet, date_pheno, esca, incidence_cepage, maille_cepage_annee, maille_commune, 
   cepage, maille)
gc()

# Enregistrement ----
# Remise à zéro :
date_pheno <- date_floraison %>% left_join(date_debourrement, by = c("maille", "cepage", "year"), 
                                           relationship = "one-to-one")
# Enregister les simulations de phénologies de 2003 à 2023 :
save(date_pheno, file = "data/phenologie/phenologie_2003_2023.RData") 


# Aggrégation esca / phénologies :
# Récupération des données ESCA
esca <- read.csv("data/ESCA/esca_filtre.csv") # Données ESCA
esca <- esca %>% filter(!cepage %in% c("sciaccarellu", "vermentinu", "niellucciu")) 
esca <- esca %>% select(!c("identifiant_parcelle_analyse", "nom_commune_analyse", "organisme_notateur"))
# Jointure  maille / commune
maille_commune <- read.csv("data/geographie/correspondance_maille_commune.csv", sep = ",")

esca <- esca %>% left_join(maille_commune[,c("ID_MAILLE", "INSEE_COM")],
                           by = c("code_commune_INSEE" = "INSEE_COM"),
                           relationship = "many-to-many", multiple = "any")

esca_pheno <- esca %>% left_join(date_pheno,
                           by = c("ID_MAILLE" = "maille", "annee" = "year", "cepage" = "cepage"),
                           relationship = "many-to-one", multiple = "any")

# Enregister la première partie des données de modélisation :
save(esca_pheno, file = "data/modelisation/esca_phenologie.RData") 

#############################################################################################"

# Simulations phénologies via DRIAS 2006-2023 ----
load("data/maille_cepage_annee.RData")
maille_cepage <- maille_cepage_annee %>% select(cepage, ID_MAILLE) %>% distinct()

# 4_5 : 
load("data/climat/previsions_4_5_2006_2023.RData") 
# Simulation des dates de débourrement et de floraison :
debourrement_4_5_passe <- data.frame(maille = NULL, cepage = NULL, year = NULL, debourrement_SUWE = NULL)
flo_4_5_passe <- data.frame(maille = NULL, cepage = NULL, year = NULL, flo_gfv = NULL, flo_SUWE = NULL)
for(maille in unique(maille_cepage$ID_MAILLE)){ # Pour chaque maille
  for(cepage in unique(maille_cepage[maille_cepage$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    data <- drias_4_5_2006_2023[drias_4_5_2006_2023$ID_MAILLE == maille,] # Données correspondantes
    deb <- debourrement(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    debourrement_4_5_passe <- rbind(debourrement_4_5_passe,
                                      data.frame(maille = maille, cepage = cepage, year = deb$year, 
                                                 debourrement_SUWE = deb$debourrement_SUWE))
    flo <- floraison(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    flo_4_5_passe <- rbind(flo_4_5_passe, data.frame(maille = maille, cepage = cepage, 
                                                       year = flo$year, flo_gfv = flo$flo_gfv,
                                                       flo_SUWE = flo$flo_SUWE))
  }
}

# Aggrégation des deux dates
pheno_4_5_passe <- flo_4_5_passe %>% left_join(debourrement_4_5_passe, by = c("maille", "cepage", "year"), 
                                           relationship = "one-to-one")
#Enregistrement sous forme de RDdata :
save(pheno_4_5_passe, file = "data/phenologie/phenologie_4_5_passe.RData")

# 8_5 : 
load("data/climat/previsions_8_5_2006_2023.RData")
# Simulation des dates de débourrement et de floraison :
debourrement_8_5_passe <- data.frame(maille = NULL, cepage = NULL, year = NULL, debourrement_SUWE = NULL)
flo_8_5_passe <- data.frame(maille = NULL, cepage = NULL, year = NULL, flo_gfv = NULL, flo_SUWE = NULL)
for(maille in unique(maille_cepage$ID_MAILLE)){ # Pour chaque maille
  for(cepage in unique(maille_cepage[maille_cepage$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    data <- drias_8_5_2006_2023[drias_8_5_2006_2023$ID_MAILLE == maille,] # Données correspondantes
    deb <- debourrement(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    debourrement_8_5_passe <- rbind(debourrement_8_5_passe,
                                    data.frame(maille = maille, cepage = cepage, year = deb$year, 
                                               debourrement_SUWE = deb$debourrement_SUWE))
    flo <- floraison(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    flo_8_5_passe <- rbind(flo_8_5_passe, data.frame(maille = maille, cepage = cepage, 
                                                     year = flo$year, flo_gfv = flo$flo_gfv,
                                                     flo_SUWE = flo$flo_SUWE))
  }
}

# Aggrégation des deux dates
pheno_8_5_passe <- flo_8_5_passe %>% left_join(debourrement_8_5_passe, by = c("maille", "cepage", "year"), 
                                               relationship = "one-to-one")
#Enregistrement sous forme de RDdata :
save(pheno_8_5_passe, file = "data/phenologie/phenologie_8_5_passe.RData")


# Simulations prévisions phénologies via DRIAS 2024-2100 ----
# 4_5 : 
load("data/climat/previsions_4_5_2024_2100.RData") 
# Simulation des dates de débourrement et de floraison :
debourrement_4_5_futur <- data.frame(maille = NULL, cepage = NULL, year = NULL, debourrement_SUWE = NULL)
flo_4_5_futur <- data.frame(maille = NULL, cepage = NULL, year = NULL, flo_gfv = NULL, flo_SUWE = NULL)
for(maille in unique(maille_cepage$ID_MAILLE)){ # Pour chaque maille
  for(cepage in unique(maille_cepage[maille_cepage$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    data <- drias_4_5[drias_4_5$ID_MAILLE == maille,] # Données correspondantes
    deb <- debourrement(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    debourrement_4_5_futur <- rbind(debourrement_4_5_futur,
                                    data.frame(maille = maille, cepage = cepage, year = deb$year, 
                                               debourrement_SUWE = deb$debourrement_SUWE))
    flo <- floraison(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    flo_4_5_futur <- rbind(flo_4_5_futur, data.frame(maille = maille, cepage = cepage, 
                                                     year = flo$year, flo_gfv = flo$flo_gfv,
                                                     flo_SUWE = flo$flo_SUWE))
  }
}

# Aggrégation des deux dates
pheno_4_5_futur <- flo_4_5_futur %>% left_join(debourrement_4_5_futur, by = c("maille", "cepage", "year"), 
                                               relationship = "one-to-one")
#Enregistrement sous forme de RDdata :
save(pheno_4_5_futur, file = "data/phenologie/phenologie_4_5_futur.RData")

# 8_5 :
load("data/climat/previsions_8_5_2024_2100.RData") 
# Simulation des dates de débourrement et de floraison :
debourrement_8_5_futur <- data.frame(maille = NULL, cepage = NULL, year = NULL, debourrement_SUWE = NULL)
flo_8_5_futur <- data.frame(maille = NULL, cepage = NULL, year = NULL, flo_gfv = NULL, flo_SUWE = NULL)
for(maille in unique(maille_cepage$ID_MAILLE)){ # Pour chaque maille
  for(cepage in unique(maille_cepage[maille_cepage$ID_MAILLE == maille,]$cepage)){ # Pour chaque cépage dans la maille
    data <- drias_8_5[drias_8_5$ID_MAILLE == maille,] # Données correspondantes
    deb <- debourrement(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    debourrement_8_5_futur <- rbind(debourrement_8_5_futur,
                                    data.frame(maille = maille, cepage = cepage, year = deb$year, 
                                               debourrement_SUWE = deb$debourrement_SUWE))
    flo <- floraison(cepage = cepage, tn = data$TN, tx = data$TX, doy = data$doy, year = data$year)
    flo_8_5_futur <- rbind(flo_8_5_futur, data.frame(maille = maille, cepage = cepage, 
                                                     year = flo$year, flo_gfv = flo$flo_gfv,
                                                     flo_SUWE = flo$flo_SUWE))
  }
}

# Aggrégation des deux dates
pheno_8_5_futur <- flo_8_5_futur %>% left_join(debourrement_8_5_futur, by = c("maille", "cepage", "year"), 
                                               relationship = "one-to-one")
#Enregistrement sous forme de RDdata :
save(pheno_8_5_futur, file = "data/phenologie/phenologie_8_5_futur.RData")
