# 
# Script Sébastien ZITO ; sebastien.zito@inrae.fr
#        Gabriel MACE ; gabriel.mace@inrae.fr
#        Chloe DELMAS ; chloe.delmas@inrae.fr
# 
#  mars 2025
# 
# 
#

# 
library(data.table)
library(terra) # données spatiales
library(RPostgres) # lecture BDDesca
library(readxl) # lecture BDDesca
library(dplyr)
library(ggplot2)
library(stringr) # manipulation de caractères
library(tidyterra) # left_join for spatvector
library(meteoRIT) # Package fait par INRAE contenant les altitudes des mailles


# 
# # -------------------------------------------------------------------------- #
# # Connect BDD ESCA ----
# # -------------------------------------------------------------------------- #
# Penser à activer le VPN globalProtect
con_PROPRE <- RPostgres::dbConnect(
  RPostgres::Postgres(),
  host = "bdd-private.plateforme-esv.fr",
  port = "5432",
  dbname = "PROPRE", # remplacer par le nom de la BDD (BRUT, PROPRE...)
  user = Sys.getenv("USERNAME_BDD"),
  password = Sys.getenv("PASSWORD_BDD")
)
#
# variables à garder (nom colonnes) dans la BBD
var2kip <- c("annee", "code_commune_INSEE", "nom_region", "cepage", "nombre_ceps_total", "nombre_morts",
             "nombre_jeunes_complants", "nombre_esca_bda_partiel_total_somme",
             "identifiant_parcelle_analyse", "age_parcelle_estime",
             "pourcentage_esca", "pourcentage_jeunes",
             "nom_departement", "nom_commune_analyse",
             "organisme_notateur", "nom_bassin_viti2")
esca <- RPostgres::dbReadTable(con_PROPRE, Id(schema = "ESCA",table = "BDD_ESCA_MAJ"))[,var2kip]
# 
# # -------------------------------------------------------------------------- #
# # Sélect Esca data ----
# # -------------------------------------------------------------------------- #
# 
# Règle de décisisons :
# 1 : limites de l'age
age_lim <- c(10,30)

# 2 : Minimum d'années dans la base ensuite
min_years <- 3 
# 
esca <- esca[!is.na(esca$age_parcelle_estime) & esca$age_parcelle_estime >= age_lim[1] 
             & esca$age_parcelle_estime <= age_lim[2],] # age 10-30 ans
count_parcelles <- esca %>% group_by(identifiant_parcelle_analyse) %>% summarise(n = n())
parcelles_bonnes <- count_parcelles[count_parcelles$n >= min_years,]$identifiant_parcelle_analyse
esca <- esca[esca$identifiant_parcelle_analyse %in% parcelles_bonnes,]
# 

# Ajout de la couleur de chaque cépage
couleur_vin <- data.frame(cepage = 
                            c("Cabernet sauvignon", "Meunier", "Merlot", "Pinot noir", "Trousseau",
                              "Gamay", "Poulsard", "Muscat de hambourg", "Grenache", "Sciaccarellu",
                              "Syrah", "Cinsault", "Carignan", "Mourvèdre", "Muscat petits grains",
                              "Niellucciu", "Fer servadou",
                              "Chardonnay", "Sauvignon blanc", "Chenin", "Cabernet franc", "Melon",
                              "Ugni blanc", "Gewurztraminer", "Riesling", "Savagnin", "Semillon",
                              "Pinot auxerrois", "Vermentinu"),
                          couleur = c(rep("Cépage noir", 17), rep("Cépage blanc", 12)))

esca <- esca %>% left_join(couleur_vin, by = "cepage")

# Transformation du nom des cépages
esca$cepage <- str_replace_all(str_to_lower(esca$cepage), "è", "e")

rm(parcelles_bonnes) ; rm(count_parcelles, couleur_vin) ; gc()
# 
# # -------------------------------------------------------------------------- #
# # Stat descriptives ----
# # -------------------------------------------------------------------------- #
# 
# 
png("graphs/stats_descriptives/nb_parcelle_annee_cepage.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  count(cepage, couleur) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(cepage,n), fill = couleur)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("Cépage noir" = "darkviolet", "Cépage blanc" = "gold")) + 
  labs(title = "Nombre de parcelle_année par cépage",
       x = "Nombre d'observations",
       y = "Cépage")
dev.off()
# 
png("graphs/stats_descriptives/nb_parcelle_annee_region.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  count(nom_region) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(nom_region,n))) + 
  geom_bar(fill = "darkgrey", stat = "identity") +
  labs(title = "Nombre de parcelle_année par région",
       x = "Nombre d'observations",
       y = "Région")
dev.off()
# 
png("graphs/stats_descriptives/nb_parcelle_par_annee.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  count(annee) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(annee, annee))) + 
  geom_bar(fill = "grey", stat = "identity") +
  labs(title = "Nombre de parcelles par année",
       x = "Nombre d'observations",
       y = "Année")
dev.off()
# 
png("graphs/stats_descriptives/nb_parcelle_annee_age.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  count(age_parcelle_estime) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(age_parcelle_estime, age_parcelle_estime))) + 
  geom_bar(fill = "grey", stat = "identity") +
  labs(title = "Nombre de parcelle_année par âge",
       x = "Nombre d'observations",
       y = "Âge")
dev.off()
# 
png("graphs/stats_descriptives/incidence_cepage.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  group_by(cepage, couleur) %>%
  summarise(moy_esca = round(mean(pourcentage_esca),2), se_esca = round(sd(pourcentage_esca)/sqrt(n()),2)) %>%
  ggplot(aes(x = moy_esca, y = reorder(cepage, moy_esca), fill = couleur)) + 
  geom_bar(stat = "identity")  +
  scale_fill_manual(values=c("Cépage noir" = "darkviolet", "Cépage blanc" = "gold")) + 
  geom_errorbar( aes(y=reorder(cepage, moy_esca), xmin=moy_esca-se_esca, xmax=moy_esca+se_esca), width=0.4) +
  geom_text(aes(x = moy_esca + se_esca + 0.5, label= moy_esca), color="black", size=4, ) +
  labs(title = "Incidence moyenne (+/- erreur type) d'esca par cépage",
       x = "Incidence moyenne",
       y = "Cépage")
dev.off()
# 
png("graphs/stats_descriptives/incidence_region.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  group_by(nom_region) %>%
  summarise(moy_esca = round(mean(pourcentage_esca),2), se_esca = round(sd(pourcentage_esca)/sqrt(n()),2)) %>%
  ggplot(aes(x = moy_esca, y = reorder(nom_region, moy_esca))) + 
  geom_bar(fill = "darkgrey", stat = "identity") +
  geom_errorbar( aes(y=reorder(nom_region, moy_esca), xmin=moy_esca-se_esca, xmax=moy_esca+se_esca), width=0.4) +
  geom_text(aes(x = moy_esca + se_esca + 0.5, label= moy_esca), color="black", size=4, ) +
  labs(title = "Incidence moyenne (+/- erreur type) d'esca par région",
       x = "Incidence moyenne",
       y = "Région")
dev.off()
# 
png("graphs/stats_descriptives/incidence_annee.png",width = 10, height = 7, units = "in", res = 100)
esca %>% 
  group_by(annee) %>%
  summarise(moy_esca = round(mean(pourcentage_esca),2), se_esca = round(sd(pourcentage_esca)/sqrt(n()),2)) %>%
  ggplot(aes(x = moy_esca, y = reorder(annee, annee))) + 
  geom_bar(fill = "grey", stat = "identity") +
  geom_errorbar( aes(y=reorder(annee, annee), xmin=moy_esca-se_esca, xmax=moy_esca+se_esca), width=0.4) +
  geom_text(aes(x = moy_esca + se_esca +0.5, label= moy_esca), color="black", size=4, ) +
  labs(title = "Incidence moyenne (+/- erreur type) d'esca par année",
       x = "Incidence moyenne",
       y = "Année")
dev.off()
# 
png("graphs/stats_descriptives/incidence_age.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  group_by(age_parcelle_estime) %>%
  summarise(moy_esca = round(mean(pourcentage_esca),2), se_esca = round(sd(pourcentage_esca)/sqrt(n()),2)) %>%
  ggplot(aes(x = moy_esca, y = reorder(age_parcelle_estime, age_parcelle_estime))) + 
  geom_bar(fill = "grey", stat = "identity") +  
  geom_errorbar( aes(y=reorder(age_parcelle_estime, age_parcelle_estime), xmin=moy_esca-se_esca, xmax=moy_esca+se_esca), width=0.4) +
  geom_text(aes(x = moy_esca + se_esca + 0.5, label= moy_esca), color="black", size=4, ) +
  labs(title = "Incidence moyenne (+/- erreur type) d'esca par âge",
       x = "Incidence moyenne",
       y = "Âge")
dev.off()
# 
png("graphs/stats_descriptives/incidence_nb_ceps.png", width = 10, height = 7, units = "in", res = 100)
esca %>% 
  group_by(nombre_ceps_total) %>%
  summarise(moy_esca = mean(pourcentage_esca)) %>%
  ggplot(aes(x = nombre_ceps_total, y = moy_esca)) + 
  geom_point() +
  labs(title = "Incidence moyenne / nombre de ceps",
       x = "Nombre ceps",
       y = "Incidence moyenne")
dev.off()
# 
# # -------------------------------------------------------------------------- #
# # BDD Communes shp ----
# # -------------------------------------------------------------------------- #
# 
# Communes shp file
# Download here : A metre à jour avec le bon lien
# 
communes <- vect("data/geographie/COMMUNE_FRMETDROM.shp")
communes <- project(communes, "EPSG:2154")
com_asso <- vect("data/geographie/COMMUNE_ASSOCIEE_OU_DELEGUEE.shp")
# 
# 69101 : Jarnioux devient commune déléguée au sein de Porte des Pierres Dorées (69159)
# 39395 : Orbagna devient commune déléguée au sein de Beaufort-Orbagna (39043) (commune nouvelle)
esca[esca$code_commune_INSEE == "69101",]$nom_commune_analyse <- "PORTEDESPIERRESDOREES"
esca[esca$code_commune_INSEE == "69101",]$code_commune_INSEE <- "69159"
esca[esca$code_commune_INSEE == "39395",]$nom_commune_analyse <- "BEAUFORTORBAGNA"
esca[esca$code_commune_INSEE == "39395",]$code_commune_INSEE <- "39043"


## Choix du au nombre fabile d'observations ou à l'impossibilité de simuler la phénologie:
esca <- esca %>% filter(!cepage %in% c("fer servadou","sciaccarellu", "vermentinu", "niellucciu")) 

# Transformation du nom des bassins viticoles :
# description_parcelles <- read.csv("data/Reg_viti_water_balance_param.csv")
# Diois : une parcelle 6 années de muscat petits grains en drôme -> Cotes-du-Rhone nord
esca$region_viticole <- esca$nom_bassin_viti2
esca[esca$nom_departement %in% c("Drome") & esca$nom_bassin_viti2 == "Diois",]$region_viticole <- "Cotes-du-Rhone nord"
# Ardeche : 07241 : Saint Germain = "Cotes-du-Rhone sud"
esca[esca$code_commune_INSEE == "07241",]$region_viticole <- "Cotes-du-Rhone sud"
esca[esca$nom_departement %in% c("Vaucluse", "Gard") & esca$nom_bassin_viti2 == "Vallée du Rhône",]$region_viticole <- "Cotes-du-Rhone sud"
write.csv(esca, "data/ESCA/esca_filtre.csv", row.names=FALSE)

# 
# identification des communes à retrouver dans la base selon le code INSEE
com2find <- as.character(unique(esca$code_commune_INSEE))
# 
# Vérifier que le code Insee contient bien 5 chiffres et ajouter un 0 devant si ce n'est pas le cas
com2find[nchar(com2find) != 5] <- paste0("0", com2find[nchar(com2find) != 5])
# 
# Extraction des communes dans la base de donnée communes
selected_com <- communes[communes$INSEE_COM %in% com2find,]
# 
# Extraction des communes dans la base de donnée communes_associee
miss_cod <-com2find[!com2find %in% communes$INSEE_COM]
selected_com_asso <- com_asso[com_asso$INSEE_CAD %in% miss_cod | com_asso$INSEE_COM %in% miss_cod,] # INSEE_CAD remplace INSEE_COD dans la BDD com_asso
#
# Vérifier que l'on à bien toutes les communes
if(length(selected_com) + length(selected_com_asso) != length(com2find)) { 
  print("Attention certaines communes n'ont pas été trouvées dans la base de données : ")
  print(com2find[!(com2find %in% selected_com$INSEE_COM | com2find %in% selected_com_asso$INSEE_CAD)])
}
# 
# Assemblage des communes et communes associee
names2kip <- c("NOM", "INSEE_COM")
selected_com <- selected_com[,names2kip]
if (nrow(selected_com_asso) >=1) {
  selected_com_asso <- selected_com_asso[,c("NOM", "INSEE_CAD")]
  names(selected_com_asso) <- names2kip
  All_com <- rbind(selected_com, selected_com_asso)
} else {
  All_com <- selected_com
}
# 
plot(All_com)
# 
rm(communes) ; rm(com_asso) ; gc()
# 
# # -------------------------------------------------------------------------- #
# #  Récup grille SAFRAN  ----
# # -------------------------------------------------------------------------- #
# 
# Get SAFRAN coordinates in Lambert II : EPSG 27572
# Réaliser une fois, puis enregistre dans data pour gagner du temps !!!
# 
# url <- paste0("https://object.files.data.gouv.fr/meteofrance/data/synchro_ftp/REF_CC/SIM/QUOT_SIM2_1980-1989.csv.gz")
# df <- fread(url, select = c("LAMBX", "LAMBY"))
# Saf_coords <- unique(paste0(df$LAMBX,"-",df$LAMBY))
# Saf_coords <- data.frame(X = as.numeric(sapply(strsplit(Saf_coords, "-"), function(x) x[1])),
#                          Y = as.numeric(sapply(strsplit(Saf_coords, "-"), function(x) x[2])))
# # Attention X et Y sont arrondis d'un facteur 100 dans la BDD 
# Saf_coords$X <- Saf_coords$X*100
# Saf_coords$Y <- Saf_coords$Y*100
# # 
# write.csv(Saf_coords, "data/SAFRAN_coords_EPSG27572.csv", row.names = F)
# # 
Saf_coords <- fread("data/geographie/SAFRAN_coords_EPSG27572.csv")
# 
#  Transformation au format spatial
coordinates <- vect(Saf_coords,geom = c("X", "Y"),crs = "EPSG:27572")  # Lambert II (EPSG:27572)
# 
#  Taille de la grille SAFRAN en m
cell_size <- 8000
# 
# Création de la grille 
safran_grid <- lapply(1:nrow(Saf_coords), function(i) {
  if (i %% 100 == 0) print(i)  
  x <- Saf_coords$X[i]
  y <- Saf_coords$Y[i]
  
  # Define the extents of the square around each centroid
  ext_xmin <- x - cell_size / 2
  ext_xmax <- x + cell_size / 2
  ext_ymin <- y - cell_size / 2
  ext_ymax <- y + cell_size / 2
  
  # Create coordinates for the vertices of the square
  square_coords <- matrix(c(
    ext_xmin, ext_ymin,  # Bottom-left corner
    ext_xmax, ext_ymin,  # Bottom-right corner
    ext_xmax, ext_ymax,  # Top-right corner
    ext_xmin, ext_ymax,  # Top-left corner
    ext_xmin, ext_ymin   # Close the polygon
  ), ncol = 2, byrow = TRUE)
  
  # Create a SpatPolygon from the square coordinates
  grid_cell <- vect(square_coords, type = "polygon", crs = "EPSG:27572")
  return(grid_cell)
})
safran_grid <- do.call(rbind, safran_grid)
# 
# Add X and Y columns
safran_grid <- cbind(safran_grid,data.frame(X = geom(coordinates)[,3], Y = geom(coordinates)[,4]))
# 
# Check projection, and put communes 
safran_grid <- project(safran_grid, "EPSG:2154")
# 
# Try if overlap
plot(safran_grid)
plot(All_com,col = "red", add = TRUE)
# 
# # -------------------------------------------------------------------------- #
# # Selection grille Esca ----
# # -------------------------------------------------------------------------- #
# 
# Selon la plus grande aire d'intersection : 
safran_grid$ID <- c(1:length(safran_grid))
intersections <- terra::intersect(safran_grid, All_com) # ici renvoie bien les polygones communes mais découpé par la grille (donc avec plus d'entité)
plot(intersections)
intersections$area <- expanse(intersections) # calc surface
max_area_per_polygon <- intersections[order(intersections$INSEE_COM, -intersections$area),] # tri selon la surface
max_area_per_polygon <- max_area_per_polygon[!duplicated(max_area_per_polygon$INSEE_COM),] # garde le 1er (plus grande surface)
selected_grid_cells <- safran_grid[safran_grid$ID %in% max_area_per_polygon$ID,]
plot(selected_grid_cells)
# 
# Selon la distance des centroïdes : 
near <- nearest(x = All_com, y = safran_grid, centroids = TRUE) # centroids les plus proches
near$from_id # communes
near$to_id # mailles
near$INSEE_COM <- All_com$INSEE_COM # récup les codes INSEE
near <- near[order(near$INSEE_COM),] # Faire en sorte que near et max_area soit dans le même ordre
plot(safran_grid[near$to_id], col = "blue")
# 
# Comparaison des deux méthodes : 
summary(max_area_per_polygon$INSEE_COM == near$INSEE_COM) # Vérifier la bonne correspondance des codes INSEE
summary(max_area_per_polygon$ID == near$to_id) # Nombre de maille différentes
# 
# Comparaison des 2 méthodes
plot(All_com[near$from_id[!max_area_per_polygon$ID == near$to_id]])
plot(safran_grid[max_area_per_polygon$ID[!max_area_per_polygon$ID == near$to_id]], add = TRUE, col = "blue")
plot(All_com[near$from_id[!max_area_per_polygon$ID == near$to_id]])
plot(safran_grid[near$to_id[!max_area_per_polygon$ID == near$to_id]], add = TRUE, col = "red")
# 
# CHOIX DE LA METHODE : Aire d'intersection
# 
# ! Attention : Les données X et Y safran sont à diviser par 100 par rapport à celles de météo France !
# Création d'un dataframe qui lie chaque commune à sa maille :
maille_commune <- data.frame(NOM_COM = max_area_per_polygon$NOM, INSEE_COM = max_area_per_polygon$INSEE_COM, 
                            ID_MAILLE = max_area_per_polygon$ID, X_maille = max_area_per_polygon$X/100, 
                            Y_maille = max_area_per_polygon$Y/100)


maille_commune <- maille_commune %>% 
  left_join(esca[, c("code_commune_INSEE", "region_viticole")] %>% unique(), 
            by = c("INSEE_COM" = "code_commune_INSEE"),
            relationship = "one-to-one")

# Ecriture des mailles SAFRAN, avec si ont les gardes ou non :
in_esca <- (1:9892 %in% maille_commune$ID_MAILLE)
safran_grid$in_esca <- in_esca
writeVector(safran_grid, "data/geographie/maille_safran.shp", overwrite = TRUE)

# Enquête du nombre de parcelles par maille ----

# Récupération contour France : 
limites_france <- vect("data/geographie/france_poly_limits.shp")
limites_france <- project(limites_france, "EPSG:2154")


ESCA_maille <- esca %>% 
  left_join(maille_commune[, c("INSEE_COM", "ID_MAILLE")], 
            by = c("code_commune_INSEE" = "INSEE_COM"),
            relationship = "many-to-one")

# 
png("graphs/stats_descriptives/densite_parcelle_maille.png", width = 10, height = 7, units = "in", res = 100)
ESCA_maille %>% distinct(identifiant_parcelle_analyse, .keep_all = TRUE) %>% 
  group_by(ID_MAILLE) %>% summarise(n = n()) %>% arrange(desc(n)) %>%
  ggplot(aes(x=n, after_stat(count))) + 
  geom_density(position = "stack") +
  xlab("Nombre de parcelle unique par maille")+
  ggtitle("Densité de parcelle unique par maille")
dev.off()
# 
parcelle_maille <- ESCA_maille %>% distinct(identifiant_parcelle_analyse, .keep_all = TRUE) %>% 
  group_by(ID_MAILLE, nom_region) %>% summarise(n = n()) %>% arrange(desc(n))
selected_grid_cells <- selected_grid_cells %>%
  left_join(parcelle_maille, by = c("ID" = "ID_MAILLE")) %>% rename(nb_parcelle = n)
# 

png("graphs/stats_descriptives/carte_parcelles_maille.png", width = 10, height = 7, units = "in", res = 100)
ggplot() +
  geom_spatvector(data = limites_france, fill = NA, color = "black") +
  geom_spatvector(data = selected_grid_cells, aes(fill = nb_parcelle))  +
  scale_fill_gradient(
    low = "blue",
    high = "red",
    na.value = "grey50") +
  labs(fill = "Nombre de parcelles")
dev.off()
# 
# Création d'un dataframe qui pour chaque cépage donne la commune et la maille où il existe : 
maille_cepage <- unique(esca[,c("code_commune_INSEE", "cepage")])
maille_cepage$code_commune_INSEE <- as.character(maille_cepage$code_commune_INSEE)
maille_cepage <- maille_cepage %>% left_join(maille_commune, by = c("code_commune_INSEE" = "INSEE_COM"))
### Enlever les communes dans ce df et supprimer les doublons maille / cépage
maille_cepage <- maille_cepage[,c("ID_MAILLE","cepage", "X_maille", "Y_maille")] %>% distinct()

write.csv2(maille_cepage, file = "data/geographie/correspondance_maille_cepage.csv", row.names = FALSE)

# Ajout de la réserve utile :----

ru_france <- vect("data/geographie/ru.shp")
classes_RU = c(25,75,125,175,225)
ru_france$RUmean <- sapply(as.factor(ru_france$classe), function(x) classes_RU[x])
ru_france$RUinf <- sapply(as.factor(ru_france$classe), function(x) classes_RU[x]-25)
ru_france$RUsup <- sapply(as.factor(ru_france$classe), function(x) classes_RU[x]+25)
ru_france <- ru_france[, c("RUmean", "RUinf", "RUsup")]

intersections <- terra::intersect(safran_grid, ru_france) # ici renvoie bien les polygones communes mais découpé par la grille (donc avec plus d'entité)
intersections$area <- expanse(intersections) # calc surface
# Limiter le nombre de NA du aux classes 9 : prendre la deuxième plus grande aire d'intersections :
intersections$area[is.na(intersections$RUmean)] <- 0 

# intersect * surf
intersections <- cbind(intersections[,names(intersections)[!names(intersections) %in% c("RUsup","RUmean","RUinf")]], 
                  do.call(cbind, lapply(c("RUsup", "RUmean", "RUinf"), function(x) intersections[[x]] * intersections$area)))

#  Aggregate by sum
agg <- aggregate(intersections[,c("ID", "area", "RUsup", "RUmean", "RUinf")], by="ID", fun=sum, na.rm=TRUE)
names(agg) <- gsub("agg_", "", names(agg))
#  Sum /surf
agg <- cbind(agg[,names(agg)[!names(agg) %in% c("RUsup","RUmean","RUinf")]], 
             do.call(cbind, lapply(c("RUsup", "RUmean", "RUinf"), function(x) agg[[x]] / agg$area)))

safran_grid[which(safran_grid$ID %in% agg$ID),]

agg$X <- safran_grid[which(safran_grid$ID %in% agg$ID),]$X
agg$Y <- safran_grid[which(safran_grid$ID %in% agg$ID),]$Y

safran_ru <- data.frame(ID_MAILLE = agg$ID, X_maille = agg$X/100, Y_maille = agg$Y/100, RU = as.numeric(agg$RUmean))

maille_commune <- maille_commune %>% left_join(safran_ru, by = c("X_maille", "Y_maille", "ID_MAILLE"))

# Ajout de l'altitude :----
alt_safran <- project(vect(meteoRIT::safran, crs = "EPSG:27572"), "EPSG:2154")
near <- nearest(y = alt_safran, x = safran_grid)
intersections <- terra::intersect(safran_grid, alt_safran)
intersections$area <- expanse(intersections) # calc surface
max_area_per_polygon <- intersections[order(intersections$ID, -intersections$area),] # tri selon la surface
max_area_per_polygon <- max_area_per_polygon[!duplicated(max_area_per_polygon$ID),] # garde le 1er (plus grande surface)
alt_safran <- data.frame(ID_MAILLE = max_area_per_polygon$ID, X_maille = max_area_per_polygon$X/100, 
                        Y_maille = max_area_per_polygon$Y/100, altitude = max_area_per_polygon$altitude)
maille_commune <- maille_commune %>% left_join(alt_safran, by = c("X_maille", "Y_maille", "ID_MAILLE"))

write.csv(maille_commune, file = "data/geographie/correspondance_maille_commune.csv", row.names = FALSE)
##

# Une fois les pixels safran identifiés, les sélectionner dans la BDD
# 

# Dataframe provisoire permettant de sélectionner les mailles de météo France
#XY_true <- data.frame(maille_ID = selected_grid_cells$ID, X = selected_grid_cells$X/100, Y = selected_grid_cells$Y/100)


# Fonction permettant de récuperer uniquement les mailles d'intérêt
select_maille <- function(XY_true, df){
  return(XY_true %>% left_join(df, by= c("X_maille" = "LAMBX", "Y_maille" = "LAMBY"), multiple = "all", 
                               relationship = "many-to-many"))
}


# Récupération données climat SAFRAN ----

# Get SAFRAN data 

# ! Attention, le lien de la dernière période change tout les mois : il faut mettre le mois précedent
# Exemple en avril 2025 : "2020-202503"

date_actuelle <- Sys.Date()
mois_prec <- as.Date(format(date_actuelle, "%Y-%m-01")) - 1 # récuperer le mois précedent
mois_prec <- format(mois_prec, "%Y%m")# Formater en "annéemois"

periods <- c("1980-1989", "1990-1999", "2000-2009", "2010-2019", paste0("2020-", mois_prec))
complet = data.frame()
# 

# 
# Increase time for downloading
options(timeout = 600)  # 600 seconds = 10 minutes
# 
start_time <- Sys.time()
# 
for(periode in periods){
  print(periode)
  if(periode == periods[5]) url <- paste0("https://object.files.data.gouv.fr/meteofrance/data/synchro_ftp/REF_CC/SIM/QUOT_SIM2_previous-",
                                              periode, ".csv.gz")
  else url <- paste0("https://object.files.data.gouv.fr/meteofrance/data/synchro_ftp/REF_CC/SIM/QUOT_SIM2_", periode, ".csv.gz")
  df <- fread(url, 
              select = c("LAMBX", "LAMBY", "DATE", "PRENEI_Q", "PRELIQ_Q", "TINF_H_Q", "TSUP_H_Q", "HU_Q", "SWI_Q"))
  df <- select_maille(maille_cepage, df)
  df <- df %>% mutate(X_maille = X_maille * 100, Y_maille = Y_maille * 100)
  complet <- rbind(complet, df)
  rm(df) ; gc()
}

# Enregistrer les données climatiques SAFRAN sous forme de Rdata
save(complet, file = "data/climat/climat_maille_cepage_jour_complet.RData") 


end_time <- Sys.time()

print(end_time -  start_time) # 14.33782 mins
