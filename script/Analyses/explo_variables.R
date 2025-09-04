# Jointure ESCA - phénologies ----
# Remise à zéro :
date_pheno <- date_floraison %>% left_join(date_debourrement, by = c("maille", "cepage", "year"), 
                                           relationship = "one-to-one")
# Enregister les simulations de phénologies de 2003 à 2023 :
save(date_pheno, file = "data/phenologie/phenologie_2003_2023.RData") 
esca <- read.csv("data/ESCA/esca_filtre.csv")
maille_commune <- read.csv("data/geographie/correspondance_maille_commune.csv", sep = ",")
# Jointure selon l'année, le cépage, et la maille :
esca <- esca %>% left_join(maille_commune[, c("INSEE_COM", "ID_MAILLE")], 
                           by = c("code_commune_INSEE" = "INSEE_COM"), 
                           relationship = "many-to-one", multiple = "all")
esca <- esca %>% left_join(date_pheno,
                           by = c("annee" = "year", "cepage" = "cepage", "ID_MAILLE" = "maille"),
                           relationship = "many-to-one")


## Corrélation incidences esca / jours juliens de phéno
cor(esca$pourcentage_esca, esca$debourrement_SUWE, use = "na.or.complete")
plot(esca$debourrement_SUWE, esca$pourcentage_esca)
cor(esca$pourcentage_esca, esca$flo_gfv, use = "na.or.complete")
plot(esca$flo_gfv, esca$pourcentage_esca)
cor(esca$pourcentage_esca, esca$flo_SUWE, use = "na.or.complete")
plot(esca$flo_SUWE, esca$pourcentage_esca)
# Pas de corrélation : peut être car effet du cépage
# -> glm 
library(lme4)
library(lmerTest)
##
rl <- lmer(pourcentage_esca ~ flo_gfv + (1|cepage), data = esca)
summary(rl) 
coef(rl) 
# Effet aléatoire du cépage juste sur l'intercept : pas d'effet de la floraison sur l'incidence d'esca
# p-value 0.130061 > 0.05.

##
rl <- lmer(pourcentage_esca ~ debourrement_SUWE + (1|cepage), data = esca)
summary(rl) 
coef(rl) 

# * intéraction entre le cépage et la phénologie sur l'incidence de l'esca :
rl <- lmerTest::lmer(pourcentage_esca ~ flo_gfv * cepage + (1|cepage), data = esca)
summary(rl) 
coef(rl) 

rl <- lmerTest::lmer(pourcentage_esca ~ flo_SUWE * cepage + (1|cepage), data = esca)
summary(rl) 
coef(rl) 

rl <- lmerTest::lmer(pourcentage_esca ~ debourrement_SUWE * cepage + (1|cepage), data = esca)
summary(rl) 
coef(rl) 
