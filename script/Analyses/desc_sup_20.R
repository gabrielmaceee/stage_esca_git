library(dplyr)
library(ggplot2)
library(glue)

load("data/modelisation/observations.RData")
observations$sup_20 <- as.factor(as.numeric((observations$pourcentage_esca >= 20)))
table(observations[(observations$pourcentage_esca >= 20), "region_viticole"])

obs_20 <- observations[(observations$pourcentage_esca >= 20),]

for(colname in c("cepage", "region_viticole", "annee", "age_parcelle_estime")){
  my_plot <- ggplot(obs_20, aes(x = obs_20 %>% pull(colname))) + 
    geom_bar(stat = 'count') +  
    coord_flip() + 
    ggtitle(paste0("Répartition des modalités de : ", colname, ", (>=20%)")) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())
  plot(my_plot)
}

obs_n_20 <- observations[(observations$pourcentage_esca < 20),]
for(colname in c("cepage", "region_viticole", "annee", "age_parcelle_estime")){
  my_plot <- ggplot(obs_n_20, aes(x = obs_n_20 %>% pull(colname))) + 
    geom_bar(stat = 'count') +  
    coord_flip() + 
    ggtitle(paste0("Répartition des modalités de : ", colname, ", (<20%)")) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())
  plot(my_plot)
}

table((obs_20 %>% group_by(ID_MAILLE) %>% summarise(n = n()))$n  )
View((obs_20 %>% group_by(ID_MAILLE) %>% summarise(n = n())))
n_20 <- obs_20 %>% group_by(ID_MAILLE) %>% summarise(n = n())

n_tt <- observations %>% group_by(ID_MAILLE) %>% summarise(n = n())

n_20 <- n_20 %>%left_join(n_tt, by = "ID_MAILLE", suffix = c("_20",".tt"))
n_20$part <- n_20$n_20 / n_20$n.tt


for(colname in c("cepage", "region_viticole", "annee", "age_parcelle_estime")){
  n_tt <- observations %>% group_by_(colname) %>% summarise(n = n())
  n_20 <- obs_20 %>% group_by_(colname) %>% summarise(n = n())
  n_tt <- n_tt %>%left_join(n_20, by = colname, suffix = c(".tt", "_20"))
  n_tt$n_20[is.na(n_tt$n_20)] <- 0
  n_tt$part <- n_tt$n_20 / n_tt$n.tt
  
  my_plot <- ggplot(n_tt, aes_string(x = colname, y = "part")) + 
    geom_bar(stat = 'identity') +  
    coord_flip() + 
    geom_text(data = n_tt, aes_string(x = colname, y = n_tt$part +0.01, label = "n.tt"), 
              position = position_dodge(width = 0.75), size = 3) +  
    ggtitle(paste0(colname, " : part d'incidence >= 20%")) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())
  
  plot(my_plot)
}

