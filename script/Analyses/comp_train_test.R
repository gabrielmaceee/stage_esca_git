library(dplyr)
load("data/modelisation/observations_parcelles.RData")


# On définit la proportion de données qu'on veut dans le jeu d'entraînement
train_ratio <- 0.8
obs_sort <- observations %>% arrange(annee)
train_indices <- obs_sort %>%
  group_by(identifiant_parcelle_analyse) %>%
  group_split() %>%
  lapply(function(group) {
    n <- nrow(group)
    if(n <= 5) seq_len(ceiling(n/2)) # ceiling : arrondi au supérieur, floor : arrondi à l'inférieur
    else seq_len(floor(n*train_ratio)) 
  })


# Construction du jeu d'entraînement et de test
# Associer les indices à chaque groupe
train <- bind_rows(
  Map(function(group, idx) group[idx, ], group_split(obs_sort, obs_sort$identifiant_parcelle_analyse), train_indices)
)

test <- anti_join(obs_sort, train, by = colnames(obs_sort))


moy_parcelle <- train %>% group_by(identifiant_parcelle_analyse) %>% summarise(moy_parcelle = mean(pourcentage_esca))
train <- train %>% left_join(moy_parcelle, by = "identifiant_parcelle_analyse")
train$moy_parcelle <- as.factor(round(train$moy_parcelle, 0))
test <- test %>% left_join(moy_parcelle, by = "identifiant_parcelle_analyse")
test$moy_parcelle <- as.factor(round(test$moy_parcelle, 0))

plot(density(train$pourcentage_esca), col = "red")
lines(density(test$pourcentage_esca), col = "blue")

mean(train$pourcentage_esca)
mean(test$pourcentage_esca)
mean(observations$pourcentage_esca)

median(train$pourcentage_esca)
median(test$pourcentage_esca)
median(observations$pourcentage_esca)

sd(train$pourcentage_esca)
sd(test$pourcentage_esca)
sd(observations$pourcentage_esca)

summary(train$pourcentage_esca)
summary(test$pourcentage_esca)

hist(train$annee, breaks = 21, xlab = "Année", main = "Compte des années jeu train")
train %>% 
  count(annee) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(annee, annee))) + 
  geom_bar(fill = "grey", stat = "identity") +
  labs(title = "Compte des années jeu train",
       x = "Nombre d'observations",
       y = "Année") + coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


test %>% 
  count(annee) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(annee, annee))) + 
  geom_bar(fill = "grey", stat = "identity") +
  labs(title = "Compte des années jeu test",
       x = "Nombre d'observations",
       y = "Année") + coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


observations %>% 
  count(annee) %>% arrange(desc(n)) %>%
  ggplot(aes(x = n, y = reorder(annee, annee))) + 
  geom_bar(fill = "grey", stat = "identity") +
  labs(title = "Compte des années jeu complet",
       x = "Nombre d'observations",
       y = "Année") + coord_flip() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


to_box <- data.frame(jeu = c(rep("train", dim(train)[1]), rep("test", dim(test)[1])),
                     esca = c(train$pourcentage_esca, test$pourcentage_esca))
boxplot(data = to_box, esca~jeu)

ggplot(to_box, aes(x=jeu, y=esca, fill=jeu)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set1")


observations %>% count(annee)
