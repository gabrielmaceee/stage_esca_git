# Prédiction de l'évolution de l’incidence de l’esca en France avec le changement climatique

## Getting started

## Outils requis

R et Python

## Origine des données

 - Observations d'esca :  Plateforme d'Épidémiosurveillance en Santé Végétale (ESV)
 - Observations météo : Météo France : SICLIMA : SAFRAN
 - Prévisions climatique : 19 modèles de climat issues de CMIP6 (GIEC)
 

## Organisation des dossiers

### data
 - climat : Données d'observations et de prévisions climatiques
 - ESCA : Données d'observations d'esca
 - geographie : Données associé aux mailles, régions, communes ...
 - modelisation : Données à utiliser pour les modélisations
 - vitadapt : Données d'observations de la parcelle vitadapt (météo et esca)


### graphs
 - carto : Cartographie d'une partie des indices éco-climatiques, et de l'incidence d'esca.
 - distribution_var_selectionnees : Distribution des variables sélectionnées pour le glm final.
 - lien_var_selectionnees_esca : Lien entre l'incidence d'esca et les variables sélectionnées pour le glm final.
 - modelisation : Pas fini : Graphiques intéractifs des valeurs de Shapley.
 - stats_descriptives : Répartition de l'incidence d'esca et des parcelles.



### modeles
Sauvegarde des modèles sélectionnés.

### R_functions
Scripts contenant des fonctions utiles pour la simulation des phénologies ("fonctions_pheno.R"), le calcul des 
valeurs de Shapley en R("shap.R"), ainsi que pour simuler les incidences d'esca futures ("functions_projections.R")

### script
 - Analyses : Contient tous les scripts utiles à l'analyse de nos données et la modélisation, organisés par 
 ordre temporel.
 - Bilan hydrique : Contient les scripts utiles à la simulation du bilan hydrique (uniquement la partie diffusable).
 - Projections : Contient les scripts utiles à la simulations des projections futures.

