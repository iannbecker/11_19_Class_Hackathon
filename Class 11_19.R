##############################
#
# Class Week 12
# 11/19//2024
# Ian Becker
#
##############################

library(rgbif)
library(dismo)
library(dplyr)
library(terra)
library(ggplot2)

#### Loading in gbif data

Species.Name <- "Vireo atricapilla"

Species.key <- name_backbone(name=Species.Name)$speciesKey
locations <- occ_search(taxonKey=Species.key,
                        hasCoordinate = TRUE,
                        limit=5000)

locations <- locations[["data"]]

#### Actual data 

vireo_data <- read.csv("vireo_data_sdm.csv")
env_data_current <- rast("env_current.grd")
env_data_forecast <- rast("env_forecast.grd")

plot(env_data_current$precip)

env_data_forecast

coords <- vireo_data %>%
  select(lon, lat)
extracted_values <- extract(env_data_current, coords)
vireo_data <- bind_cols(vireo_data, extracted_values)

#### Building SDM

logistic <- glm(present ~ tmin + precip, 
                family = binomial(link = "logit"),
                data = vireo_data)

summary(logistic)

predictions <- terra::predict(env_data_current, logistic, type = "response")
plot(predictions, main = "Pred Current", col=rev(heat.colors(25)), ext = extent(-140, -50, 20, 60))

predictions2 <- terra::predict(env_data_forecast, logistic, type = "response")
plot(predictions2, main = "Pred Future", col=rev(heat.colors(25)), ext = extent(-140, -50, 20, 60))

#### Hackathon 

install.packages("cubing")

library(cubing)

aCube <- getCubieCube("Superflip")

plot(aCube)

is.solvable(aCube, split = TRUE)

cubing::solver(aCube)

solution <- solver(aCube)

