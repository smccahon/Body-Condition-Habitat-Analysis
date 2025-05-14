#----------------------------------------------#
#      All Species Fat Habitat Analysis        #
#          Created 2025-04-11                  #
#          Modified 2025-05-13                 #
#----------------------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(trtools)
library(lme4)
library(car)
library(viridis)
library(lubridate)

options(digits = 3)

# read data
birds <- read.csv("Body_Condition_Habitat_Analysis_2025-03-31.csv")

# ...make new columns ----
# neonicotinoid detection column
birds$Detection <- ifelse(birds$OverallNeonic > 0, 
                          "Detection", "Non-detection")

# ...reorder and manipulate relevant factor variables ----
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Detection <- as.factor(birds$Detection)

birds$AgCategory <- factor(birds$AgCategory,
                           levels = c("Low", "Moderate", "High"))

birds$DominantCrop <- factor(birds$DominantCrop,
                             levels = c("Grassland",
                                        "Soybean",
                                        "Wheat",
                                        "Canola"))

birds$Event <- factor(birds$Event,
                      levels = c("Fall 2021",
                                 "Spring 2022",
                                 "Spring 2023",
                                 "Fall 2023"))

birds$Site <- as.factor(birds$Site)

birds$Permanence <- factor(birds$Permanence,
                           levels = c("Temporary", "Seasonal", "Semipermanent", 
                                      "Permanent"))

birds$MigStatus <- factor(birds$MigStatus,
                          levels = c("Resident", "Migratory"))

# standardize time to something more simple ------------------------------------
birds$Time <- strptime(birds$Time, format = "%H:%M")
birds$Time <- as.POSIXct(birds$Time, tz = "America/Chicago")
attributes(birds$Time)$tzone

# Calculate seconds since midnight (start of the day)
birds$seconds_since_midnight <- as.numeric(difftime(birds$Time, 
                                                    floor_date(birds$Time, "day"), 
                                                    units = "secs"))

birds$sin <- sin(2 * pi * birds$seconds_since_midnight / (24 * 3600))
birds$cos <- cos(2 * pi * birds$seconds_since_midnight / (24 * 3600))

birds$formatted_time <- format(as.POSIXct(birds$seconds_since_midnight, 
                                          origin = "1970-01-01", tz = "UTC"), 
                               "%H:%M")

# only include birds with fat
birds <- birds %>% 
  filter(!is.na(Fat))

birds$Fat <- factor(birds$Fat,
                    levels = c("0", "1", "2", "3", "4", "5", "6"))

birds <- birds %>%
  mutate(Fat.G = case_when(
    Fat == 0 | Fat == 1 ~ "Low",
    Fat == 2 | Fat == 3 | Fat == 4 | Fat == 5 | Fat == 6 ~ "High"       
  ))

birds$Fat.G <- factor(birds$Fat.G,
                      levels = c("Low", "High"))

# Subset data to omit spring (multicollinearity between fat and spring)
birds <- birds[birds$Event %in% c("Fall 2021", "Fall 2023"), ]

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# Test for Correlations--------------------------------------------------------- 

# subset data
sample <- birds[, c("PercentAg",
                       "Percent_Total_Veg",
                       "Event",
                       "Sex",
                       "Biomass",
                       "Diversity",
                       "Permanence",
                       "Percent_Exposed_Shoreline",
                       "Detection",
                       "seconds_since_midnight",
                       "Site",
                       "AgCategory",
                       "SPEI",
                       "Julian",
                       "DominantCrop",
                       "NearestCropDistance_m",
                       "Dist_Closest_Wetland_m",
                       "Max_Flock_Size",
                       "MigStatus"
)]

# convert categorical to numeric for correlation matrix
sample$Event <- as.numeric(sample$Event)
sample$Permanence <- as.numeric(sample$Permanence)
sample$AgCategory <- as.numeric(sample$AgCategory)
sample$DominantCrop <- as.numeric(sample$DominantCrop)
sample$Sex <- as.numeric(sample$Sex)
sample$Detection <- as.numeric(sample$Detection)
sample$Site <- as.numeric(sample$Site)
sample$MigStatus <- as.numeric(sample$MigStatus)

cor(sample)
# correlations > 0.6:
# % ag and ag category
# % ag and nearest crop distance
# % exposed shoreline and % veg
# Event & SPEI
# % exposed shoreline and permanence
# permanence and SPEI
# permanence and julian
# permanence and exposed shoreline
# exposed shoreline and SPEI
# exposed shoreline and Julian
# SPEI and julian


# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- glmer(Fat.G ~ PercentAg + (1 | Species), family = "binomial", data = birds.cs)
m2 <- glmer(Fat.G ~ DominantCrop + (1 | Species), family = "binomial", data = birds.cs)
m3 <- glmer(Fat.G ~ AgCategory + (1 | Species), family = "binomial", data = birds.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m3)

# Ag Category by far the best but perfect separation, use % ag

# temporal --> yes include interaction
m1 <- glmer(Fat.G ~ Julian + (1 | Species), family = "binomial", data = birds.cs)
m2 <- glmer(Fat.G ~ Julian * MigStatus + (1 | Species), family = "binomial", data = birds.cs)
m3 <- glmer(Fat.G ~ MigStatus + (1 | Species), family = "binomial", data = birds.cs)
m4 <- glmer(Fat.G ~ Julian + MigStatus + (1 | Species), family = "binomial", data = birds.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m2)

# transformation for time needed? no
plot(birds$seconds_since_midnight, birds$Fat.G)

# interactions between season and capture time? yes, keep in models
m1 <- glmer(Fat.G ~ Event + seconds_since_midnight + (1 | Species), family = "binomial", 
          data = birds.cs)
m2 <- glmer(Fat.G ~ Event * seconds_since_midnight + (1 | Species), family = "binomial", 
          data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m2) # interaction significant, keep in model

# is species needed as a fixed effect?
m <- glmer(Fat.G ~ (1 | Species), family = "binomial", data = birds.cs)

# Model Selection with interactions (Stage 1) ----------------------------------

# agriculture
m1 <- glmer(Fat.G ~ PercentAg + (1 | Species), family = "binomial", 
            data = birds.cs)

# vegetation
m2 <- glmer(Fat.G ~ Percent_Total_Veg + (1 | Species), family = "binomial", 
            data = birds.cs)

# habitat
m3 <- glmer(Fat.G ~ Permanence + (1 | Species), family = "binomial", 
              data = birds.cs)
m4 <- glmer(Fat.G ~ Percent_Exposed_Shoreline + (1 | Species), family = "binomial", 
            data = birds.cs)
m5 <- glmer(Fat.G ~ Dist_Closest_Wetland_m + (1 | Species), family = "binomial", 
            data = birds.cs)

# weather
m6 <- glmer(Fat.G ~ SPEI + (1 | Species), family = "binomial", data = birds.cs)

# life history
m7 <- glmer(Fat.G ~ Sex + (1 | Species), family = "binomial", data = birds.cs)
m8 <- glmer(Fat.G ~ MigStatus + (1 | Species), family = "binomial", data = birds.cs)

# temporal
m9 <- glmer(Fat.G ~ Julian + (1 | Species), family = "binomial", 
            data = birds.cs)
m10 <- glmer(Fat.G ~ seconds_since_midnight + (1 | Species), family = "binomial", 
            data = birds.cs)
m11 <- glmer(Fat.G ~ Event + (1 | Species), family = "binomial", data = birds.cs)

# interactions
m12 <- glmer(Fat.G ~ Event * seconds_since_midnight + (1 | Species), 
             family = "binomial", data = birds.cs)

m13 <- glmer(Fat.G ~ Julian * MigStatus + (1 | Species), 
             family = "binomial", data = birds.cs)

# flock
m14 <- glmer(Fat.G ~ Max_Flock_Size + (1 | Species), family = "binomial", 
             data = birds.cs)

# null
m15 <- glmer(Fat.G ~ (1 | Species), family = "binomial", data = birds.cs)


model_names <- paste0("m", 1:15)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results (informative parameters across models)
confint(m13) # julian * migratory status important (migrants lose fat with time)

confint(m9) # julian date important

summary(m6)
confint(m6) # SPEI significant (fat is higher in wetter years)...correlated with date

summary(m7) 
confint(m7) # sex not important

summary(m12)
confint(m12) # event * time important

summary(m3)
confint(m3) # Fat is lower in permanent wetlands compared to seasonal


# stage 2 doesn't make sense because all variables are highly correlated




# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- glmer(Fat.G ~ Biomass + (1 | Species), data = birds.sub, family = "binomial")

summary(m)
confint(m) # no effect of biomass

m <- glmer(Fat.G ~ Diversity + (1 | Species), data = birds.sub, family = "binomial")

summary(m)
confint(m) # no effect of diversity










