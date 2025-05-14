#----------------------------------------------#
# All Species Pectoral Muscle Habitat Analysis #
#          Created 2025-04-11                  #
#         Modified 2025-04-11                  #
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

# Only include birds with pectoral muscles 
birds <- birds %>%
  filter(!is.na(PecSizeBest))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data except for response
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("PecSizeBest"), scale))

# Only include sites with at least three individuals
birds.s <- birds %>% 
  group_by(Site) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data except for response
birds.s.cs <- birds.s %>%
  mutate(across(where(is.numeric) & !matches("PecSizeBest"), scale))


# Test for Correlations--------------------------------------------------------- 

# subset data
sample <- birds.cs[, c("PercentAg",
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
# % ag and dominant crop
# % veg and % exposed shoreline
# event and SPEI
# event and julian
# SPEI and julian

# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lmer(PecSizeBest ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PecSizeBest ~ DominantCrop + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PecSizeBest ~ AgCategory + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)

# Ag Category is best but not significant so doesn't really matter
# I reran analysis with ag category and results did not change from % ag

# temporal
m1 <- lmer(PecSizeBest ~ Julian + MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(PecSizeBest ~ Julian * MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m3 <- lmer(PecSizeBest ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(PecSizeBest ~ Julian + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m2)

# model with interaction is MUCH more informative, keep in models

# transformation for time needed? no
plot(birds$seconds_since_midnight, birds$PecSizeBest)

# interactions between season and capture time? yes, keep in models
m1 <- lmer(PecSizeBest ~ Event + seconds_since_midnight + (1|Species), 
           data = birds.cs, 
           REML = FALSE)
m2 <- lmer(PecSizeBest ~ Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, 
           REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m2)

# Model Selection with interactions (Stage 1) ----------------------------------

# agriculture
m1 <- lmer(PecSizeBest ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(PecSizeBest ~ Percent_Total_Veg + (1|Species), 
           data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(PecSizeBest ~ Permanence + (1|Species), 
           data = birds.cs, REML = FALSE)
m4 <- lmer(PecSizeBest ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(PecSizeBest ~ Dist_Closest_Wetland_m + (1|Species), 
           data = birds.cs, REML = FALSE)

# weather
m6 <- lmer(PecSizeBest ~ SPEI + (1|Species), 
           data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(PecSizeBest ~ Sex + (1|Species), 
           data = birds.cs, REML = FALSE)
m8 <- lmer(PecSizeBest ~ MigStatus + (1|Species), 
           data = birds.cs, REML = FALSE)
m9 <- lmer(PecSizeBest ~ MigStatus * Julian + (1|Species), 
           data = birds.cs, REML = FALSE)

# temporal
m10 <- lmer(PecSizeBest ~ Julian + (1|Species), 
           data = birds.cs, REML = FALSE)
m11 <- lmer(PecSizeBest ~ Event * seconds_since_midnight + (1|Species), data = birds.cs,
            REML = FALSE)
m12 <- lmer(PecSizeBest ~ Event + (1|Species), data = birds.cs,
            REML = FALSE)
m13 <- lmer(PecSizeBest ~ seconds_since_midnight + (1|Species), data = birds.cs,
            REML = FALSE)

# flock
m14 <- lmer(PecSizeBest ~ Max_Flock_Size + (1|Species), data = birds.cs, REML = FALSE)

# null
m15 <- lmer(PecSizeBest ~ (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:15)

models <- mget(model_names)

aictab(models, modnames = model_names)

# event * time, julian * migratory status, max flock size
# can't have Julian*MigStatus as informed null due to correlations

# informative parameters to move onto stage 2 ----------------------------------
# event * time, julian * migratory status

plot(birds$Max_Flock_Size, birds$PecSizeBest)

# Model Selection with Site as Random Effect -----------------------------------

# agriculture
m1 <- lmer(PecSizeBest ~ PercentAg + (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# vegetation
m2 <- lmer(PecSizeBest ~ Percent_Total_Veg + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)

# habitat
m3 <- lmer(PecSizeBest ~ Permanence + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)
m4 <- lmer(PecSizeBest ~ Percent_Exposed_Shoreline + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)
m5 <- lmer(PecSizeBest ~ Dist_Closest_Wetland_m + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)

# weather
m6 <- lmer(PecSizeBest ~ SPEI + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)

# life history
m7 <- lmer(PecSizeBest ~ Sex + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)
m8 <- lmer(PecSizeBest ~ MigStatus + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)
m9 <- lmer(PecSizeBest ~ MigStatus * Julian + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)

# temporal
m10 <- lmer(PecSizeBest ~ Julian + (1|Species) + (1|Site), 
            data = birds.s.cs, REML = FALSE)
m11 <- lmer(PecSizeBest ~ Event * seconds_since_midnight + (1|Species) + (1|Site), data = birds.s.cs,
            REML = FALSE)
m12 <- lmer(PecSizeBest ~ Event + (1|Species) + (1|Site), data = birds.s.cs,
            REML = FALSE)
m13 <- lmer(PecSizeBest ~ seconds_since_midnight + (1|Species) + (1|Site), data = birds.s.cs,
            REML = FALSE)

# flock
m14 <- lmer(PecSizeBest ~ Max_Flock_Size + (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# null
m15 <- lmer(PecSizeBest ~ (1 | Species) + (1|Site), data = birds.s.cs, REML = FALSE)

model_names <- paste0("m", 1:15)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m1 <- lmer(PecSizeBest ~ MigStatus * Julian + (1|Species), 
           data = birds.s.cs, REML = FALSE)

m2 <- lmer(PecSizeBest ~ MigStatus * Julian + (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)

# Calculate AICc values for both models
AICc_m1 <- AICc(m1)
AICc_m2 <- AICc(m2)

# Step 1: Calculate delta AICc (ΔAICc)
delta_AICc_m1 <- AICc_m1 - min(AICc_m1, AICc_m2)
delta_AICc_m2 <- AICc_m2 - min(AICc_m1, AICc_m2)

# Step 2: Calculate the exponentiated delta AICc values
exp_delta_m1 <- exp(-delta_AICc_m1 / 2)
exp_delta_m2 <- exp(-delta_AICc_m2 / 2)

# Step 3: Sum of all exponentiated delta AICc values
sum_exp_delta <- exp_delta_m1 + exp_delta_m2

# Step 4: Calculate model weights
weight_m1 <- exp_delta_m1 / sum_exp_delta
weight_m2 <- exp_delta_m2 / sum_exp_delta

# Step 5: Create a summary data frame
model_summary <- data.frame(
  Model = c("m1", "m2"),
  AICc = c(AICc_m1, AICc_m2),
  Delta_AICc = c(delta_AICc_m1, delta_AICc_m2),
  Weight = c(weight_m1, weight_m2)
)

# Step 6: Display the summary table
print(model_summary)



# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- lmer(PecSizeBest ~ Biomass + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lmer(PecSizeBest ~ Diversity + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of diversity


