#-----------------------------------------#
#  Lesser Yellowlegs Fat Habitat Analysis #
#          Linear regression              #
#          Created 2025-04-30             #
#         Modified 2025-04-30             #
#-----------------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(trtools)
library(lme4)
library(car)
library(viridis)
library(lubridate)
library(nlme)

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
                                 "Fall 2023"))

birds$Site <- as.factor(birds$Site)

birds$Age <- factor(birds$Age,
                    levels = c("Juvenile", "Adult"))

# combine temporary and seasonal wetlands
birds$Permanence[birds$Permanence %in% c("Temporary", "Seasonal")] <- 
  "Temporary/Seasonal"

birds$Permanence <- factor(birds$Permanence,
                           levels = c("Temporary/Seasonal", "Semipermanent", 
                                      "Permanent"))

# ...subset and filter data ----
leye <- subset(birds, Species %in% c("Lesser Yellowlegs"))

# ...standardize time  ------------------------------------------------------------
leye$Time <- strptime(leye$Time, format = "%H:%M")
leye$Time <- as.POSIXct(leye$Time, tz = "America/Chicago")
attributes(leye$Time)$tzone

# Calculate seconds since midnight (start of the day)
leye$seconds_since_midnight <- as.numeric(difftime(leye$Time, 
                                                   floor_date(leye$Time, "day"), 
                                                   units = "secs"))

leye$sin <- sin(2 * pi * leye$seconds_since_midnight / (24 * 3600))
leye$cos <- cos(2 * pi * leye$seconds_since_midnight / (24 * 3600))

leye$formatted_time <- format(as.POSIXct(leye$seconds_since_midnight, 
                                         origin = "1970-01-01", tz = "UTC"), 
                              "%H:%M")

# ...standardize data except for response ----
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("Fat"), scale))

# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lmer(Fat ~ PercentAg + (1|Site), REML = FALSE, data = leye.cs)

# vegetation
m2 <- lmer(Fat ~ Percent_Total_Veg + (1|Site), REML = FALSE, data = leye.cs)

# habitat
m3 <- lmer(Fat ~ Permanence + (1|Site), REML = FALSE, data = leye.cs)
m4 <- lmer(Fat ~ Percent_Exposed_Shoreline + (1|Site), REML = FALSE, data = leye.cs)
m5 <- lmer(Fat ~ Dist_Closest_Wetland_m + (1|Site), REML = FALSE, data = leye.cs)

# weather
m6 <- lmer(Fat ~ SPEI + (1|Site), REML = FALSE, data = leye.cs)

# life history
m7 <- lmer(Fat ~ Age + (1|Site), REML = FALSE, data = leye.cs)
m8 <- lmer(Fat ~ Sex + (1|Site), REML = FALSE, data = leye.cs)

# temporal
m9 <- lmer(Fat ~ Julian + (1|Site), REML = FALSE, data = leye.cs)
m10 <- lmer(Fat ~ seconds_since_midnight + (1|Site), REML = FALSE, data = leye.cs)
m11 <- lmer(Fat ~ Event + (1|Site), REML = FALSE, data = leye.cs)
m12 <- lmer(Fat ~ Event * seconds_since_midnight + (1|Site), REML = FALSE,
          data = leye.cs)

# flock
m13 <- lmer(Fat ~ Max_Flock_Size + (1|Site), REML = FALSE, data = leye.cs)

# null
m14 <- lmer(Fat ~ 1 + (1|Site), REML = FALSE, data = leye.cs)


model_names <- paste0("m", 1:14)

models <- mget(model_names)

aictab(models, modnames = model_names)

# informative parameters:
# julian, age

# Add Julian as informed null (Report in Thesis as Stage 1) ----

# agriculture
m1 <- lmer(Fat ~ PercentAg + Julian + (1|Site), REML = FALSE, data = leye.cs)

# vegetation
m2 <- lmer(Fat ~ Percent_Total_Veg + Julian + (1|Site), REML = FALSE, data = leye.cs)

# habitat
m3 <- lmer(Fat ~ Permanence + Julian + (1|Site), REML = FALSE, data = leye.cs)
m4 <- lmer(Fat ~ Percent_Exposed_Shoreline + Julian + (1|Site), REML = FALSE, data = leye.cs)
m5 <- lmer(Fat ~ Dist_Closest_Wetland_m + Julian + (1|Site), REML = FALSE, data = leye.cs)

# weather
m6 <- lmer(Fat ~ SPEI + Julian + (1|Site), REML = FALSE, data = leye.cs)

# life history
m7 <- lmer(Fat ~ Age + Julian + (1|Site), REML = FALSE, data = leye.cs)
m8 <- lmer(Fat ~ Sex + Julian + (1|Site), REML = FALSE, data = leye.cs)

# temporal
m9 <- lmer(Fat ~ Julian + (1|Site), REML = FALSE, data = leye.cs)
m10 <- lmer(Fat ~ seconds_since_midnight + Julian + (1|Site), REML = FALSE, data = leye.cs)
m11 <- lmer(Fat ~ Event + Julian + (1|Site), REML = FALSE, data = leye.cs)
m12 <- lmer(Fat ~ Event * seconds_since_midnight + Julian + (1|Site), REML = FALSE,
          data = leye.cs)

# flock
m13 <- lmer(Fat ~ Max_Flock_Size + Julian + (1|Site), REML = FALSE, data = leye.cs)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m7)
confint(m7)

summary(m6)
confint(m6)

# results change with random effect... age and Julian is top model, SPEI not significant

# exclude sites with only 1 bird
leye.m <- leye.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()

# run reduced dataset with site as a random effect ----
# agriculture
m1 <- lmer(Fat ~ PercentAg + Julian + (1|Site), REML = FALSE, data = leye.m)

# vegetation
m2 <- lmer(Fat ~ Percent_Total_Veg + Julian + (1|Site), REML = FALSE, data = leye.m)

# habitat
m3 <- lmer(Fat ~ Permanence + Julian + (1|Site), REML = FALSE, data = leye.m)
m4 <- lmer(Fat ~ Percent_Exposed_Shoreline + Julian + (1|Site), REML = FALSE, data = leye.m)
m5 <- lmer(Fat ~ Dist_Closest_Wetland_m + Julian + (1|Site), REML = FALSE, data = leye.m)

# weather
m6 <- lmer(Fat ~ SPEI + Julian + (1|Site), REML = FALSE, data = leye.m)

# life history
m7 <- lmer(Fat ~ Age + Julian + (1|Site), REML = FALSE, data = leye.m)
m8 <- lmer(Fat ~ Sex + Julian + (1|Site), REML = FALSE, data = leye.m)

# temporal
m9 <- lmer(Fat ~ Julian + (1|Site), REML = FALSE, data = leye.m)
m10 <- lmer(Fat ~ seconds_since_midnight + Julian + (1|Site), REML = FALSE, data = leye.m)
m11 <- lmer(Fat ~ Event + Julian + (1|Site), REML = FALSE, data = leye.m)
m12 <- lmer(Fat ~ Event * seconds_since_midnight + Julian + (1|Site), REML = FALSE,
            data = leye.m)

# flock
m13 <- lmer(Fat ~ Max_Flock_Size + Julian + (1|Site), REML = FALSE, data = leye.m)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# informative parameters
# age and Julian (top model and only one within 2 delta AICc)
# SPEI + Julian now significant (but not within 2 delta AICc)

# rerun reduced dataset without site as a random effect ----

# agriculture
m1 <- lm(Fat ~ PercentAg + Julian, data = leye.m)

# vegetation
m2 <- lm(Fat ~ Percent_Total_Veg + Julian, data = leye.m)

# habitat
m3 <- lm(Fat ~ Permanence + Julian, data = leye.m)
m4 <- lm(Fat ~ Percent_Exposed_Shoreline + Julian, data = leye.m)
m5 <- lm(Fat ~ Dist_Closest_Wetland_m + Julian, data = leye.m)

# weather
m6 <- lm(Fat ~ SPEI + Julian, data = leye.m)

# life history
m7 <- lm(Fat ~ Age + Julian, data = leye.m)
m8 <- lm(Fat ~ Sex + Julian, data = leye.m)

# temporal
m9 <- lm(Fat ~ Julian, data = leye.m)
m10 <- lm(Fat ~ seconds_since_midnight + Julian, data = leye.m)
m11 <- lm(Fat ~ Event + Julian, data = leye.m)
m12 <- lm(Fat ~ Event * seconds_since_midnight + Julian,
          data = leye.m)

# flock
m13 <- lm(Fat ~ Max_Flock_Size + Julian, data = leye.m)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# SPEI + Julian top model and only one within 2 delta AICc
m1 <- lm(Fat ~ SPEI + Julian, data = leye.m)
m2 <- lmer(Fat ~ SPEI + Julian + (1|Site), data = leye.m, REML = FALSE)

# Is there enough support to include the random effect? No ----

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
birds.sub <- subset(leye.cs, !is.na(Biomass))
birds.sub.or <- subset(birds.sub, Biomass < 3)

# fat (outlier needed to be removed)
m <- lm(Fat ~ Biomass, data = birds.sub.or)

summary(m)
confint(m) # negative effect of biomass on fat...after removing outlier, no effect

m <- lm(Fat ~ Diversity, data = birds.sub.or)

summary(m)
confint(m) # no effect of diversity

# mass 
m <- lm(Mass ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass on mass

m <- lm(Mass ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of biomass on mass

# pectoral muscle 
birds.sub.p <- subset(birds.sub, !is.na(PecSizeBest))

m <- lm(PecSizeBest ~ Biomass, data = birds.sub.p)

summary(m)
confint(m) # no effect of biomass on mass

m <- lm(PecSizeBest ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of biomass on mass

