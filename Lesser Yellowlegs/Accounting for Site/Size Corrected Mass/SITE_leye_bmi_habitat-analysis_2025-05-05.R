#-----------------------------------------#
#  Lesser Yellowlegs BMI Habitat Analysis #
#          Linear regression              #
#          Created 2025-04-10             #
#         Modified 2025-05-05             #
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

# standardize time to something more simple ------------------------------------
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

# ---------------------------------------------------------------------------- #

# Perform PCA: Females ####
# Peig and Green 2009 & Bajracharya 2022

# Subset the dataset to only include females
leye_female <- subset(leye, Sex == "Female")

# Subset the dataset to only include tarsus, wing length, and bill length
leye_female <- leye_female[, c("Wing", "Culmen", "DiagTarsus")]

# Run PCA on the subsetted data
pca_result <- prcomp(leye_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_female <- pca_result$x[, 1]

# Regress PC1 against body mass
m_female <- lm(PC1_female ~ Mass, data = subset(leye, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.f <- resid(m_female)

# Add size-corrected mass back into the full datasest
leye$sc.mass <- NA  # Initialize with NA values
leye[leye$Sex == "Female", "sc.mass"] <- sc.mass.f


# Perform PCA: Males ####

# Subset the dataset to only include males
leye_male <- subset(leye, Sex == "Male")

# Subset the dataset to only include tarsus, wing length, and bill length
leye_male <- leye_male[, c("Wing", "Culmen", "DiagTarsus")]

# Run PCA on the subsetted data
pca_result <- prcomp(leye_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_male <- pca_result$x[, 1]

# Regress PC1 against body mass
m_male <- lm(PC1_male ~ Mass, data = subset(leye, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.m <- resid(m_male)

# Add size-corrected mass back into the full datasest
leye[leye$Sex == "Male", "sc.mass"] <- sc.mass.m

# Add size-corrected mass back into the full datasest
leye[leye$Sex == "Male", "sc.mass"] <- sc.mass.m

# Add size-corrected mass back into the full datasest
leye[leye$Sex == "Female", "sc.mass"] <- sc.mass.f



# ...standardize data except for response ----
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("sc.mass"), scale))

# Test for Correlations--------------------------------------------------------- 

# subset data
sample <- leye[, c("PercentAg",
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
                   "Age",
                   "DominantCrop",
                   "NearestCropDistance_m",
                   "Dist_Closest_Wetland_m",
                   "Max_Flock_Size"
)]

# convert categorical to numeric for correlation matrix
sample$Event <- as.numeric(sample$Event)
sample$Permanence <- as.numeric(sample$Permanence)
sample$AgCategory <- as.numeric(sample$AgCategory)
sample$DominantCrop <- as.numeric(sample$DominantCrop)
sample$Sex <- as.numeric(sample$Sex)
sample$Detection <- as.numeric(sample$Detection)
sample$Site <- as.numeric(sample$Site)
sample$Age <- as.numeric(sample$Age)

cor(sample)
# correlations > 0.6:
# % ag & ag category
# % ag & dominant crop
# % ag & nearest crop distance
# % ag and max flock size
# % total veg & dominant crop
# flock size and dominant crop
# permanence and SPEI
# shoreline and SPEI
# detection & julian

# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lm(sc.mass ~ PercentAg, data = leye.cs)
m2 <- lm(sc.mass ~ DominantCrop, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a better predictor

# Does time need a transformation? no
plot(leye.cs$sc.mass, leye.cs$seconds_since_midnight)

# temporal
m1 <- lm(sc.mass ~ Julian * Event, data = leye.cs)
m2 <- lm(sc.mass ~ Julian + Event, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model with interaction is MUCH better

# is time * event necessary? no
m1 <- lm(sc.mass ~ seconds_since_midnight * Event, data = leye.cs)
m2 <- lm(sc.mass ~ seconds_since_midnight + Event, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lm(sc.mass ~ PercentAg, data = leye.cs)

# vegetation
m2 <- lm(sc.mass ~ Percent_Total_Veg, data = leye.cs)

# habitat
m3 <- lm(sc.mass ~ Permanence, data = leye.cs)
m4 <- lm(sc.mass ~ Percent_Exposed_Shoreline, data = leye.cs)
m5 <- lm(sc.mass ~ Dist_Closest_Wetland_m, data = leye.cs)

# weather
m6 <- lm(sc.mass ~ SPEI, data = leye.cs)

# life history
m7 <- lm(sc.mass ~ Age, data = leye.cs)
m8 <- lm(sc.mass ~ Sex, data = leye.cs)

# temporal
m9 <- lm(sc.mass ~ Julian, data = leye.cs)
m10 <- lm(sc.mass ~ seconds_since_midnight, data = leye.cs)
m11 <- lm(sc.mass ~ Event, data = leye.cs)

# flock
m12 <- lm(sc.mass ~ Max_Flock_Size, data = leye.cs)

# null
m13 <- lm(sc.mass ~ 1, data = leye.cs)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results --> permanence top model but not significant (next is null)
summary(m3)
confint(m3)

# Random effect necessary? ----

# exclude sites with only 1 bird
leye.m <- leye.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()

# run reduced dataset with site as a random effect ----
# agriculture
m1 <- lmer(sc.mass ~ PercentAg + (1|Site), REML = FALSE, data = leye.m)

# vegetation
m2 <- lmer(sc.mass ~ Percent_Total_Veg + (1|Site), REML = FALSE, data = leye.m)

# habitat
m3 <- lmer(sc.mass ~ Permanence + (1|Site), REML = FALSE, data = leye.m)
m4 <- lmer(sc.mass ~ Percent_Exposed_Shoreline + (1|Site), REML = FALSE, data = leye.m)
m5 <- lmer(sc.mass ~ Dist_Closest_Wetland_m + (1|Site), REML = FALSE, data = leye.m)

# weather
m6 <- lmer(sc.mass ~ SPEI + (1|Site), REML = FALSE, data = leye.m)

# life history
m7 <- lmer(sc.mass ~ Age + (1|Site), REML = FALSE, data = leye.m)
m8 <- lmer(sc.mass ~ Sex + (1|Site), REML = FALSE, data = leye.m)

# temporal
m9 <- lmer(sc.mass ~ Julian + (1|Site), REML = FALSE, data = leye.m)
m10 <- lmer(sc.mass ~ seconds_since_midnight + (1|Site), REML = FALSE, data = leye.m)
m11 <- lmer(sc.mass ~ Event + (1|Site), REML = FALSE, data = leye.m)

# flock
m12 <- lmer(sc.mass ~ Max_Flock_Size +(1|Site), REML = FALSE, data = leye.m)

# null
m13 <- lmer(sc.mass ~ 1 + (1|Site), data = leye.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# informative parameters --- none and model is overfitting
confint(m3)

# rerun analysis without random effect and with reduced dataset ----
# agriculture
m1 <- lm(sc.mass ~ PercentAg, data = leye.m)

# vegetation
m2 <- lm(sc.mass ~ Percent_Total_Veg, data = leye.m)

# habitat
m3 <- lm(sc.mass ~ Permanence, data = leye.m)
m4 <- lm(sc.mass ~ Percent_Exposed_Shoreline, data = leye.m)
m5 <- lm(sc.mass ~ Dist_Closest_Wetland_m, data = leye.m)

# weather
m6 <- lm(sc.mass ~ SPEI, data = leye.m)

# life history
m7 <- lm(sc.mass ~ Age, data = leye.m)
m8 <- lm(sc.mass ~ Sex, data = leye.m)

# temporal
m9 <- lm(sc.mass ~ Julian, data = leye.m)
m10 <- lm(sc.mass ~ seconds_since_midnight, data = leye.m)
m11 <- lm(sc.mass ~ Event, data = leye.m)

# flock
m12 <- lm(sc.mass ~ Max_Flock_Size, data = leye.m)

# null
m13 <- lm(sc.mass ~ 1, data = leye.m)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m1 <- lm(sc.mass ~ Permanence, data = leye.m)
m2 <- lmer(sc.mass ~ Permanence + (1|Site), data = leye.m, REML = FALSE)

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
