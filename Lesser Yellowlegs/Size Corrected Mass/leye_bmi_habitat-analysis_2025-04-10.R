#-----------------------------------------#
#  Lesser Yellowlegs BMI Habitat Analysis #
#          Linear regression              #
#          Created 2025-04-10             #
#         Modified 2025-04-10             #
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

#standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Perform PCA: Females ####
# Peig and Green 2009 & Bajracharya 2022

# Subset the dataset to only include females
leye_female <- subset(leye.cs, Sex == "Female")

# Subset the dataset to only include tarsus, wing length, and bill length
leye_female <- leye_female[, c("Wing", "Culmen", "DiagTarsus")]

# Run PCA on the subsetted data
pca_result <- prcomp(leye_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_female <- pca_result$x[, 1]

# Regress PC1 against body mass
m_female <- lm(PC1_female ~ Mass, data = subset(leye.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.f <- resid(m_female)

# Add size-corrected mass back into the full datasest
leye.cs$sc.mass <- NA  # Initialize with NA values
leye.cs[leye.cs$Sex == "Female", "sc.mass"] <- sc.mass.f


# Perform PCA: Males ####

# Subset the dataset to only include males
leye_male <- subset(leye.cs, Sex == "Male")

# Subset the dataset to only include tarsus, wing length, and bill length
leye_male <- leye_male[, c("Wing", "Culmen", "DiagTarsus")]

# Run PCA on the subsetted data
pca_result <- prcomp(leye_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_male <- pca_result$x[, 1]

# Regress PC1 against body mass
m_male <- lm(PC1_male ~ Mass, data = subset(leye.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.m <- resid(m_male)

# Add size-corrected mass back into the full datasest
leye.cs[leye.cs$Sex == "Male", "sc.mass"] <- sc.mass.m

# Add size-corrected mass back into the full datasest
leye[leye.cs$Sex == "Male", "sc.mass"] <- sc.mass.m

# Add size-corrected mass back into the full datasest
leye[leye.cs$Sex == "Female", "sc.mass"] <- sc.mass.f

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
m1 <- lm(Uric ~ PercentAg, data = leye)
m2 <- lm(Uric ~ DominantCrop, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a better predictor

# Does time need a transformation? no
plot(leye$Uric, leye$seconds_since_midnight)

# temporal
m1 <- lm(Uric ~ Julian * Event, data = leye)
m2 <- lm(Uric ~ Julian + Event, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model without interaction is better

# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lm(Uric ~ PercentAg, data = leye.cs)

# vegetation
m2 <- lm(Uric ~ Percent_Total_Veg, data = leye.cs)

# habitat
m3 <- lm(Uric ~ Permanence, data = leye.cs)
m4 <- lm(Uric ~ Percent_Exposed_Shoreline, data = leye.cs)
m5 <- lm(Uric ~ Dist_Closest_Wetland_m, data = leye.cs)

# weather
m6 <- lm(Uric ~ SPEI, data = leye.cs)

# life history
m7 <- lm(Uric ~ Age, data = leye.cs)
m8 <- lm(Uric ~ Sex, data = leye.cs)

# temporal
m9 <- lm(Uric ~ Julian, data = leye.cs)
m10 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m11 <- lm(Uric ~ Event, data = leye.cs)

# flock
m12 <- lm(Uric ~ Max_Flock_Size, data = leye.cs)

# null
m13 <- lm(Uric ~ 1, data = leye.cs)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results --> permanence top model but not significant (next is null)
summary(m6)
confint(m6)
