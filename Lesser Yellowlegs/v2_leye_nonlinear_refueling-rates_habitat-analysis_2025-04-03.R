#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#             Non-linear regression                  #
#               Created 2025-04-01                   #
#              Modified 2025-04-03                   #
#----------------------------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(trtools)
library(lme4)
library(car)
library(viridis)


# read data
birds <- read.csv("Body_Condition_Habitat_Analysis_2025-03-31.csv")

# ...make new columns ----
# neonicotinoid detection column
birds$Detection <- ifelse(birds$OverallNeonic > 0, 
                          "Detection", "Non-detection")

# convert capture time to number of minutes after sunrise
birds <- birds %>%
  mutate(
    hour_of_day = as.numeric(format(strptime(Time, "%H:%M"), "%H")),  # Extract hour (0-23)
    minute_of_day = as.numeric(format(strptime(Time, "%H:%M"), "%M")),  # Extract minute (0-59)
    time_in_minutes = hour_of_day * 60 + minute_of_day  # Total time in minutes
  ) %>% 
  mutate(
    hour_of_day_s = as.numeric(format(strptime(Sunrise, "%H:%M"), "%H")),  # Extract hour (0-23)
    minute_of_day_s = as.numeric(format(strptime(Sunrise, "%H:%M"), "%M")),  # Extract minute (0-59)
    time_in_minutes_s = hour_of_day_s * 60 + minute_of_day_s  # Total time in minutes
  ) %>% 
  mutate(
    ts.sunrise = time_in_minutes - time_in_minutes_s)

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

leye$Event <- factor(leye$Event, 
                     levels = c("Fall 2021", "Spring 2022", "Fall 2023"),
                     labels = c("Fall 2021", "Spring 2022", "Fall 2023"))

# filter birds that only contain metabolite information (n = 29)
leye <- leye %>% 
  filter(!is.na(Tri) & !is.na(Beta))

# is there enough variation in buffer presence?
table(leye$Buffered) # nope
table(birds$Buffered) # also nope

# Perform PCA ------------------------------------------------------------------

# Subset the dataset to only include 'Tri' and 'Beta'
leye_subset <- leye[, c("Tri", "Beta")]

# Remove rows with NAs for PCA
leye_subset_clean <- leye_subset[complete.cases(leye_subset), ]

# Run PCA on the cleaned data
pca_result <- prcomp(leye_subset_clean, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye$PC1 <- NA  # Initialize with NA values
leye$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(leye)

# Merge PCA scores back to the original LEYE unstandardized dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye$PC1 <- NA  # Initialize with NA values
leye$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

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
                    "ts.sunrise",
                    "Site",
                    "AgCategory",
                    "SPEI",
                    "Julian",
                    "Age",
                    "DominantCrop",
                    "NearestCropDistance_m",
                    "Dist_Closest_Wetland_m",
                    "Max_Flock_Size",
                    "ts.sunrise"
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

# ...covariates with > 0.6 correlation ----
# PercentAg & AgCategory (0.84)
# DominantCrop & PercentAg (0.72)
# AgCategory & DominantCrop (0.64)
# Percent_Total_Veg & Percent_Exposed_Shoreline (-0.61)
# Julian & SPEI (-0.67)

# which correlated covariates should I drop? -----------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg, data = leye)
m2 <- lm(PC1 ~ DominantCrop, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a much better predictor

# weather and time
m1 <- lm(PC1 ~ SPEI, data = leye)
m2 <- lm(PC1 ~ Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# if included in the same model, Julian better predictor

# habitat and vegetation
m1 <- lm(PC1 ~ Percent_Total_Veg, data = leye)
m2 <- lm(PC1 ~ Percent_Exposed_Shoreline, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# if included in the same model, % exposed shoreline better predictor

# which data variable to use for LEYE? -----------------------------------------
m1 <- lm(PC1 ~ DaysIntoSeason, data = leye)
m2 <- lm(PC1 ~ Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Julian is better

# Modeling time (linear) & Percent Ag (linear)----------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + ts.sunrise, data = leye)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + ts.sunrise, data = leye)

# habitat
m3 <- lm(PC1 ~ Permanence + ts.sunrise, data = leye)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + ts.sunrise, data = leye)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + ts.sunrise, data = leye)

# weather
m6 <- lm(PC1 ~ SPEI + ts.sunrise, data = leye)

# life history
m7 <- lm(PC1 ~ Age + ts.sunrise, data = leye)
m8 <- lm(PC1 ~ Sex + ts.sunrise, data = leye)

# temporal
m9 <- lm(PC1 ~ Julian + ts.sunrise, data = leye)
m10 <- lm(PC1 ~ Event + ts.sunrise, data = leye)
m11 <- lm(PC1 ~ ts.sunrise, data = leye)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + ts.sunrise, data = leye)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# important parameters: permanence, dist to wetland, surrounding ag,
#                       shorebird abundance, time

confint(m4)


# Modeling time (linear) & Percent Ag (quadratic)-------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + I(PercentAg^2) + ts.sunrise, data = leye)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + ts.sunrise, data = leye)

# habitat
m3 <- lm(PC1 ~ Permanence + ts.sunrise, data = leye)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + ts.sunrise, data = leye)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + ts.sunrise, data = leye)

# weather
m6 <- lm(PC1 ~ SPEI + ts.sunrise, data = leye)

# life history
m7 <- lm(PC1 ~ Age + ts.sunrise, data = leye)
m8 <- lm(PC1 ~ Sex + ts.sunrise, data = leye)

# temporal
m9 <- lm(PC1 ~ Julian + ts.sunrise, data = leye)
m10 <- lm(PC1 ~ Event + ts.sunrise, data = leye)
m11 <- lm(PC1 ~ ts.sunrise, data = leye)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + ts.sunrise, data = leye)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# important parameters: surrounding ag (quadratic), time, permanence,
#                       dist to closest wetland, shorebird abundance


confint(m11)


# Modeling (time^3) & Percent Ag (linear)---------------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# habitat
m3 <- lm(PC1 ~ Permanence + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# weather
m6 <- lm(PC1 ~ SPEI + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# life history
m7 <- lm(PC1 ~ Age + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m8 <- lm(PC1 ~ Sex + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# temporal
m9 <- lm(PC1 ~ Julian + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m10 <- lm(PC1 ~ Event + ts.sunrise + I(ts.sunrise^2) +
            I(ts.sunrise^3), data = leye)
m11 <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2) +
            I(ts.sunrise^3), data = leye)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + ts.sunrise + I(ts.sunrise^2) +
            I(ts.sunrise^3), data = leye)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m1)
confint(m1)

# important parameters: time, (dist to wetland marginally significant)

# Modeling (time^3) & Percent Ag (quadratic)------------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + I(PercentAg^2) + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# habitat
m3 <- lm(PC1 ~ Permanence + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# weather
m6 <- lm(PC1 ~ SPEI + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# life history
m7 <- lm(PC1 ~ Age + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m8 <- lm(PC1 ~ Sex + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)

# temporal
m9 <- lm(PC1 ~ Julian + ts.sunrise + I(ts.sunrise^2) +
           I(ts.sunrise^3), data = leye)
m10 <- lm(PC1 ~ Event + ts.sunrise + I(ts.sunrise^2) +
            I(ts.sunrise^3), data = leye)
m11 <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2) +
            I(ts.sunrise^3), data = leye)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + ts.sunrise + I(ts.sunrise^2) +
            I(ts.sunrise^3), data = leye)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# important parameters: 


confint(m10)

