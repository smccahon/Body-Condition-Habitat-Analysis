#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#             Non-linear regression                  #
#               Created 2025-04-03                   #
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

# filter birds that only contain metabolite information (n = 29)
leye <- leye %>% 
  filter(!is.na(Tri) & !is.na(Beta))

# is there enough variation in buffer presence?
# table(leye$Buffered) # nope
# table(birds$Buffered) # also nope

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

# should i be standardizing my data?
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

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

# ...covariates with > 0.6 correlation ----
# All Birds:
# PercentAg & AgCategory (0.84)
# DominantCrop & PercentAg (0.72)
# AgCategory & DominantCrop (0.64)
# Percent_Total_Veg & Percent_Exposed_Shoreline (-0.61)
# Julian & SPEI (-0.67)

# LEYE refueling specifically (n = 29):
# PercentVeg + Dominant Crop (-0.71)
# PercentAg & AgCategory (0.80)
# PercentAg + Dominant Crop (0.62)
# NearestCropDistance + PercentAg (-0.656)
# permanence and exposed shoreline (-0.62)
# Ag Category and permanence (-0.94)
# SPEI & %exposed shoreline (-0.63)
# Ag category & nearest crop distance (-0.73)
# SPEI and permanence (-0.94)
# Max flock size and dominant crop (0.69)


# which correlated covariates should I drop? REDO THIS WHOLE SECTION FOR LEYE---

# agriculture
m1 <- lm(PC1 ~ PercentAg, data = leye)
m2 <- lm(PC1 ~ DominantCrop, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a much better predictor


# which data variable to use for LEYE? -----------------------------------------
m1 <- lm(PC1 ~ DaysIntoSeason, data = leye)
m2 <- lm(PC1 ~ Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Julian is better

# Modeling time & Percent Ag (linear)----------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + seconds_since_midnight +
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) +
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + seconds_since_midnight +
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) +
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# habitat
m3 <- lm(PC1 ~ Permanence + seconds_since_midnight + 
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) +
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + seconds_since_midnight + 
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + seconds_since_midnight +  
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) + 
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# weather
m6 <- lm(PC1 ~ SPEI + seconds_since_midnight + 
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) + 
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# life history
m7 <- lm(PC1 ~ Age + seconds_since_midnight + 
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) + 
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)
m8 <- lm(PC1 ~ Sex + seconds_since_midnight + 
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) + 
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# temporal
m9 <- lm(PC1 ~ Julian + seconds_since_midnight +            
           sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
           cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)
m10 <- lm(PC1 ~ Event + seconds_since_midnight +            
            sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
            cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)
m11 <- lm(PC1 ~ seconds_since_midnight +            
            sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
            cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + seconds_since_midnight +           
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
          cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)


# STANDARDIZED ##
# agriculture
m1 <- lm(PC1 ~ PercentAg + seconds_since_midnight + sin + cos, data = leye.cs)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + seconds_since_midnight +
           sin +
           cos, data = leye.cs)

# habitat
m3 <- lm(PC1 ~ Permanence + seconds_since_midnight + 
           sin +
           cos, data = leye.cs)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + seconds_since_midnight + 
           sin +             
           cos, data = leye.cs)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + seconds_since_midnight +  
           sin + 
           cos, data = leye.cs)

# weather
m6 <- lm(PC1 ~ SPEI + seconds_since_midnight + 
           sin + 
           cos, data = leye.cs)

# life history
m7 <- lm(PC1 ~ Age + seconds_since_midnight + 
           sin + 
           cos, data = leye.cs)
m8 <- lm(PC1 ~ Sex + seconds_since_midnight + 
           sin + 
           cos, data = leye.cs)

# temporal
m9 <- lm(PC1 ~ Julian + seconds_since_midnight +            
           sin +             
           cos, data = leye.cs)
m10 <- lm(PC1 ~ Event + seconds_since_midnight +            
            sin +             
            cos, data = leye.cs)
m11 <- lm(PC1 ~ seconds_since_midnight +            
            sin +             
            cos, data = leye.cs)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + seconds_since_midnight +           
            sin +             
            cos, data = leye.cs)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# important parameters: time only

confint(m8)
plot(predict(m11),rstudent(m11))


# Modeling time & Percent Ag (quadratic)-------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + I(PercentAg^2) + seconds_since_midnight +
           sin + cos, data = leye.cs)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + seconds_since_midnight +
           sin + cos, data = leye.cs)

# habitat
m3 <- lm(PC1 ~ Permanence + seconds_since_midnight +
           sin + cos, data = leye.cs)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + seconds_since_midnight +
           sin + cos, data = leye.cs)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + seconds_since_midnight +
           sin + cos, data = leye.cs)

# weather
m6 <- lm(PC1 ~ SPEI + seconds_since_midnight +
           sin + cos, data = leye.cs)

# life history
m7 <- lm(PC1 ~ Age + seconds_since_midnight + sin + cos, data = leye.cs)
m8 <- lm(PC1 ~ Sex + seconds_since_midnight + sin + cos, data = leye.cs)

# temporal
m9 <- lm(PC1 ~ Julian + seconds_since_midnight + sin + cos, data = leye.cs)
m10 <- lm(PC1 ~ Event + seconds_since_midnight + sin + cos, data = leye.cs)
m11 <- lm(PC1 ~ seconds_since_midnight + sin + cos, data = leye.cs)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + seconds_since_midnight +
            sin + cos, data = leye.cs)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# important parameters: surrounding ag (quadratic), time, permanence,
#                       dist to closest wetland, shorebird abundance


confint(m5)


# Modeling time & Percent Ag (exponential)--------------------------------------

## CAN'T FIND GOOD STARTING VALUES!!
lm_model <- lm(PC1 ~ seconds_since_midnight + sin + cos, data = leye.cs)

m1 <- nls(PC1 ~ beta0 + beta1 * -exp(PercentAg/a), 
          data = leye.cs,
          start = list(beta0 = 2, beta1 = 2, a = 40),
          control = nls.control(maxiter = 500))

m1 <- nls(PC1 ~ beta0 + beta1 * -exp(PercentAg/a) + 
            beta2 * seconds_since_midnight + 
            beta3 * sin +
            beta4 * cos, 
          data = leye.cs,
          start = list(beta0 = 2, beta1 = 2, a = 40, 
                       beta2 = 0, 
                       beta3 = 0,
                       beta4 = 0),
          control = nls.control(maxiter = 1000))

m1 <- nls(PC1 ~ beta0 + beta1 * -exp(PercentAg/a) +
            beta2 * seconds_since_midnight +
  beta3 * sin(2 * pi * seconds_since_midnight / (24 * 3600)) +
  beta4 * cos(2 * pi * seconds_since_midnight / (24 * 3600)),
  data = leye,
  start = list(beta0 = 2, beta1 = 2, a = 40, beta2 = 0, beta3 = 0,
              beta4 = 0),
  control = nls.control(maxiter = 100))

summary(m1)
