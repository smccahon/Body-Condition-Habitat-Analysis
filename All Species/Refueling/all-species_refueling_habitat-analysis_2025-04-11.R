#----------------------------------------------#
# All Species Refueling Rates Habitat Analysis #
#            Created 2025-04-11                #
#           Modified 2025-04-11                #
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

# filter birds that only contain metabolite information (n = 29)
birds <- birds %>% 
  filter(!is.na(Tri) & !is.na(Beta))


# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric), scale))

# ---------------------------------------------------------------------------- #

# Perform PCA ####

# Subset the dataset to only include 'Tri' and 'Beta'
birds_subset <- birds.cs[, c("Tri", "Beta")]

# Remove rows with NAs for PCA
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ]

# Run PCA on the cleaned data
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
birds.cs$PC1 <- NA  # Initialize with NA values
birds.cs$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds.cs[complete.cases(birds_subset), "PC1"] <- pca_scores[, 1]
birds.cs[complete.cases(birds_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(birds.cs)

# Merge PCA scores back to the original birds unstandardized dataset
# Create a data frame to store the PCA scores for the rows with no missing data
birds$PC1 <- NA  # Initialize with NA values
birds$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds[complete.cases(birds_subset), "PC1"] <- pca_scores[, 1]
birds[complete.cases(birds_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
# print(birds)

# View the first few rows of the updated data frame
# head(birds)

# ---------------------------------------------------------------------------- #

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
# % ag and dominant crop
# % veg and % exposed shoreline
# event and julian
# SPEI and julian

# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lmer(PC1 ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(PC1 ~ DominantCrop + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(PC1 ~ AgCategory + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)

# ag category is best (not significant), but % ag is still within 2 delta AICc

# temporal
m1 <- lmer(PC1 ~ Julian + MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(PC1 ~ Julian * MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m3 <- lmer(PC1 ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(PC1 ~ Julian + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model with interaction is more informative, but not much more so than 
# model without interaction. Drop interaction

# transformation for time needed? visually no ----------------------------------
plot(birds$seconds_since_midnight, birds$PC1)

# interactions between season and capture time? no
m1 <- lmer(PC1 ~ Event + seconds_since_midnight + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(PC1 ~ Event * seconds_since_midnight + (1|Species), data = birds.cs, 
           REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lmer(PC1 ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(PC1 ~ Percent_Total_Veg + (1|Species), data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(PC1 ~ Permanence + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(PC1 ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(PC1 ~ Dist_Closest_Wetland_m + (1|Species),  data = birds.cs, 
           REML = FALSE)

# weather
m6 <- lmer(PC1 ~ SPEI + (1|Species), 
           data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(PC1 ~ Sex + (1|Species), 
           data = birds.cs, REML = FALSE)
m8 <- lmer(PC1 ~ MigStatus + (1|Species), 
           data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(PC1 ~ Julian + (1|Species), 
           data = birds.cs, REML = FALSE)
m10 <- lmer(PC1 ~ Event + (1|Species), data = birds.cs,
            REML = FALSE)
m11 <- lmer(PC1 ~ seconds_since_midnight + (1 | Species), 
            data = birds.cs, REML = FALSE)

# flock
m12 <- lmer(PC1 ~ Max_Flock_Size + (1|Species), data = birds.cs, REML = FALSE)

# null
m13 <- lmer(PC1 ~ 1 + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
confint(m5) # distance to wetland top model but not informative. next is null

confint(m3)
confint(m11)

plot(birds$Time, birds$PC1)

# Variables to move on to chapter 2: dist. to wetland --------------------------
# top model is Dist. to closest wetland (but it's not informative)



# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- lmer(PC1 ~ Biomass + (1|Species), data = birds.sub, 
           REML = FALSE)

summary(m)
confint(m) # no effect of biomass

m <- lmer(PC1 ~ Diversity + (1|Species), data = birds.sub, 
          REML = FALSE)

summary(m)
confint(m) # no effect of diversity













