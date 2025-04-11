#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#               Created 2025-03-31                   #
#              Modified 2025-03-31                   #
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

# Transform Data ---------------------------------------------------------------

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

# Standardize continuous variables
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# Perform PCA ------------------------------------------------------------------

# Subset the dataset to only include 'Tri' and 'Beta'
leye_subset <- leye.cs[, c("Tri", "Beta")]

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
leye.cs$PC1 <- NA  # Initialize with NA values
leye.cs$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye.cs[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye.cs[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(leye.cs)

# Merge PCA scores back to the original LEYE unstandardized dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye$PC1 <- NA  # Initialize with NA values
leye$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
# print(leye)
# head(leye)

# Test for Correlations--------------------------------------------------------- 

# subset data
sample <- birds[, c("PercentAg",
                   "Percent_Total_Veg",
                   "Event",
                   "Sex",
                   "Biomass",
                   "Diversity",
                   "Permanence",
                   "MaxBufferWidth_m",
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
                   "PercentBufferAroundWetland",
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
# PercentAg & AgCategory (0.84)
# DominantCrop & PercentAg (0.72)
# AgCategory & DominantCrop (0.64)
# Percent_Total_Veg & Percent_Exposed_Shoreline (-0.61)
# PercentBuffer & MaxBufferWidth (0.68)
# Julian & SPEI (-0.67)

# ...drop correlated covariates ----
# agriculture
m1 <- lm(PC1 ~ PercentAg, data = leye.cs)
m2 <- lm(PC1 ~ AgCategory, data = leye.cs)
m3 <- lm(PC1 ~ DominantCrop, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)
# PercentAg is the most informative covariate

# vegetation
m1 <- lm(PC1 ~ Percent_Total_Veg, data = leye.cs)
m2 <- lm(PC1 ~ Percent_Exposed_Shoreline, data = leye.cs)
 
model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)
# If in the same model, exposed shoreline is a better predictor

# buffer
m1 <- lm(PC1 ~ PercentBufferAroundWetland, data = leye.cs)
m2 <- lm(PC1 ~ MaxBufferWidth_m, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)
# max buffer width is a better predictor

# Julian & SPEI
m1 <- lm(PC1 ~ Julian, data = leye.cs)
m2 <- lm(PC1 ~ SPEI, data = leye.cs)
m3 <- lm(PC1 ~ DaysIntoSeason, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)
# If in the same model, Julian is a better predictor

# Model Selection (Without Site)------------------------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg + I(PercentAg^2), data = leye.cs)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg, data = leye.cs)
m3 <- lm(PC1 ~ MaxBufferWidth_m, data = leye.cs)

# habitat
m4 <- lm(PC1 ~ Permanence, data = leye.cs)
m5 <- lm(PC1 ~ Percent_Exposed_Shoreline, data = leye.cs)
m6 <- lm(PC1 ~ Dist_Closest_Wetland_m, data = leye.cs)

# weather
m7 <- lm(PC1 ~ SPEI, data = leye.cs)

# life history
m8 <- lm(PC1 ~ Age, data = leye.cs)
m9 <- lm(PC1 ~ Sex, data = leye.cs)

# temporal
m10 <- lm(PC1 ~ Julian, data = leye.cs)
m11 <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2) + I(ts.sunrise^3), data = leye.cs)
m12 <- lm(PC1 ~ Event, data = leye.cs)

# flock effect
m13 <- lm(PC1 ~ Max_Flock_Size, data = leye.cs)

# null
m14 <- lm(PC1 ~ 1, data = leye.cs)

model_names <- paste0("m", 1:14)

models <- mget(model_names)

aictab(models, modnames = model_names)

# significant predictors:
# permanence, time, % ag

# Model Selection (Top Covariates)----------------------------------------------
m1 <- lm(PC1 ~ PercentAg, data = leye.cs)
m2 <- lm(PC1 ~ Permanence, data = leye.cs)
m3 <- lm(PC1 ~ ts.sunrise, data = leye.cs)

m4 <- lm(PC1 ~ PercentAg + Permanence, data = leye.cs)
m5 <- lm(PC1 ~ PercentAg + ts.sunrise, data = leye.cs)
m6 <- lm(PC1 ~ ts.sunrise + Permanence, data = leye.cs)

m7 <- lm(PC1 ~ PercentAg + Permanence + ts.sunrise, data = leye.cs)

m8 <- lm(PC1 ~ 1, data= leye.cs)

model_names <- paste0("m", 1:8)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m5)
summary(m6)

confint(m5) 

# should time be informed null?? include time in every model

# Plot models ------------------------------------------------------------------
# ...time & permanence ----
m <- lm(PC1 ~ ts.sunrise + Permanence, data = leye)

d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise), 
                                              max(leye$ts.sunrise), 
                                              length.out = 1000),  
                 Permanence = unique(leye$Permanence))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$yhat <- predictions$fit

d$lower_CI <- d$yhat - 1.96 * predictions$se.fit
d$upper_CI <- d$yhat + 1.96 * predictions$se.fit

ggplot(d, aes(x = ts.sunrise, y = yhat, col = Permanence)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI,
                  fill = Permanence), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "Time of Capture Since Sunrise (min)", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2, 
           label = "Sunrise", angle = 90, 
           vjust = -0.5, hjust = 0.5,
           size = 5) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  scale_color_viridis_d(alpha = 1, begin = 0, end = 1) +
  scale_fill_viridis_d(alpha = 1, begin = 0, end = 1) +
  geom_point(data = leye, aes(x = ts.sunrise, y = PC1), size = 3)

# RELATIONSHIP WITH TIME IS NONLINEAR

# ...time and % agriculture ----
m <- lm(PC1 ~ ts.sunrise + PercentAg, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                  max(leye$PercentAg), 
                                  length.out = 1000),  
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d, se.fit = TRUE)

d$yhat <- predictions$fit

d$lower_CI <- d$yhat - 1.96 * predictions$se.fit
d$upper_CI <- d$yhat + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, color = NA) + 
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), size = 2) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1)

# % AGRICULTURE MAY ALSO BE NONLINEAR

# Model Diagnostics ------------------------------------------------------------
m5 <- lm(PC1 ~ PercentAg + ts.sunrise, data = leye.cs)
m6 <- lm(PC1 ~ ts.sunrise + Permanence, data = leye.cs)

par(mfrow = c(2,2))
plot(m5)
plot(m6)

par(mfrow = c(1,1))
plot(predict(m5), rstudent(m5))

