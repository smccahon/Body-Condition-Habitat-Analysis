#----------------------------------------------#
# All Species Refueling Rates Habitat Analysis #
#            Created 2025-05-05                #
#           Modified 2025-05-05                #
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
library(nlme)

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


# ---------------------------------------------------------------------------- #

# Perform PCA ####

# Subset the dataset to only include 'Tri' and 'Beta'
birds_subset <- birds[, c("Tri", "Beta")]

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
birds$PC1 <- NA  # Initialize with NA values
birds$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds[complete.cases(birds_subset), "PC1"] <- pca_scores[, 1]
birds[complete.cases(birds_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(birds)

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

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("PC1"), scale))


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

# cyclical or linear?
m1 <- lmer(PC1 ~ seconds_since_midnight + 
             sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
             cos(2 * pi * seconds_since_midnight / (24 * 3600)) + (1|Species),
           REML = FALSE, data = birds)

m2 <- lmer(PC1 ~ seconds_since_midnight + (1|Species), REML = FALSE, data = birds)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names) # linear is much better 
                                       # but i dont know if thats appropriate

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


# Is random effect of species necessary? I guess not, but inferences don't change ----

# agriculture
m1 <- lm(PC1 ~ PercentAg, data = birds.cs)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg, data = birds.cs)

# habitat
m3 <- lm(PC1 ~ Permanence, data = birds.cs)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline, 
           data = birds.cs)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m,  data = birds.cs)

# weather
m6 <- lm(PC1 ~ SPEI, 
           data = birds.cs)

# life history
m7 <- lm(PC1 ~ Sex, 
           data = birds.cs)
m8 <- lm(PC1 ~ MigStatus, 
           data = birds.cs)

# temporal
m9 <- lm(PC1 ~ Julian, data = birds.cs)
m10 <- lm(PC1 ~ Event, data = birds.cs)
m11 <- lm(PC1 ~ seconds_since_midnight, data = birds.cs)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size, data = birds.cs)

# null
m13 <- lm(PC1 ~ 1, data = birds.cs)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# plot relationships ----
m <- lmer(PC1 ~ PercentAg + (1 | Species), REML = FALSE, data = birds)

d <- expand.grid(PercentAg = seq(min(birds$PercentAg), 
                                              max(birds$PercentAg), 
                                              length.out = 1000),
                 Species = unique(birds$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "All Species Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = birds, aes(x = PercentAg, y = PC1, color = Species)) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1)


m <- lmer(PC1 ~ seconds_since_midnight + (1 | Species), REML = FALSE, data = birds)

d <- expand.grid(seconds_since_midnight = seq(min(birds$seconds_since_midnight), 
                                 max(birds$seconds_since_midnight), 
                                 length.out = 1000),
                 Species = unique(birds$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = seconds_since_midnight, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Time", 
       y = "All Species Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = birds, aes(x = seconds_since_midnight, y = PC1, 
                               color = Species)) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  scale_x_time(labels = scales::time_format("%H:%M"),
               breaks = seq(0, 86400, by = 7200)) 
  # scale_color_viridis_d(alpha = 1, begin = 1, end = 0)

# ...time as cyclical ----

m <- lmer(PC1 ~ seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)) + (1|Species),
          REML = FALSE, data = birds)

d <- expand.grid(seconds_since_midnight = seq(min(birds$seconds_since_midnight), 
                                              max(birds$seconds_since_midnight), 
                                              length.out = 1000),
                 Species = unique(birds$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

# time
ggplot(d, aes(x = seconds_since_midnight, y = fit)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "Capture Time", 
       y = "All Species Fattening Index",
       col = "% Surrounding Agriculture") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = birds, aes(x = seconds_since_midnight, y = PC1, 
                                   col = PercentAg), size = 3) +
  scale_x_time(labels = scales::time_format("%H:%M"),
               breaks = seq(0, 86400, by = 7200)) +
  scale_color_viridis(alpha = 1, begin = 1, end = 0) +
  geom_text(data = birds, aes(x = seconds_since_midnight, y = PC1, 
                                  label = sprintf("%.0f", PercentAg)), 
            size = 3, vjust = -0.5) 



# Exclude early morning birds and Rerun Analysis ----
birds_late <- subset(birds, formatted_time > "08:00")
birds_late.cs <- birds_late %>%
  mutate(across(where(is.numeric), scale))


# agriculture
m1 <- lmer(PC1 ~ PercentAg + (1|Species), data = birds_late.cs, REML = FALSE)

# vegetation
m2 <- lmer(PC1 ~ Percent_Total_Veg + (1|Species), data = birds_late.cs, REML = FALSE)

# habitat
m3 <- lmer(PC1 ~ Permanence + (1|Species), data = birds_late.cs, REML = FALSE)
m4 <- lmer(PC1 ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds_late.cs, REML = FALSE)
m5 <- lmer(PC1 ~ Dist_Closest_Wetland_m + (1|Species),  data = birds_late.cs, 
           REML = FALSE)

# weather
m6 <- lmer(PC1 ~ SPEI + (1|Species), 
           data = birds_late.cs, REML = FALSE)

# life history
m7 <- lmer(PC1 ~ Sex + (1|Species), 
           data = birds_late.cs, REML = FALSE)
m8 <- lmer(PC1 ~ MigStatus + (1|Species), 
           data = birds_late.cs, REML = FALSE)

# temporal
m9 <- lmer(PC1 ~ Julian + (1|Species), 
           data = birds_late.cs, REML = FALSE)
m10 <- lmer(PC1 ~ Event + (1|Species), data = birds_late.cs,
            REML = FALSE)
m11 <- lmer(PC1 ~ seconds_since_midnight + (1 | Species), 
            data = birds_late.cs, REML = FALSE)

# flock
m12 <- lmer(PC1 ~ Max_Flock_Size + (1|Species), data = birds_late.cs, REML = FALSE)

# null
m13 <- lmer(PC1 ~ 1 + (1 | Species), data = birds_late.cs, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names) # null is top model

confint(m4)
confint(m5)

m <- lmer(PC1 ~ seconds_since_midnight + (1|Species),
          REML = FALSE, data = birds_late)

d <- expand.grid(seconds_since_midnight = seq(min(birds_late$seconds_since_midnight), 
                                              max(birds_late$seconds_since_midnight), 
                                              length.out = 1000),
                 Species = unique(birds_late$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

# time
ggplot(d, aes(x = seconds_since_midnight, y = fit)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "Capture Time", 
       y = "All Species Fattening Index",
       col = "% Surrounding Agriculture") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = birds_late, aes(x = seconds_since_midnight, y = PC1, 
                               col = PercentAg), size = 3) +
  scale_x_time(labels = scales::time_format("%H:%M"),
               breaks = seq(0, 86400, by = 7200)) +
  scale_color_viridis(alpha = 1, begin = 1, end = 0) +
  geom_text(data = birds_late, aes(x = seconds_since_midnight, y = PC1, 
                              label = sprintf("%.0f", PercentAg)), 
            size = 3, vjust = -0.5) 

# agriculture
m <- lmer(PC1 ~ PercentAg + (1|Species),
          REML = FALSE, data = birds_late)

d <- expand.grid(PercentAg = seq(min(birds_late$PercentAg), 
                                              max(birds_late$PercentAg), 
                                              length.out = 1000),
                 Species = unique(birds_late$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "All Species Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = birds_late, aes(x = PercentAg, y = PC1), size = 3)

# SPEI
m <- lmer(PC1 ~ SPEI + (1|Species),
          REML = FALSE, data = birds_late)

d <- expand.grid(SPEI = seq(min(birds_late$SPEI), 
                                 max(birds_late$SPEI), 
                                 length.out = 1000),
                 Species = unique(birds_late$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "Drought Condition Index", 
       y = "All Species Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = birds_late, aes(x = SPEI, y = PC1), size = 3)




# PLOT RESIDUALS AGAINST AG ----
m <- lmer(PC1 ~ seconds_since_midnight + (1|Species), REML = FALSE, 
          data = birds.cs)

birds.cs$residuals <- resid(m)

ggplot(birds.cs, aes(x = PercentAg, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 0.75) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_classic() +
  labs(x = "% Surrounding Agriculture",
       y = "Residuals of Fattening Index ~ Time",
       title = "Unexplained Variation in Fattening Index vs. % Agriculture") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

summary(lm(residuals ~ PercentAg, data = birds.cs)) # Adj. R2 = -0.0102

# agriculture
m1 <- lmer(PC1 ~ PercentAg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(PC1 ~ Percent_Total_Veg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(PC1 ~ Permanence  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)
m4 <- lmer(PC1 ~ Percent_Exposed_Shoreline  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m5 <- lmer(PC1 ~ Dist_Closest_Wetland_m  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(PC1 ~ SPEI  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(PC1 ~ Sex  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m8 <- lmer(PC1 ~ MigStatus  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(PC1 ~ Julian  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m10 <- lmer(PC1 ~ Event  + (1|Species) + (1|Site), data = birds.m,
            REML = FALSE)
m11 <- lmer(PC1 ~ seconds_since_midnight + (1 | Species) + (1|Site), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(PC1 ~ Max_Flock_Size  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# null
m13 <- lmer(PC1 ~ 1 + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names) # results do not change


# run reduced dataset with site as a random effect ----
# agriculture
m1 <- lmer(PC1 ~ PercentAg + (1|Species), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(PC1 ~ Percent_Total_Veg + (1|Species), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(PC1 ~ Permanence + (1|Species), data = birds.m, REML = FALSE)
m4 <- lmer(PC1 ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds.m, REML = FALSE)
m5 <- lmer(PC1 ~ Dist_Closest_Wetland_m + (1|Species),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(PC1 ~ SPEI + (1|Species), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(PC1 ~ Sex + (1|Species), 
           data = birds.m, REML = FALSE)
m8 <- lmer(PC1 ~ MigStatus + (1|Species), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(PC1 ~ Julian + (1|Species), 
           data = birds.m, REML = FALSE)
m10 <- lmer(PC1 ~ Event + (1|Species), data = birds.m,
            REML = FALSE)
m11 <- lmer(PC1 ~ seconds_since_midnight + (1 | Species), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(PC1 ~ Max_Flock_Size + (1|Species), data = birds.m, REML = FALSE)

# null
m13 <- lmer(PC1 ~ 1 + (1 | Species), data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m1 <- lmer(PC1 ~ Dist_Closest_Wetland_m  + (1|Species),  data = birds.m, 
           REML = FALSE)
m2 <- lmer(PC1 ~ Dist_Closest_Wetland_m  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

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




# any relationships with beta alone?----

# filter birds that only contain metabolite information (n = 89)
birds <- birds %>% 
  filter(!is.na(Beta))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# relabel species names
birds$Species <- factor(birds$Species, 
                        levels = c("Killdeer", "Least Sandpiper", "Lesser Yellowlegs",
                                   "Longbilled Dowitcher", "Pectoral Sandpiper",
                                   "Semipalmated Sandpiper", "Willet",
                                   "Wilsons Phalarope"),
                        labels = c("Killdeer", "Least Sandpiper", "Lesser Yellowlegs",
                                   "Long-billed Dowitcher", "Pectoral Sandpiper",
                                   "Semipalmated Sandpiper", "Willet",
                                   "Wilson's Phalarope"))

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("Beta"), scale))


# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lmer(Beta ~ PercentAg + seconds_since_midnight + (1|Species), data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(Beta ~ Percent_Total_Veg + seconds_since_midnight + (1|Species), data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(Beta ~ Permanence + seconds_since_midnight + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Beta ~ Percent_Exposed_Shoreline + seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + (1|Species),  data = birds.cs, 
           REML = FALSE)

# weather
m6 <- lmer(Beta ~ SPEI + seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(Beta ~ Sex + seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)
m8 <- lmer(Beta ~ MigStatus + seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(Beta ~ Julian + seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)
m10 <- lmer(Beta ~ Event + seconds_since_midnight + (1|Species), data = birds.cs,
            REML = FALSE)
m11 <- lmer(Beta ~ seconds_since_midnight + (1 | Species), 
            data = birds.cs, REML = FALSE)

# flock
m12 <- lmer(Beta ~ Max_Flock_Size + seconds_since_midnight + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names) # time + % ag top model and only within 2 AICc

confint(m1)

plot(birds$Beta, birds)

summary(m1)
confint(m1)

m <- lmer(Beta ~ PercentAg + (1|Species) + seconds_since_midnight, data = birds,
          REML = FALSE)

d <- expand.grid(PercentAg = seq(min(birds$PercentAg), 
                                              max(birds$PercentAg), 
                                              length.out = 1000),
                 seconds_since_midnight = mean(birds$seconds_since_midnight),
                 Species = unique(birds$Species))

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit


ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "All Species Beta-Hydroxybutyrate Levels") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = birds, aes(x = PercentAg, y = Beta, color = Species), size = 3)

plot(m) # severe heteroscedasticity

# weighted least squares
birds.cs %>% summarize(variance = var(Beta), weight = 1/var(Beta))

# put results back into data frame
birds.cs <- birds.cs %>% 
  group_by(Species) %>% 
  mutate(w = 1/var(Beta)) %>% 
  ungroup()

m.wls <- lmer(Beta ~ PercentAg + (1|Species) + seconds_since_midnight,
              data = birds.cs,
              REML = FALSE,
              weights = w)

# fixed; much better
plot(predict(m.wls),rstudent(m.wls)) # MUCH better

summary(m.wls)
confint(m.wls)

confint(m.wls, method = "Wald")


# weighted least squares
birds %>% summarize(variance = var(Beta), weight = 1/var(Beta))

# put results back into data frame
birds <- birds %>% 
  group_by(Species) %>% 
  mutate(w = 1/var(Beta)) %>% 
  ungroup()



m.wls <- lmer(Beta ~ PercentAg + (1|Species) + seconds_since_midnight,
              data = birds,
              REML = FALSE,
              weights = w)


# Wald (asymptotic confidence intervals), based on fixed effect estimates
confint(m.wls, method = "Wald")

d <- expand.grid(PercentAg = seq(min(birds$PercentAg), 
                                 max(birds$PercentAg), 
                                 length.out = 1000),
                 seconds_since_midnight = mean(birds$seconds_since_midnight),
                 Species = unique(birds$Species))

predictions <- predict(m.wls, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit


ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "BHB Levels (mmol/L)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = birds, aes(x = PercentAg, y = Beta, color = Species), 
             size = 3) +
  scale_color_viridis_d(begin = 0, end = 1, alpha = 0.8)
  

plot(m.wls)

# is random effect of site necessary? ----
# exclude sites with only 1 bird
birds.m <- birds.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()


# agriculture
m1 <- lmer(Beta ~ PercentAg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(Beta ~ Percent_Total_Veg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(Beta ~ Permanence  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)
m4 <- lmer(Beta ~ Percent_Exposed_Shoreline  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m5 <- lmer(Beta ~ Dist_Closest_Wetland_m  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(Beta ~ SPEI  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(Beta ~ Sex  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m8 <- lmer(Beta ~ MigStatus  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(Beta ~ Julian  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m10 <- lmer(Beta ~ Event  + (1|Species) + (1|Site), data = birds.m,
            REML = FALSE)
m11 <- lmer(Beta ~ seconds_since_midnight + (1 | Species) + (1|Site), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(Beta ~ Max_Flock_Size  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# null
m13 <- lmer(Beta ~ 1 + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names) # results do not change, some models have singular fit warning


# run reduced dataset with site as a random effect ----
# agriculture
m1 <- lmer(Beta ~ PercentAg + (1|Species), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(Beta ~ Percent_Total_Veg + (1|Species), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(Beta ~ Permanence + (1|Species), data = birds.m, REML = FALSE)
m4 <- lmer(Beta ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds.m, REML = FALSE)
m5 <- lmer(Beta ~ Dist_Closest_Wetland_m + (1|Species),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(Beta ~ SPEI + (1|Species), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(Beta ~ Sex + (1|Species), 
           data = birds.m, REML = FALSE)
m8 <- lmer(Beta ~ MigStatus + (1|Species), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(Beta ~ Julian + (1|Species), 
           data = birds.m, REML = FALSE)
m10 <- lmer(Beta ~ Event + (1|Species), data = birds.m,
            REML = FALSE)
m11 <- lmer(Beta ~ seconds_since_midnight + (1 | Species), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(Beta ~ Max_Flock_Size + (1|Species), data = birds.m, REML = FALSE)

# null
m13 <- lmer(Beta ~ 1 + (1 | Species), data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m1 <- lmer(Beta ~ PercentAg  + (1|Species),  data = birds.m, 
           REML = FALSE)
m2 <- lmer(Beta ~ PercentAg  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

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



# any relationships with Tri alone? ----

# filter birds that only contain metabolite information (n = 89)
birds <- birds %>% 
  filter(!is.na(Tri))

# relabel species names
birds$Species <- factor(birds$Species, 
                        levels = c("Killdeer", "Least Sandpiper", "Lesser Yellowlegs",
                                   "Longbilled Dowitcher", "Pectoral Sandpiper",
                                   "Semipalmated Sandpiper", "Willet",
                                   "Wilsons Phalarope"),
                        labels = c("Killdeer", "Least Sandpiper", "Lesser Yellowlegs",
                                   "Long-billed Dowitcher", "Pectoral Sandpiper",
                                   "Semipalmated Sandpiper", "Willet",
                                   "Wilson's Phalarope"))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("Tri"), scale))


# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lmer(Tri ~ PercentAg +  (1|Species), data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(Tri ~ Percent_Total_Veg +  (1|Species), data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(Tri ~ Permanence +  (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Tri ~ Percent_Exposed_Shoreline +  (1|Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(Tri ~ Dist_Closest_Wetland_m +  (1|Species),  data = birds.cs, 
           REML = FALSE)

# weather
m6 <- lmer(Tri ~ SPEI +  (1|Species), 
           data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(Tri ~ Sex +  (1|Species), 
           data = birds.cs, REML = FALSE)
m8 <- lmer(Tri ~ MigStatus +  (1|Species), 
           data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(Tri ~ Julian +  (1|Species), 
           data = birds.cs, REML = FALSE)
m10 <- lmer(Tri ~ Event +  (1|Species), data = birds.cs,
            REML = FALSE)
m11 <- lmer(Tri ~ seconds_since_midnight + (1 | Species), 
            data = birds.cs, REML = FALSE)

# flock
m12 <- lmer(Tri ~ Max_Flock_Size + (1|Species), data = birds.cs, REML = FALSE)

# null
m13 <- lmer(Tri ~ 1 + (1|Species),data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names) # no significant parameters

confint(m10)


# is random effect of site necessary? ----
# exclude sites with only 1 bird
birds.m <- birds.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()


# agriculture
m1 <- lmer(Tri ~ PercentAg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(Tri ~ Percent_Total_Veg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(Tri ~ Permanence  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)
m4 <- lmer(Tri ~ Percent_Exposed_Shoreline  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m5 <- lmer(Tri ~ Dist_Closest_Wetland_m  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(Tri ~ SPEI  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(Tri ~ Sex  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m8 <- lmer(Tri ~ MigStatus  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(Tri ~ Julian  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m10 <- lmer(Tri ~ Event  + (1|Species) + (1|Site), data = birds.m,
            REML = FALSE)
m11 <- lmer(Tri ~ seconds_since_midnight + (1 | Species) + (1|Site), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(Tri ~ Max_Flock_Size  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# null
m13 <- lmer(Tri ~ 1 + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names) # results do not change, most models have singular fit warning


# run reduced dataset with site as a random effect ----
# agriculture
m1 <- lmer(Tri ~ PercentAg + (1|Species), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(Tri ~ Percent_Total_Veg + (1|Species), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(Tri ~ Permanence + (1|Species), data = birds.m, REML = FALSE)
m4 <- lmer(Tri ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds.m, REML = FALSE)
m5 <- lmer(Tri ~ Dist_Closest_Wetland_m + (1|Species),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(Tri ~ SPEI + (1|Species), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(Tri ~ Sex + (1|Species), 
           data = birds.m, REML = FALSE)
m8 <- lmer(Tri ~ MigStatus + (1|Species), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(Tri ~ Julian + (1|Species), 
           data = birds.m, REML = FALSE)
m10 <- lmer(Tri ~ Event + (1|Species), data = birds.m,
            REML = FALSE)
m11 <- lmer(Tri ~ seconds_since_midnight + (1 | Species), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(Tri ~ Max_Flock_Size + (1|Species), data = birds.m, REML = FALSE)

# null
m13 <- lmer(Tri ~ 1 + (1 | Species), data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m1 <- lmer(Tri ~ Event  + (1|Species),  data = birds.m, 
           REML = FALSE)
m2 <- lmer(Tri ~ Event  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

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

m <- lmer(PC1 ~ Biomass + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lmer(PC1 ~ Diversity + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of diversity
