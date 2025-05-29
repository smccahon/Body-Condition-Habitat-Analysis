#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#               linear regression                    #
#               Created 2025-04-03                   #
#              Modified 2025-04-11                   #
#----------------------------------------------------#

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

# ...standardize data except for response ----
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("PC1"), scale))


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

# leye$PercentAg <- scale(leye$PercentAg)
cor(leye$PercentAg, leye$sin)
cor(leye$PercentAg, leye$cos)

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


# which correlated covariates should I drop? -----------------------------------
# REDO THIS WHOLE SECTION FOR LEYE if you're doing stage 2 ----

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

# Plot PC1 ~ PercentAg + Time --------------------------------------------------
leye$formatted_time <- format(as.POSIXct(leye$seconds_since_midnight, 
                                         origin = "1970-01-01", tz = "UTC"), 
                              "%H:%M")




m <- lm(PC1 ~ PercentAg + seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length.out = 1000),
                 seconds_since_midnight = mean(leye$seconds_since_midnight))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "% Surrounding Ag within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1, 
                              col = formatted_time), size = 3) +
  scale_color_viridis_d(alpha = 1, begin = 1, end = 0) +
  geom_text(data = leye, aes(x = PercentAg, y = PC1, 
                             label = Site), 
            size = 3, vjust = -0.5) 


# % ag and time are correlated...which is a better fit? ------------------------

m1 <- lm(PC1 ~ PercentAg, data = leye.cs)
m2 <- lm(PC1 ~ seconds_since_midnight + sin + cos, data = leye.cs)
m3 <- lm(PC1 ~ 1, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)


confint(m1)

# exclude data prior to 08:00 and assess relationships -------------------------

# format subset data
leye$formatted_time <- format(as.POSIXct(leye$seconds_since_midnight, 
                                         origin = "1970-01-01", tz = "UTC"), 
                              "%H:%M")


leye_late <- subset(leye, formatted_time > "08:00")
leye_late.cs <- leye_late %>%
  mutate(across(where(is.numeric), scale))

m <- lm(PC1 ~ PercentAg,
        data = leye_late.cs)

m <- lm(PC1 ~ PercentAg + seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye_late)

d <- expand.grid(seconds_since_midnight = seq(min(leye_late$seconds_since_midnight), 
                                 max(leye_late$seconds_since_midnight), 
                                 length.out = 1000),
                 PercentAg = mean(leye_late$PercentAg))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

# time
ggplot(d, aes(x = seconds_since_midnight, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Capture Time", 
       y = "Lesser Yellowlegs Fattening Index",
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
  geom_point(data = leye_late, aes(x = seconds_since_midnight, y = PC1, 
                              col = PercentAg), size = 3) +
  scale_x_time(labels = scales::time_format("%H:%M"),
               breaks = seq(0, 86400, by = 7200)) +
  scale_color_viridis(alpha = 1, begin = 1, end = 0) +
  # geom_text(data = leye_late, aes(x = seconds_since_midnight, y = PC1,
  #                                 label = Site))
   geom_text(data = leye_late, aes(x = seconds_since_midnight, y = PC1, 
                              label = sprintf("%.0f", PercentAg)), 
             size = 3, vjust = -0.5) 

confint(m)


# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(leye.cs, !is.na(Biomass))

m <- lm(PC1 ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lm(PC1 ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of diversity




# PLOT RESIDUALS AGAINST AG ----
m <- lm(PC1 ~ seconds_since_midnight + sin + cos, data = leye.cs)

leye.cs$residuals <- resid(m)

ggplot(leye.cs, aes(x = PercentAg, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 0.75) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_classic() +
  labs(x = "% Surrounding Agriculture",
       y = "Residuals of Fattening Index ~ Time",
       title = "Unexplained Variation in Fattening Index vs. % Agriculture (Lesser Yellowlegs)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

summary(lm(residuals ~ PercentAg, data = leye.cs)) # Adj. R2 = -0.0241



# are random effects necessary? ----
leye.m <- leye.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()

m1 <- lmer(PC1 ~ PercentAg + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)

# vegetation
m2 <- lmer(PC1 ~ Percent_Total_Veg + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)

# habitat
m3 <- lmer(PC1 ~ Permanence + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)
m4 <- lmer(PC1 ~ Percent_Exposed_Shoreline + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)
m5 <- lmer(PC1 ~ Dist_Closest_Wetland_m + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)

# weather
m6 <- lmer(PC1 ~ SPEI + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)

# life history
m7 <- lmer(PC1 ~ Age + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)
m8 <- lmer(PC1 ~ Sex + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)

# temporal
m9 <- lmer(PC1 ~ Julian + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)
m10 <- lmer(PC1 ~ Event + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)
m11 <- lmer(PC1 ~ seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)

# flock
m12 <- lmer(PC1 ~ Max_Flock_Size + seconds_since_midnight + sin + cos +(1|Site), REML = FALSE, data = leye.m)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names) # results did not change

# run reduced dataset without site ----

m1 <- lm(PC1 ~ PercentAg + seconds_since_midnight + sin + cos, data = leye.m)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg + seconds_since_midnight +
           sin +
           cos, data = leye.m)

# habitat
m3 <- lm(PC1 ~ Permanence + seconds_since_midnight + 
           sin +
           cos, data = leye.m)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline + seconds_since_midnight + 
           sin +             
           cos, data = leye.m)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m + seconds_since_midnight +  
           sin + 
           cos, data = leye.m)

# weather
m6 <- lm(PC1 ~ SPEI + seconds_since_midnight + 
           sin + 
           cos, data = leye.m)

# life history
m7 <- lm(PC1 ~ Age + seconds_since_midnight + 
           sin + 
           cos, data = leye.m)
m8 <- lm(PC1 ~ Sex + seconds_since_midnight + 
           sin + 
           cos, data = leye.m)

# temporal
m9 <- lm(PC1 ~ Julian + seconds_since_midnight +            
           sin +             
           cos, data = leye.m)
m10 <- lm(PC1 ~ Event + seconds_since_midnight +            
            sin +             
            cos, data = leye.m)
m11 <- lm(PC1 ~ seconds_since_midnight +            
            sin +             
            cos, data = leye.m)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size + seconds_since_midnight +           
            sin +             
            cos, data = leye.m)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)


confint(m1)


# Is there enough support to include the random effect? No ----
m1 <- lm(PC1 ~ seconds_since_midnight + sin + cos, data = leye.m)
m2 <- lmer(PC1 ~ seconds_since_midnight + sin + cos + (1|Site),
           data = leye.m, REML = FALSE)


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
