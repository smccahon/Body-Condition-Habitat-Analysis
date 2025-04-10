#-----------------------------------------#
#  Lesser Yellowlegs Fat Habitat Analysis #
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

# Subset data to only include fall 2021 and fall 2023
leye <- leye[leye$Event %in% c("Fall 2021", "Fall 2023"), ]

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

# Convert Fat to two levels (low and high)
leye <- leye %>%
  mutate(Fat.G = case_when(
    Fat == 0 | Fat == 1 | Fat == 2 ~ "Low",       
    Fat == 3 | Fat == 4 | Fat == 5 ~ "High"))

leye$Fat.G <- factor(leye$Fat.G,
                     levels = c("Low", "High"))

leye$Fat.G_binary <- ifelse(leye$Fat.G == "Low", 0, 1)

#standardize data
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
# correlations > 0.6:
# % ag and ag category
# % ag and dominant crop
# % ag and nearest crop distance
# % veg and dominant crop
# % veg and exposed shoreline
# dominant crop and exposed shoreline
# year and SPEI
# year and Julian
# permanence & SPEI
# permanence & Julian
# SPEI and julian
# julian and detection

# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- glm(Fat.G ~ PercentAg, family = "binomial", data = leye.cs)

# vegetation
m2 <- glm(Fat.G ~ Percent_Total_Veg, family = "binomial", data = leye.cs)

# habitat
m3 <- glm(Fat.G ~ Permanence, family = "binomial", data = leye.cs)
m4 <- glm(Fat.G ~ Percent_Exposed_Shoreline, family = "binomial", data = leye.cs)
m5 <- glm(Fat.G ~ Dist_Closest_Wetland_m, family = "binomial", data = leye.cs)

# weather
m6 <- glm(Fat.G ~ SPEI, family = "binomial", data = leye.cs)

# life history
m7 <- glm(Fat.G ~ Age, family = "binomial", data = leye.cs)
m8 <- glm(Fat.G ~ Sex, family = "binomial", data = leye.cs)

# temporal
m9 <- glm(Fat.G ~ Julian, family = "binomial", data = leye.cs)
m10 <- glm(Fat.G ~ seconds_since_midnight, family = "binomial", data = leye.cs)
m11 <- glm(Fat.G ~ Event, family = "binomial", data = leye.cs)

# flock
m12 <- glm(Fat.G ~ Max_Flock_Size, family = "binomial", data = leye.cs)

# null
m13 <- glm(Fat.G ~ 1, family = "binomial", data = leye.cs)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m13)
confint(m13)

# informative parameters -------------------------------------------------------
# SPEI, julian, event, permanence


# Plot SPEI ----
m <- glm(Fat.G ~ SPEI, data = leye, family = "binomial")

d <- expand.grid(SPEI = seq(min(leye$SPEI),
                                 max(leye$SPEI),
                                 length = 1000)) 

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

leye$Fat.G <- as.numeric(factor(leye$Fat.G, levels = c("Low", "High"))) - 1

ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Drought Condition Index", 
       y = "P(High Fat) in Lesser Yellowlegs") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_jitter(data = leye, aes(x = SPEI, y = Fat.G), 
              size = 2, width = 0.1, height = 0.05)
  
# Plot permanence ----
m <- glm(Fat.G ~ Permanence, data = leye, family = "binomial")

d <- expand.grid(Permanence = unique(leye$Permanence))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

leye$Fat.G <- as.numeric(factor(leye$Fat.G, levels = c("Low", "High"))) - 1

ggplot(d, aes(x = Permanence, y = fit)) +
  geom_point(size = 5, col = "black") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 1) +
  theme_classic() +
  labs(x = "Wetland Permanence", 
       y = "Model Estimated P(High Fat) in Lesser Yellowlegs") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top")


# Plot event ----
m <- glm(Fat.G ~ Event, data = leye, family = "binomial")

d <- expand.grid(Event = unique(leye$Event))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

leye$Fat.G <- as.numeric(factor(leye$Fat.G, levels = c("Low", "High"))) - 1

ggplot(d, aes(x = Event, y = fit)) +
  geom_point(size = 5, col = "black") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 1) +
  theme_classic() +
  labs(x = "Year", 
       y = "Model Estimated P(High Fat) in Lesser Yellowlegs") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  scale_x_discrete(labels = c("Fall 2021" = "2021",
                              "Fall 2023" = "2023"))

