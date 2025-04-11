#-----------------------------------------#
# Lesser Yellowlegs Mass Habitat Analysis #
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
# detection & julian


# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lm(Mass ~ PercentAg, data = leye)
m2 <- lm(Mass ~ DominantCrop, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a better predictor

# temporal
m1 <- lm(Mass ~ Julian * Event, data = leye)
m2 <- lm(Mass ~ Julian + Event, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# interaction is much better...but model is unstable. drop interaction

# time and event?
m1 <- lm(Mass ~ seconds_since_midnight * Event, data = leye)
m2 <- lm(Mass ~ seconds_since_midnight + Event, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model with interaction is much for informative

# is date important within season? yes but it's not making sense...too complex
leye.fall <- subset(leye, Event == "Fall 2023")
leye.fall <- subset(leye, Event == "Fall 2021")
leye.spring <- subset(leye, Event == "Spring 2022")

ggplot(leye.spring, aes(x = Julian, y = Mass)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# Model Selection (Stage 1) Event * Time as Informed Null-----------------------

# agriculture
m1 <- lm(Mass ~ PercentAg + Event * seconds_since_midnight, data = leye.cs)

# vegetation
m2 <- lm(Mass ~ Percent_Total_Veg + Event * seconds_since_midnight, data = leye.cs)

# habitat
m3 <- lm(Mass ~ Permanence + Event * seconds_since_midnight, data = leye.cs)
m4 <- lm(Mass ~ Percent_Exposed_Shoreline + Event * seconds_since_midnight, 
         data = leye.cs)
m5 <- lm(Mass ~ Dist_Closest_Wetland_m + Event * seconds_since_midnight, 
         data = leye.cs)

# weather
m6 <- lm(Mass ~ SPEI + Event * seconds_since_midnight, data = leye.cs)

# life history
m7 <- lm(Mass ~ Age + Event * seconds_since_midnight, data = leye.cs)
m8 <- lm(Mass ~ Sex + Event * seconds_since_midnight, data = leye.cs)

# temporal
m9 <- lm(Mass ~ Julian + Event * seconds_since_midnight, data = leye.cs)
m10 <- lm(Mass ~ seconds_since_midnight*Event, data = leye.cs)

# flock
m11 <- lm(Mass ~ Max_Flock_Size + Event * seconds_since_midnight, data = leye.cs)

# informed null
m12 <- lm(Mass ~ Event * seconds_since_midnight, data = leye.cs)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m7) # age in the top model but not significant
confint(m7)

summary(m1) # percent ag not significant
confint(m1)

# CONCLUSION -------------------------------------------------------------------
# Event * Time most informative, age is in top model (but not significant)

# Variables to move forward to chapter 2: age + event * time -------------------

#------------------------------------------------------------------------------#

# PREVIOUS ANALYSIS NOT TAKING EVENT*TIME INTO CONSIDERATION -------------------

# informative parameters 
# Event:Fall 2023 (m11): (B = -1.16, CI: -1.75, -0.57)
# Age (m7): (B = 0.88, CI: 0.32, 1.43)
# Max Flock Size (m12): (B = -0.332, CI: -0.59, -0.07)

# looks like a real relationship -- suggestive of competition?
ggplot(leye, aes(x = Max_Flock_Size, y = Mass)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

# Model Selection (Stage 2) 
m1 <- lm(Mass ~ Age, data = leye.cs)
m2 <- lm(Mass ~ Event, data = leye.cs)
m3 <- lm(Mass ~ Max_Flock_Size, data = leye.cs)

m4 <- lm(Mass ~ Age + Event, data = leye.cs)
m5 <- lm(Mass ~ Age + Max_Flock_Size, data = leye.cs)
m6 <- lm(Mass ~ Max_Flock_Size + Event, data = leye.cs)

m7 <- lm(Mass ~ Age + Event + Max_Flock_Size, data = leye.cs)

model_names <- paste0("m", 1:7)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m5) # age and flock size informative
confint(m7) # only event informative
confint(m2) # event informative
confint(m4) # event informative
confint(m6) # event informative
confint(m1) # age informative

# Variables to move forward to chapter 2
# age, flock size, event

ggplot(leye, aes(x = Event, y = Mass)) + geom_boxplot() +
  theme_classic() +
  labs(x = "Sampling Event", 
       y = "Lesser Yellowlegs Body Mass (g)") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top")

