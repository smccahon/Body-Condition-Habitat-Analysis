#----------------------------------------------------#
# Lesser Yellowlegs Pectoral Muscle Habitat Analysis #
#                linear regression                   #
#               Created 2025-04-07                   #
#              Modified 2025-05-29                   #
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
birds <- read.csv("Body_Condition_Habitat_Analysis_2025-05-29.csv")

# ...make new columns ----

# ...reorder and manipulate relevant factor variables ----
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

# make detection columns a factor
# ** note: did not look at neonics in inverts because there were only two detections
birds$PlasmaDetection <- as.factor(birds$PlasmaDetection)

birds$WaterNeonicDetection <- as.factor(birds$WaterNeonicDetection)

birds$AnyDetection <- as.factor(birds$AnyDetection)

birds$WaterOrInvertDetection <- as.factor(birds$WaterOrInvertDetection)

birds$InvertPesticideDetection <- as.factor(birds$InvertPesticideDetection)

# categorize other factor variables
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

# subset birds that have pectoral muscle
leye <- leye %>% 
  filter(!is.na(PecSizeBest))

# standardize data except for response
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("PecSizeBest"), scale))

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
                   "DaysIntoSeason",
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
# % ag & % exposed shoreline
# % ag & ag category
# % ag & dominant crop
# % ag & nearest crop distance
# % ag and max flock size
# % total veg & dominant crop
# sampling event & detection
# sampling event & SPEI
# sampling event & julian
# permanence & SPEI
# detection & julian
# ag category & nearest crop distance
# dominant crop & max flock size
# days into season & % ag
# days into season & sampling event
# days into season & julian
# days into season & dist. to closest wetland

# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lm(PecSizeBest ~ PercentAg, data = leye)
m2 <- lm(PecSizeBest ~ DominantCrop, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a much better predictor

# temporal
m1 <- lm(PecSizeBest ~ Julian * Event, data = leye)
m2 <- lm(PecSizeBest ~ Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model without interaction is the most informative

# age & % ag interaction
m1 <- lm(PecSizeBest ~ Age * PercentAg, data = leye)
m2 <- lm(PecSizeBest ~ Age + PercentAg, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model without interaction is the most informative

# time transformation necessary? no
plot(leye$seconds_since_midnight, leye$PecSizeBest) # not really any pattern
m1 <- lm(PecSizeBest ~ seconds_since_midnight + sin + cos, data = leye.cs)
m2 <- lm(PecSizeBest ~ seconds_since_midnight, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# time and event interaction? yes
# not necessary to be informed null though...all within 2 delta AICc
m1 <- lm(PecSizeBest ~ seconds_since_midnight * Event, data = leye.cs)
m2 <- lm(PecSizeBest ~ seconds_since_midnight + Event, data = leye.cs)
m3 <- lm(PecSizeBest ~ 1, data = leye.cs)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)


# detection or event? detection is misleading...
plot(leye$Detection, leye$PecSizeBest)
plot(leye$Event, leye$PecSizeBest)
m1 <- lm(PecSizeBest ~ Detection, data = leye)
m2 <- lm(PecSizeBest ~ Event, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m1)

# interaction between ag and drought needed?
m1 <- lm(PecSizeBest ~ PercentAg, data = leye)
m2 <- lm(PecSizeBest ~ SPEI, data = leye)
m3 <- lm(PecSizeBest ~ SPEI + PercentAg, data = leye)
m4 <- lm(PecSizeBest ~ SPEI * PercentAg, data = leye)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model Selection --------------------------------------------------------------

# agriculture
m1 <- lm(PecSizeBest ~ PercentAg, data = leye.cs)

# vegetation
m2 <- lm(PecSizeBest ~ Percent_Total_Veg, data = leye.cs)

# habitat
m3 <- lm(PecSizeBest ~ Permanence, data = leye.cs)
m4 <- lm(PecSizeBest ~ Percent_Exposed_Shoreline, data = leye.cs)
m5 <- lm(PecSizeBest ~ Dist_Closest_Wetland_m, data = leye.cs)

# weather
m6 <- lm(PecSizeBest ~ SPEI, data = leye.cs)

# life history
m7 <- lm(PecSizeBest ~ Age, data = leye.cs)
m8 <- lm(PecSizeBest ~ Sex, data = leye.cs)

# temporal
m9 <- lm(PecSizeBest ~ Julian, data = leye.cs)
m10 <- lm(PecSizeBest ~ seconds_since_midnight, data = leye.cs)
m11 <- lm(PecSizeBest ~ Event, data = leye.cs)

# flock
m12 <- lm(PecSizeBest ~ Max_Flock_Size, data = leye.cs)

# null
m13 <- lm(PecSizeBest ~ 1, data = leye.cs)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
confint(m7) # age is important
confint(m1) # ag is important
confint(m9) # date is important (Julian)

deviance(m3)

# Stage 2: Multiple combinations with informative parameters--------------------
# excluding variables with correlations
m1 <- lm(PecSizeBest ~ PercentAg, data = leye.cs)
m2 <- lm(PecSizeBest ~ Age, data = leye.cs)
m3 <- lm(PecSizeBest ~ Julian, data = leye.cs)

m4 <- lm(PecSizeBest ~ PercentAg + Age, data = leye.cs)
m5 <- lm(PecSizeBest ~ PercentAg + Julian, data = leye.cs)
m6 <- lm(PecSizeBest ~ Age + Julian, data = leye.cs)

m7 <- lm(PecSizeBest ~ PercentAg + Julian + Age, data = leye.cs)

model_names <- paste0("m", 1:7)

models <- mget(model_names)

aictab(models, modnames = model_names)

deviance(m1)


# results
confint(m6) # age and percent ag is the top model

deviance(m3)

# plot top model ---------------------------------------------------------------
m <- lm(PecSizeBest ~ Age + PercentAg, data = leye)

d <- expand.grid(Age = c("Juvenile", "Adult"),
                 PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = (PercentAg), y = fit, color = Age)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Age), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_bw() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Age") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                              color = Age), size = 2) +
  scale_color_manual(values = c("Juvenile" = "#009E73", 
                                "Adult" = "#CC79A7")) +
  scale_fill_manual(values = c("Juvenile" = "#009E73",  
                               "Adult" = "#CC79A7")) 

# coop meeting presentation
d <- expand.grid(Age = c("Adult"),
                 PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = (PercentAg), y = fit, color = Age)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Age), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Age") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "none") +
  scale_color_manual(values = c("Adult" = "firebrick2")) +
  scale_fill_manual(values = c("Adult" = "firebrick2")) 


# juveniles only
leye.juv <- subset(leye, Age == "Juvenile")

m <- lm(PecSizeBest ~ PercentAg, data = leye.juv)

d <- expand.grid(PercentAg = seq(min(leye.juv$PercentAg),
                                 max(leye.juv$PercentAg),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score]))) +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye.juv, aes(x = PercentAg, y = PecSizeBest), size = 2)

confint(m) # for juveniles alone, not significant

# Model Diagnostics ------------------------------------------------------------
m <- lm(PecSizeBest ~ Age + PercentAg, data = leye)
plot(predict(m),rstudent(m)) # great, no concerns

par(mfrow = c(2,2))
plot(m) # perfect; satisfies all assumptions

# CONCLUSIONS & NOTES-----------------------------------------------------------
# First Stage: Age, % Ag, and Date (Julian & Days Into Season) were important
# Second Stage: Tried possible combinations of top covariates, top model was 
# Pectoral Muscle Size ~ Age + % Surrounding Agriculture. This model also was
# the only one within 2 delta AICc and both parameters have an effect. Pectoral
# muscle size and % ag are negatively correlated, and adults have larger pectoral
# muscle sizes.

# COVARIATES TO MOVE FORWARD ---------------------------------------------------
# Age, % Surrounding Ag, and Julian. 

### RECHECK FOR CORRELATIONS

# ARE NEONICS ANY MORE INFORMATIVE? --------------------------------------------
m1 <- lm(PecSizeBest ~ PercentAg, data = leye)
m2 <- lm(PecSizeBest ~ Age, data = leye)
m3 <- lm(PecSizeBest ~ Julian, data = leye)
m4 <- lm(PecSizeBest ~ Detection, data = leye)

m5 <- lm(PecSizeBest ~ PercentAg + Age, data = leye)
m6 <- lm(PecSizeBest ~ PercentAg + Detection, data = leye)

m7 <- lm(PecSizeBest ~ Age + Julian, data = leye)
m8 <- lm(PecSizeBest ~ Age + Julian, data = leye)

m9 <- lm(PecSizeBest ~ Julian + Detection, data = leye)

m10 <- lm(PecSizeBest ~ PercentAg + Age + Detection, data = leye)
m11 <- lm(PecSizeBest ~ Age + Julian + Detection, data = leye)

m12 <- lm(PecSizeBest ~ 1, data = leye)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
confint(m5) # age and % ag important
confint(m6) # detection and % ag important
confint(m10) # none are important

# plot top models --------------------------------------------------------------
m <- lm(PecSizeBest ~ Detection + PercentAg, data = leye)

d <- expand.grid(Detection = c("Detection", "Non-detection"),
                 PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

ggplot(d, aes(x = PercentAg, y = fit, color = Detection)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Detection), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Plasma Neonicotinoid Detection") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                              color = Detection), size = 2) +
  scale_color_manual(values = c("Detection" = "palevioletred2", 
                                "Non-detection" = "steelblue1")) +
  scale_fill_manual(values = c("Detection" = "palevioletred2", 
                               "Non-detection" = "steelblue1"))

# is detection confounded by date? ---------------------------------------------
# absolutely: detections were 100% in spring 2022
ggplot(leye, aes(x = Detection, y = PecSizeBest, color = Event)) + 
  geom_boxplot() + theme_classic()

# is percent ag confounded by date? --------------------------------------------
# yes...
m <- lm(PecSizeBest ~ PercentAg, data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

# days into season
ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Days into Trapping Season") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                              color = DaysIntoSeason), size = 2) +
  scale_color_viridis(alpha = 1, begin = 0, end = 1) +
  geom_text(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                             label = sprintf("%.0f", DaysIntoSeason)), 
            size = 3, vjust = -0.5)

# julian
ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Julian Date") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                              color = Julian), size = 2) +
  scale_color_viridis(alpha = 1, begin = 0, end = 1) +
  geom_text(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                             label = sprintf("%.0f", Julian)), 
            size = 3, vjust = -0.5)

# season
ggplot(d, aes(x = PercentAg, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Sampling Event") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                              color = Event), size = 3) +
  scale_color_viridis_d(alpha = 1, begin = 0, end = 0.8)

# is days into season informative? yes
m <- lm(PecSizeBest ~ DaysIntoSeason, data = leye.fall)

d <- expand.grid(DaysIntoSeason = seq(min(leye.fall$DaysIntoSeason),
                                 max(leye.fall$DaysIntoSeason),
                                 length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

# days into season
ggplot(d, aes(x = DaysIntoSeason, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Number of Days into the Trapping Season", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "% Surrounding Agriculture") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye.fall, aes(x = DaysIntoSeason, y = PecSizeBest,
                              color = PercentAg), size = 2) + 
  scale_color_viridis(alpha = 1, begin = 0, end = 1) +
  geom_text(data = leye.fall, aes(x = DaysIntoSeason, y = PecSizeBest, 
                                  label = sprintf("%.0f", PercentAg)), 
            size = 3, vjust = -0.5)

confint(m)

cor(leye.fall$PercentAg, leye.fall$Julian)

# is Julian date informative? yes
m <- lm(PecSizeBest ~ Julian, data = leye)

d <- expand.grid(Julian = seq(min(leye$Julian),
                                      max(leye$Julian),
                                      length = 1000)) 

d <- cbind(d, predict(m, newdata = d, interval = "confidence"))

# days into season
ggplot(d, aes(x = Julian, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Julian Date", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
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
  geom_point(data = leye, aes(x = Julian, y = PecSizeBest,
                              col = PercentAg), size = 2) +
  scale_color_viridis(alpha = 1, begin = 0, end = 1) +
  geom_text(data = leye, aes(x = Julian, y = PecSizeBest, 
                             label = sprintf("%.0f", PercentAg)), 
            size = 3, vjust = -0.5)

confint(m)

# CONCLUSION AND NOTES WITH NEONICS --------------------------------------------
# mention all birds in spring had detections but pectoral muscle sizes were 
# higher in spring

# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(leye.cs, !is.na(Biomass))

m <- lm(PecSizeBest ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lm(PecSizeBest ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of diversity

#------------------------------------------------------------------------------#

# do neonics explain any further variation in pectoral muscle size? ----
# informative covariates: age, % ag, Julian


# summary statistics----
table(leye$PlasmaDetection) # n: 24, y: 17 (n = 41)
table(leye$WaterNeonicDetection) # n: 32, y: 8 (n = 40)
table(leye$AnyDetection) # n: 12, y: 29 (n = 41)
table(leye$WaterOrInvertDetection) # n: 18, y: 23 (n = 41)
table(leye$InvertPesticideDetection) # n: 15, y: 15 (n = 30)

mean(leye$OverallNeonic, na.rm = TRUE) # 2.39 ug/L
sd(leye$OverallNeonic, na.rm = TRUE) # 10.2 ug/L

# water neonic detection --> neonics not informative (n = 40)
leye.clean.water <- leye.cs[complete.cases(leye.cs[, c("PecSizeBest", "Age", 
                                                        "PercentAg", 
                                                        "WaterNeonicDetection", 
                                                        "Julian")]), ]

m1 <- lm(PecSizeBest ~ Age + PercentAg, data = leye.clean.water)
m2 <- lm(PecSizeBest ~ PercentAg, data = leye.clean.water)
m3 <- lm(PecSizeBest ~ Age, data = leye.clean.water)
m4 <- lm(PecSizeBest ~ WaterNeonicDetection, data = leye.clean.water)
m5 <- lm(PecSizeBest ~ Julian, data = leye.clean.water)
m6 <- lm(PecSizeBest ~ Age + WaterNeonicDetection, data = leye.clean.water)
m7 <- lm(PecSizeBest ~ WaterNeonicDetection + PercentAg, data = leye.clean.water)
m8 <- lm(PecSizeBest ~ Age + PercentAg + WaterNeonicDetection, data = leye.clean.water)
m9 <- lm(PecSizeBest ~ Age + PercentAg + WaterNeonicDetection +
           Julian, data = leye.clean.water)
m10 <- lm(PecSizeBest ~ Age + Julian, data = leye.clean.water)
m11 <- lm(PecSizeBest ~ Age + PercentAg + Julian, data = leye.clean.water)
m12 <- lm(PecSizeBest ~ Julian + PercentAg, data = leye.clean.water)
m13 <- lm(PecSizeBest ~ Julian + WaterNeonicDetection, data = leye.clean.water)
m14 <- lm(PecSizeBest ~ Julian + WaterNeonicDetection + PercentAg, data = leye.clean.water)
m15 <- lm(PecSizeBest ~ Age + Julian + WaterNeonicDetection, data = leye.clean.water)


# invertebrate pesticide detection --> neonics not informative
leye.clean.invert <- leye.cs[complete.cases(leye.cs[, c("PecSizeBest", "Age", 
                                                        "PercentAg", 
                                                        "InvertPesticideDetection", 
                                                        "Julian")]), ]

m1 <- lm(PecSizeBest ~ Age + PercentAg, data = leye.clean.invert)
m2 <- lm(PecSizeBest ~ PercentAg, data = leye.clean.invert)
m3 <- lm(PecSizeBest ~ Age, data = leye.clean.invert)
m4 <- lm(PecSizeBest ~ InvertPesticideDetection, data = leye.clean.invert)
m5 <- lm(PecSizeBest ~ Julian, data = leye.clean.invert)
m6 <- lm(PecSizeBest ~ Age + InvertPesticideDetection, data = leye.clean.invert)
m7 <- lm(PecSizeBest ~ InvertPesticideDetection + PercentAg, data = leye.clean.invert)
m8 <- lm(PecSizeBest ~ Age + PercentAg + InvertPesticideDetection, data = leye.clean.invert)
m9 <- lm(PecSizeBest ~ Age + PercentAg + InvertPesticideDetection +
           Julian, data = leye.clean.invert)
m10 <- lm(PecSizeBest ~ Age + Julian, data = leye.clean.invert)
m11 <- lm(PecSizeBest ~ Age + PercentAg + Julian, data = leye.clean.invert)
m12 <- lm(PecSizeBest ~ Julian + PercentAg, data = leye.clean.invert)
m13 <- lm(PecSizeBest ~ Julian + InvertPesticideDetection, data = leye.clean.invert)
m14 <- lm(PecSizeBest ~ Julian + InvertPesticideDetection + PercentAg, data = leye.clean.invert)
m15 <- lm(PecSizeBest ~ Age + Julian + InvertPesticideDetection, data = leye.clean.invert)


# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
leye.clean.waterorinvert <- leye.cs[complete.cases(leye.cs[, c("PecSizeBest", "Age", 
                                                        "PercentAg", 
                                                        "WaterOrInvertDetection", 
                                                        "Julian")]), ]

m1 <- lm(PecSizeBest ~ Age + PercentAg, data = leye.clean.waterorinvert)
m2 <- lm(PecSizeBest ~ PercentAg, data = leye.clean.waterorinvert)
m3 <- lm(PecSizeBest ~ Age, data = leye.clean.waterorinvert)
m4 <- lm(PecSizeBest ~ WaterOrInvertDetection, data = leye.clean.waterorinvert)
m5 <- lm(PecSizeBest ~ Julian, data = leye.clean.waterorinvert)
m6 <- lm(PecSizeBest ~ Age + WaterOrInvertDetection, data = leye.clean.waterorinvert)
m7 <- lm(PecSizeBest ~ WaterOrInvertDetection + PercentAg, data = leye.clean.waterorinvert)
m8 <- lm(PecSizeBest ~ Age + PercentAg + WaterOrInvertDetection, data = leye.clean.waterorinvert)
m9 <- lm(PecSizeBest ~ Age + PercentAg + WaterOrInvertDetection +
           Julian, data = leye.clean.waterorinvert)
m10 <- lm(PecSizeBest ~ Age + Julian, data = leye.clean.waterorinvert)
m11 <- lm(PecSizeBest ~ Age + PercentAg + Julian, data = leye.clean.waterorinvert)
m12 <- lm(PecSizeBest ~ Julian + PercentAg, data = leye.clean.waterorinvert)
m13 <- lm(PecSizeBest ~ Julian + WaterOrInvertDetection, data = leye.clean.waterorinvert)
m14 <- lm(PecSizeBest ~ Julian + WaterOrInvertDetection + PercentAg, data = leye.clean.waterorinvert)
m15 <- lm(PecSizeBest ~ Age + Julian + WaterOrInvertDetection, data = leye.clean.waterorinvert)


# shorebird plasma detection --> positive relationship with pecs...influenced by season
leye.clean.plasma <- leye.cs[complete.cases(leye.cs[, c("PecSizeBest", "Age", 
                                                               "PercentAg", 
                                                               "PlasmaDetection", 
                                                               "Julian",
                                                        "Event")]), ]

m1 <- lm(PecSizeBest ~ Age + PercentAg, data = leye.clean.plasma)
m2 <- lm(PecSizeBest ~ PercentAg, data = leye.clean.plasma)
m3 <- lm(PecSizeBest ~ Age, data = leye.clean.plasma)
m4 <- lm(PecSizeBest ~ PlasmaDetection, data = leye.clean.plasma)
m5 <- lm(PecSizeBest ~ Julian, data = leye.clean.plasma)
m6 <- lm(PecSizeBest ~ Age + PlasmaDetection, data = leye.clean.plasma)
m7 <- lm(PecSizeBest ~ PlasmaDetection + PercentAg, data = leye.clean.plasma)
m8 <- lm(PecSizeBest ~ Age + PercentAg + PlasmaDetection, data = leye.clean.plasma)
m9 <- lm(PecSizeBest ~ Age + PercentAg + PlasmaDetection +
           Julian, data = leye.clean.plasma)
m10 <- lm(PecSizeBest ~ Age + Julian, data = leye.clean.plasma)
m11 <- lm(PecSizeBest ~ Age + PercentAg + Julian, data = leye.clean.plasma)
m12 <- lm(PecSizeBest ~ Julian + PercentAg, data = leye.clean.plasma)
m13 <- lm(PecSizeBest ~ Julian + PlasmaDetection, data = leye.clean.plasma)
m14 <- lm(PecSizeBest ~ Julian + PlasmaDetection + PercentAg, data = leye.clean.plasma)
m15 <- lm(PecSizeBest ~ Age + Julian + PlasmaDetection, data = leye.clean.plasma)

# correlations & explanations for opposite result
plot(leye.clean.plasma$PlasmaDetection, leye.clean.plasma$PecSizeBest)
plot(leye.clean.plasma$Julian, leye.clean.plasma$PecSizeBest)
plot(leye.clean.plasma$Event, leye.clean.plasma$PecSizeBest)

table(leye.clean.plasma$PlasmaDetection, leye.clean.plasma$Event)
leye$PlasmaDetection <- as.numeric(leye$PlasmaDetection)
leye$Event <- as.numeric(leye$Event)

cor(leye$PlasmaDetection, leye$Event)
cor(leye$PlasmaDetection, leye$Julian) # highly influenced by date

leye$Event <- as.factor(leye$Event)

leye.cs$Event <- relevel(leye.cs$Event, ref = "Fall 2023")
m <- lm(PecSizeBest ~ Event, data= leye.cs) # not significant

summary(m)

t.test(data = leye, PecSizeBest~Event) # marginal difference


# any detection (plasma or environmental) --> neonics not informative
leye.clean.all <- leye.cs[complete.cases(leye.cs[, c("PecSizeBest", "Age", 
                                                        "PercentAg", 
                                                        "PlasmaDetection", 
                                                        "Julian")]), ]

m1 <- lm(PecSizeBest ~ Age + PercentAg, data = leye.clean.all)
m2 <- lm(PecSizeBest ~ PercentAg, data = leye.clean.all)
m3 <- lm(PecSizeBest ~ Age, data = leye.clean.all)
m4 <- lm(PecSizeBest ~ AnyDetection, data = leye.clean.all)
m5 <- lm(PecSizeBest ~ Julian, data = leye.clean.all)
m6 <- lm(PecSizeBest ~ Age + AnyDetection, data = leye.clean.all)
m7 <- lm(PecSizeBest ~ AnyDetection + PercentAg, data = leye.clean.all)
m8 <- lm(PecSizeBest ~ Age + PercentAg + AnyDetection, data = leye.clean.all)
m9 <- lm(PecSizeBest ~ Age + PercentAg + AnyDetection +
           Julian, data = leye.clean.all)
m10 <- lm(PecSizeBest ~ Age + Julian, data = leye.clean.all)
m11 <- lm(PecSizeBest ~ Age + PercentAg + Julian, data = leye.clean.all)
m12 <- lm(PecSizeBest ~ Julian + PercentAg, data = leye.clean.all)
m13 <- lm(PecSizeBest ~ Julian + AnyDetection, data = leye.clean.all)
m14 <- lm(PecSizeBest ~ Julian + AnyDetection + PercentAg, data = leye.clean.all)
m15 <- lm(PecSizeBest ~ Age + Julian + AnyDetection, data = leye.clean.all)

### ...AIC 
models <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15)
model.sel(models)

# model summaries:
summary(m8)
confint(m8)

