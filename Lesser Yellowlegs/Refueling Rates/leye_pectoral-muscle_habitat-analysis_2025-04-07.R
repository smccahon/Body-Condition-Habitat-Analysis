#----------------------------------------------------#
# Lesser Yellowlegs Pectoral Muscle Habitat Analysis #
#             Non-linear regression                  #
#               Created 2025-04-07                   #
#              Modified 2025-04-07                   #
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

leye$formatted_time <- format(as.POSIXct(leye$seconds_since_midnight, 
                                         origin = "1970-01-01", tz = "UTC"), 
                              "%H:%M")

# subset birds that have pectoral muscle
leye <- leye %>% 
  filter(!is.na(PecSizeBest))

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
m1 <- lm(PecSizeBest ~ Julian, data = leye)
m2 <- lm(PecSizeBest ~ DaysIntoSeason, data = leye)
m3 <- lm(PecSizeBest ~ Event, data = leye)


model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# days into season is best, but does it make the most sense? answer is yes
# however, Julian is less correlated with other variables
plot(leye$Julian, leye$PecSizeBest)
plot(leye$Event, leye$PecSizeBest)
plot(leye$DaysIntoSeason, leye$PecSizeBest)

# Model Selection --------------------------------------------------------------

# ...does time need to be included in every model? answer is no ----------------
m11 <- lm(PecSizeBest ~ seconds_since_midnight +            
            sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
            cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

summary(m11)
confint(m11)

# agriculture
m1 <- lm(PecSizeBest ~ PercentAg, data = leye)

# vegetation
m2 <- lm(PecSizeBest ~ Percent_Total_Veg, data = leye)

# habitat
m3 <- lm(PecSizeBest ~ Permanence, data = leye)
m4 <- lm(PecSizeBest ~ Percent_Exposed_Shoreline, data = leye)
m5 <- lm(PecSizeBest ~ Dist_Closest_Wetland_m, data = leye)

# weather
m6 <- lm(PecSizeBest ~ SPEI, data = leye)

# life history
m7 <- lm(PecSizeBest ~ Age, data = leye)
m8 <- lm(PecSizeBest ~ Sex, data = leye)

# temporal
m9 <- lm(PecSizeBest ~ Julian, data = leye)
m10 <- lm(PecSizeBest ~ Event, data = leye)
m11 <- lm(PecSizeBest ~ seconds_since_midnight +            
            sin(2 * pi * seconds_since_midnight / (24 * 3600)) +             
            cos(2 * pi * seconds_since_midnight / (24 * 3600)), data = leye)

# flock
m12 <- lm(PecSizeBest ~ Max_Flock_Size, data = leye)

# null
m13 <- lm(PecSizeBest ~ 1, data = leye)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
confint(m7) # age is important
confint(m1) # ag is important
confint(m9) # date is important (both Julian & Days Into Season)

# NOTE:
# I selected Julian because it is not correlated with % Ag and is still informative

# Stage 2: Multiple combinations with informative parameters--------------------
# excluding variables with correlations
m1 <- lm(PecSizeBest ~ PercentAg, data = leye)
m2 <- lm(PecSizeBest ~ Age, data = leye)
m3 <- lm(PecSizeBest ~ Julian, data = leye)

m4 <- lm(PecSizeBest ~ Julian + Age, data = leye)
m5 <- lm(PecSizeBest ~ Julian + PercentAg, data = leye)

m6 <- lm(PecSizeBest ~ Age + PercentAg, data = leye)

m7 <- lm(PecSizeBest ~ Julian + PercentAg + Age, data = leye)

model_names <- paste0("m", 1:7)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
confint(m6) # age and percent ag is the top model

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
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Age") +
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
                              color = Age), size = 2) +
  scale_color_manual(values = c("Juvenile" = "palevioletred2", 
                                "Adult" = "steelblue1")) +
  scale_fill_manual(values = c("Juvenile" = "palevioletred2", 
                               "Adult" = "steelblue1"))

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
