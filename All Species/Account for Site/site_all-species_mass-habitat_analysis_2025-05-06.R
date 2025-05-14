#-----------------------------------#
# All Species Mass Habitat Analysis #
#      Linear regression            #
#       Created 2025-05-06          #
#      Modified 2025-05-13          #
#-----------------------------------#

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
                                 "Spring 2023",
                                 "Fall 2023"))

birds$Site <- as.factor(birds$Site)

birds$Permanence <- factor(birds$Permanence,
                           levels = c("Temporary", "Seasonal", "Semipermanent", 
                                      "Permanent"))

birds$MigStatus <- factor(birds$MigStatus,
                           levels = c("Resident", "Migratory"))

# Logarithmic transformation of mass
birds <- birds %>% 
  mutate(LogMass = log(Mass))

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

# Only include birds with mass 
birds <- birds %>%
  filter(!is.na(Mass))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data except for response
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("LogMass"), scale))

# Only include sites with at least three individuals
birds.s <- birds %>% 
  group_by(Site) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data except for response
birds.s.cs <- birds.s %>%
  mutate(across(where(is.numeric) & !matches("LogMass"), scale))


# Test for Correlations--------------------------------------------------------- 

# subset data
sample <- birds.cs[, c("PercentAg",
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
# ag category and dominant crop
# Julian & SPEI

# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lmer(LogMass ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(LogMass ~ DominantCrop + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(LogMass ~ AgCategory + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)

# dominant crop is best (not significant), but % ag is still within 2 delta AICc

# temporal
m1 <- lmer(LogMass ~ Julian + MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(LogMass ~ Julian * MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m3 <- lmer(LogMass ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(LogMass ~ Julian + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model without interaction is more informative, don't include interaction

# transformation for time needed? no
plot(birds$seconds_since_midnight, birds$LogMass)

# interactions between season and capture time? yes, informed null
m1 <- lmer(LogMass ~ Event + seconds_since_midnight + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(LogMass ~ Event * seconds_since_midnight + (1|Species), data = birds.cs, 
           REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model Selection with Informed Null (Stage 1) ---------------------------------

# agriculture
m1 <- lmer(LogMass ~ PercentAg + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(LogMass ~ Percent_Total_Veg + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(LogMass ~ Permanence + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)
m4 <- lmer(LogMass ~ Percent_Exposed_Shoreline + Event * seconds_since_midnight + 
             (1|Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(LogMass ~ Dist_Closest_Wetland_m + Event * seconds_since_midnight + 
             (1|Species), 
           data = birds.cs, REML = FALSE)

# weather
m6 <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(LogMass ~ Sex + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)
m8 <- lmer(LogMass ~ MigStatus + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(LogMass ~ Julian + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)
m10 <- lmer(LogMass ~ Event*seconds_since_midnight + (1|Species), data = birds.cs,
            REML = FALSE)

# flock
m11 <- lmer(LogMass ~ Max_Flock_Size + Event * seconds_since_midnight + 
              (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:11)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m6)
confint(m6) # birds have higher mass in wetter wetlands

summary(m8)
confint(m8) # mass lower in migrants

summary(m9)
confint(m9) # mass decreases with Julian date...
 
summary(m5)
confint(m5) # dist. to closest wetland no effect

summary(m3)
confint(m3) # mass is higher in birds that use seasonal versus temporal wetlands

# Model Selection with Informed Null (Stage 2) ---------------------------------
# Julian and SPEI are correlated -- cannot be in the same model

# univariate + informed null
m1 <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m2 <- lmer(LogMass ~ MigStatus + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m3 <- lmer(LogMass ~ Permanence + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m4 <- lmer(LogMass ~ Julian + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m5 <- lmer(LogMass ~ Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

# additive (2)
m6 <- lmer(LogMass ~ SPEI + MigStatus + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m7 <- lmer(LogMass ~ SPEI + Permanence + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m8 <- lmer(LogMass ~ MigStatus + Permanence + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m9 <- lmer(LogMass ~ MigStatus + Julian + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

m10 <- lmer(LogMass ~ Permanence + Julian + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

# additive (3)
m11 <- lmer(LogMass ~ SPEI + MigStatus + Permanence + Event * seconds_since_midnight + 
              (1|Species), data = birds.cs, REML = FALSE)

m12 <- lmer(LogMass ~ Julian + MigStatus + Permanence + Event * seconds_since_midnight + 
              (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m6)
confint(m6)

# Model Selection with Site as Random Effect -----------------------------------

# agriculture
m1 <- lmer(LogMass ~ PercentAg + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# vegetation
m2 <- lmer(LogMass ~ Percent_Total_Veg + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# habitat
m3 <- lmer(LogMass ~ Permanence + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

m4 <- lmer(LogMass ~ Percent_Exposed_Shoreline + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

m5 <- lmer(LogMass ~ Dist_Closest_Wetland_m + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# weather
m6 <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# life history
m7 <- lmer(LogMass ~ Sex + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

m8 <- lmer(LogMass ~ MigStatus + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

# temporal
m9 <- lmer(LogMass ~ Julian + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

m10 <- lmer(LogMass ~ Event*seconds_since_midnight + 
              (1|Species) + (1|Site), data = birds.s.cs,
            REML = FALSE)

# flock
m11 <- lmer(LogMass ~ Max_Flock_Size + Event * seconds_since_midnight + 
              (1|Species) + (1|Site), data = birds.s.cs, REML = FALSE)

model_names <- paste0("m", 1:11)

models <- mget(model_names)

aictab(models, modnames = model_names) # results don't change and singular fit warning

# Is there enough support to include the random effect? No ----
m1 <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + 
             (1|Species), 
           data = birds.s.cs, REML = FALSE)
m2 <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + 
             (1|Species) + (1|Site), 
           data = birds.s.cs, REML = FALSE)

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



# Plot drought -----------------------------------------------------------------
m <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + (1|Species), 
           data = birds, REML = FALSE)

d <- expand.grid(SPEI = seq(min(birds$SPEI), 
                                 max(birds$SPEI), 
                                 length = 1000),
                 Event = c("Fall 2023"),                    
                 seconds_since_midnight = mean(birds$seconds_since_midnight),
                 Species = unique(birds$Species)) 

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(size = 0.8) + 
   geom_ribbon(aes(ymin = lwr, ymax = upr), 
               alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Drought Index (SPEI)", 
       y = "Log Body Mass (g)",
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
  geom_point(data = birds, aes(x = SPEI, y = LogMass,
                              color = Event), size = 2) +
  scale_color_viridis_d(alpha = 1, begin = 0, end = 0.9)

summary(m)
confint(m)

# what does this really mean?
mean(birds$SPEI) # -0.04267045
sd(birds$SPEI) # 1.221839

sd(birds.cs$SPEI)

# a one unit increase in raw SPEI (+0.818) is associated with an estimated 5.5%
# increase in body mass (g) holding other predictors constant.

# a 1 SD increase in SPEI (or 1.22 unit increase in raw SPEI) is associated 
# with a 6.7% increase in body mass (scaled)

# A 1 SD increase in SPEI (approximately a 1.22-unit increase on the raw scale) 
# was associated with a 6.7% increase in body mass, indicating that birds 
# had higher fat reserves under wetter-than-average conditions.

# The Standardized Precipitation-Evapotranspiration Index (SPEI) was 
# standardized (mean = −0.043, SD = 1.222) prior to model fitting. 
# Thus, model coefficients for SPEI represent the expected change in log 
# body mass associated with a 1 standard deviation increase in SPEI, 
# equivalent to a 1.22-unit increase on the original scale of the index.

# Variables to move forward to chapter 2 ---------------------------------------
# SPEI, Julian, MigStatus, permanence



# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- lmer(LogMass ~ Biomass + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lmer(LogMass ~ Diversity + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of diversity

plot(birds$LogMass, birds$Biomass)



birds <- read.csv("Body_Condition_Habitat_Analysis_2025-03-31.csv")

birds.tri <- birds %>% 
  filter(!is.na(Tri)) %>%
  mutate(across(where(is.numeric) & !matches("Tri"), scale)) %>% 
  filter(!is.na(Biomass))

birds.beta <- birds %>% 
  filter(!is.na(Beta)) %>%
  mutate(across(where(is.numeric) & !matches("Beta"), scale)) %>% 
  filter(!is.na(Biomass))

birds.uric <- birds %>% 
  filter(!is.na(Uric)) %>%
  mutate(across(where(is.numeric) & !matches("Uric"), scale)) %>% 
  filter(!is.na(Biomass))

birds.fat <- birds %>% 
  filter(!is.na(Fat)) %>%
  mutate(across(where(is.numeric) & !matches("Fat"), scale)) %>% 
  filter(!is.na(Biomass))

m <- lmer(Tri ~ Biomass + (1 | Species), data = birds.tri)

summary(m)
confint(m) # no effect of biomass

m <- lmer(Tri ~ Diversity + (1 | Species), data = birds.tri)

summary(m)
confint(m) # no effect of diversity


m <- lmer(Beta ~ Biomass + (1 | Species), data = birds.beta)

summary(m)
confint(m) # no effect of biomass

m <- lmer(Beta ~ Diversity + (1 | Species), data = birds.beta)

summary(m)
confint(m) # no effect of diversity


m <- lmer(Uric ~ Biomass + (1 | Species), data = birds.uric)

summary(m)
confint(m) # no effect of biomass

m <- lmer(Uric ~ Diversity + (1 | Species), data = birds.uric)

summary(m)
confint(m) # no effect of diversity

