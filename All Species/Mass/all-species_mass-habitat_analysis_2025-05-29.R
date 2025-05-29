#-----------------------------------#
# All Species Mass Habitat Analysis #
#      Linear regression            #
#       Created 2025-04-11          #
#      Modified 2025-05-29          #
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
birds <- read.csv("Body_Condition_Habitat_Analysis_2025-05-29.csv")

# make detection columns a factor
# ** note: did not look at neonics in inverts because there were only two detections
birds$PlasmaDetection <- as.factor(birds$PlasmaDetection)

birds$WaterNeonicDetection <- as.factor(birds$WaterNeonicDetection)

birds$AnyDetection <- as.factor(birds$AnyDetection)

birds$WaterOrInvertDetection <- as.factor(birds$WaterOrInvertDetection)

birds$InvertPesticideDetection <- as.factor(birds$InvertPesticideDetection)

# ...reorder and manipulate relevant factor variables ----
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

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
  mutate(across(where(is.numeric) & !matches("Mass"), scale))


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
                   "PlasmaDetection",
                   "seconds_since_midnight",
                   "Site",
                   "AgCategory",
                   "SPEI",
                   "Julian",
                   "DominantCrop",
                   "NearestCropDistance_m",
                   "Dist_Closest_Wetland_m",
                   "Max_Flock_Size",
                   "MigStatus",
                   "InvertPesticideDetection",
                   "WaterOrInvertDetection",
                   "AnyDetection",
                   "WaterNeonicDetection"
)]

# convert categorical to numeric for correlation matrix
sample$Event <- as.numeric(sample$Event)
sample$Permanence <- as.numeric(sample$Permanence)
sample$AgCategory <- as.numeric(sample$AgCategory)
sample$DominantCrop <- as.numeric(sample$DominantCrop)
sample$Sex <- as.numeric(sample$Sex)
sample$PlasmaDetection <- as.numeric(sample$PlasmaDetection)
sample$WaterNeonicDetection <- as.numeric(sample$WaterNeonicDetection)
sample$AnyDetection <- as.numeric(sample$AnyDetection)
sample$InvertPesticideDetection <- as.numeric(sample$InvertPesticideDetection)
sample$WaterOrInvertDetection <- as.numeric(sample$WaterOrInvertDetection)
sample$Site <- as.numeric(sample$Site)
sample$MigStatus <- as.numeric(sample$MigStatus)

cor(sample, use = "pairwise.complete.obs")

# correlations > 0.6:
# % ag and ag category
# % ag and dominant crop
# ag category and dominant crop
# Julian & SPEI
# invert pesticide detection & site
# water neonic & event (-0.628)


# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lmer(LogMass ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(LogMass ~ DominantCrop + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(LogMass ~ AgCategory + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)

# interaction between agriculture and SPEI: does not improve model fit
m1 <- lmer(LogMass ~ PercentAg +  Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(LogMass ~ PercentAg * SPEI + Event * seconds_since_midnight + 
             (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

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

# Plot drought -----------------------------------------------------------------
m <- lmer(LogMass ~ SPEI + Event * seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

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


# do neonics explain any further variation of body mass than event * time? ----
# informative covariates: SPEI & event*time (informed null)


# summary statistics----
table(birds$PlasmaDetection) # n: 109, y: 60 (n = 169)
table(birds$WaterNeonicDetection) # n: 149, y: 25 (n = 174)
table(birds$AnyDetection) # n: 52, y: 124 (n = 176)
table(birds$WaterOrInvertDetection) # n: 84, y: 92 (n = 176)
table(birds$InvertPesticideDetection) # n: 55, y: 67 (n = 122)

mean(birds$OverallNeonic, na.rm = TRUE) # 8.76 ug/L
sd(birds$OverallNeonic, na.rm = TRUE) # 79.4 ug/L

# water neonic detection --> neonics not informative
birds.clean.water <- birds.cs[!is.na(birds.cs$WaterNeonicDetection), ] #n = 174

m1 <- lmer(LogMass ~ Event * seconds_since_midnight + (1|Species), data = birds.clean.water,
           REML = FALSE)
m2 <- lmer(LogMass ~ Event * seconds_since_midnight + WaterNeonicDetection + 
             (1|Species), data = birds.clean.water, REML = FALSE)
m3 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species), 
           data = birds.clean.water, REML = FALSE)
m4 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species) +
             WaterNeonicDetection, data = birds.clean.water, REML = FALSE)

# invertebrate pesticide detection --> neonics not informative
birds.clean.invert <- birds.cs[!is.na(birds.cs$InvertPesticideDetection), ] #n = 122

m1 <- lmer(LogMass ~ Event * seconds_since_midnight + (1|Species), data = birds.clean.water,
           REML = FALSE)
m2 <- lmer(LogMass ~ Event * seconds_since_midnight + InvertPesticideDetection + 
             (1|Species), data = birds.clean.water, REML = FALSE)
m3 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species), 
           data = birds.clean.water, REML = FALSE)
m4 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species) +
             InvertPesticideDetection, data = birds.clean.water, REML = FALSE)

# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
birds.clean.waterorinvert <- birds.cs[!is.na(birds.cs$WaterOrInvertDetection), ] #n = 176

m1 <- lmer(LogMass ~ Event * seconds_since_midnight + (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)
m2 <- lmer(LogMass ~ Event * seconds_since_midnight + WaterOrInvertDetection + 
             (1|Species), data = birds.clean.waterorinvert, REML = FALSE)
m3 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species), 
           data = birds.clean.waterorinvert, REML = FALSE)
m4 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species) +
             WaterOrInvertDetection, data = birds.clean.waterorinvert, REML = FALSE)


# shorebird plasma detection --> neonics not informative
birds.clean.plasma <- birds.cs[!is.na(birds.cs$PlasmaDetection), ] #n = 169

m1 <- lmer(LogMass ~ Event * seconds_since_midnight + (1|Species), data = birds.clean.plasma,
           REML = FALSE)
m2 <- lmer(LogMass ~ Event * seconds_since_midnight + PlasmaDetection + 
             (1|Species), data = birds.clean.plasma, REML = FALSE)
m3 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species), 
           data = birds.clean.plasma, REML = FALSE)
m4 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species) +
             PlasmaDetection, data = birds.clean.plasma, REML = FALSE)

# any detection (plasma or environmental) --> neonics not informative
birds.clean.any <- birds.cs[!is.na(birds.cs$AnyDetection), ] #n = 176

m1 <- lmer(LogMass ~ Event * seconds_since_midnight + (1|Species), data = birds.clean.any,
           REML = FALSE)
m2 <- lmer(LogMass ~ Event * seconds_since_midnight + AnyDetection + 
             (1|Species), data = birds.clean.any, REML = FALSE)
m3 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species), 
           data = birds.clean.any, REML = FALSE)
m4 <- lmer(LogMass ~ Event * seconds_since_midnight + SPEI + (1|Species) +
             AnyDetection, data = birds.clean.any, REML = FALSE)

### ...AIC 
models <- list(m1, m2, m3, m4)
model.sel(models)

# model summaries:
summary(m4)
confint(m4)

