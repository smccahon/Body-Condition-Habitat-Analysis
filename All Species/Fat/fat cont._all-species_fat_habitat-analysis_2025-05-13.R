#----------------------------------------------#
#      All Species Fat Habitat Analysis        #
#          Created 2025-04-11                  #
#          Modified 2025-05-29                 #
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
  filter(!is.na(Fat))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data except for response
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("Fat"), scale))


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

cor(sample)
# correlations > 0.6:
# % ag and ag category
# % ag and nearest crop distance
# % exposed shoreline and % veg
# Event & SPEI
# % exposed shoreline and permanence
# permanence and SPEI
# permanence and julian
# permanence and exposed shoreline
# exposed shoreline and SPEI
# exposed shoreline and Julian
# SPEI and julian


# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lmer(Fat ~ PercentAg + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ DominantCrop + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Fat ~ AgCategory + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m3)

# Ag Category the best 

# temporal --> yes include interaction
m1 <- lmer(Fat ~ Julian + (1 | Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ Julian * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Fat ~ MigStatus + (1 | Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Fat ~ Julian + MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m2)

# interaction needed between ag and SPEI?
m1 <- lmer(Fat ~ PercentAg + Julian * MigStatus + (1 | Species), 
           data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ SPEI + Julian * MigStatus + (1 | Species), 
           data = birds.cs, REML = FALSE)
m3 <- lmer(Fat ~ PercentAg + SPEI + Julian * MigStatus + (1 | Species), 
           data = birds.cs, REML = FALSE)
m4 <- lmer(Fat ~ PercentAg * SPEI + Julian * MigStatus + (1 | Species), 
           data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)



# transformation for time needed? no
plot(birds$seconds_since_midnight, birds$Fat)

# interactions between season and capture time? yes, keep in models
m1 <- lmer(Fat ~ Event + seconds_since_midnight + (1 | Species), 
            data = birds.cs, REML = FALSE)
m2 <- lmer(Fat ~ Event * seconds_since_midnight + (1 | Species), 
            data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m2) # interaction significant, keep in model

# Model Selection with interactions (Stage 1) ----------------------------------

# agriculture
m1 <- lmer(Fat ~ PercentAg + (1 | Species), 
            data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(Fat ~ Percent_Total_Veg + (1 | Species), 
            data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(Fat ~ Permanence + (1 | Species), 
            data = birds.cs, REML = FALSE)
m4 <- lmer(Fat ~ Percent_Exposed_Shoreline + (1 | Species), 
            data = birds.cs, REML = FALSE)
m5 <- lmer(Fat ~ Dist_Closest_Wetland_m + (1 | Species), 
            data = birds.cs, REML = FALSE)

# weather
m6 <- lmer(Fat ~ SPEI + (1 | Species), data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(Fat ~ Sex + (1 | Species), data = birds.cs, REML = FALSE)
m8 <- lmer(Fat ~ MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(Fat ~ Julian + (1 | Species), 
            data = birds.cs, REML = FALSE)
m10 <- lmer(Fat ~ seconds_since_midnight + (1 | Species), 
             data = birds.cs, REML = FALSE)
m11 <- lmer(Fat ~ Event + (1 | Species), data = birds.cs, REML = FALSE)

# interactions
m12 <- lmer(Fat ~ Event * seconds_since_midnight + (1 | Species),
            data = birds.cs, REML = FALSE)

m13 <- lmer(Fat ~ Julian * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

# flock
m14 <- lmer(Fat ~ Max_Flock_Size + (1 | Species), 
             data = birds.cs, REML = FALSE)

# null
m15 <- lmer(Fat ~ (1 | Species), data = birds.cs, REML = FALSE)


model_names <- paste0("m", 1:15)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results (informative parameters across models)

# Migratory Status * Julian VERY Significant 
# Try as informed null with everything but permanence, SPEI, and exp. shoreline
# same results...no informed null then. nothing but julian * Migratory Status important ----

# Model Selection with Julian * Migratory Status as Informed Null (Stage 1) ----

# agriculture
m1 <- lmer(Fat ~ PercentAg + Julian * MigStatus + (1 | Species), 
           data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(Fat ~ Percent_Total_Veg + Julian * MigStatus +  (1 | Species), 
           data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(Fat ~ Permanence + (1 | Species), 
           data = birds.cs, REML = FALSE)
m4 <- lmer(Fat ~ Percent_Exposed_Shoreline + (1 | Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(Fat ~ Dist_Closest_Wetland_m + Julian * MigStatus + (1 | Species), 
           data = birds.cs, REML = FALSE)

# weather
m6 <- lmer(Fat ~ SPEI + (1 | Species), data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(Fat ~ Sex + Julian * MigStatus +  (1 | Species), data = birds.cs, REML = FALSE)
m8 <- lmer(Fat ~ MigStatus * Julian + (1 | Species), data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(Fat ~ seconds_since_midnight + Julian * MigStatus + (1 | Species), 
            data = birds.cs, REML = FALSE)
m10 <- lmer(Fat ~ Event + Julian * MigStatus + (1 | Species), data = birds.cs, REML = FALSE)

# interactions
m11 <- lmer(Fat ~ Event * seconds_since_midnight + Julian * MigStatus + (1 | Species),
            data = birds.cs, REML = FALSE)

# flock
m12 <- lmer(Fat ~ Max_Flock_Size + Julian * MigStatus + (1 | Species), 
            data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results (informative parameters across models)

# Model Selection with Site as Random Effect -----------------------------------

# agriculture
m1 <- lmer(Fat ~ PercentAg + (1 | Species) + (1 | Site), 
           data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(Fat ~ Percent_Total_Veg + (1 | Species) + (1 | Site), 
           data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(Fat ~ Permanence + (1 | Species) + (1 | Site), 
           data = birds.cs, REML = FALSE)
m4 <- lmer(Fat ~ Percent_Exposed_Shoreline + (1 | Species) + (1 | Site), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(Fat ~ Dist_Closest_Wetland_m + (1 | Species) + (1 | Site), 
           data = birds.cs, REML = FALSE)

# weather
m6 <- lmer(Fat ~ SPEI + (1 | Species) + (1 | Site), data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(Fat ~ Sex + (1 | Species) + (1 | Site), data = birds.cs, REML = FALSE)
m8 <- lmer(Fat ~ MigStatus + (1 | Species) + (1 | Site), data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(Fat ~ Julian + (1 | Species) + (1 | Site), 
           data = birds.cs, REML = FALSE)
m10 <- lmer(Fat ~ seconds_since_midnight + (1 | Species) + (1 | Site), 
            data = birds.cs, REML = FALSE)
m11 <- lmer(Fat ~ Event + (1 | Species) + (1 | Site), data = birds.cs, REML = FALSE)

# interactions
m12 <- lmer(Fat ~ Event * seconds_since_midnight + (1 | Species) + (1 | Site),
            data = birds.cs, REML = FALSE)

m13 <- lmer(Fat ~ Julian * MigStatus + (1 | Species) + (1 | Site), data = birds.cs, REML = FALSE)

# flock
m14 <- lmer(Fat ~ Max_Flock_Size + (1 | Species) + (1 | Site), 
            data = birds.cs, REML = FALSE)

# null
m15 <- lmer(Fat ~ (1 | Species) + (1 | Site), data = birds.cs, REML = FALSE)


model_names <- paste0("m", 1:15)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m1 <- lmer(Fat ~ MigStatus * Julian + (1|Species), 
           data = birds.s.cs, REML = FALSE)

m2 <- lmer(Fat ~ MigStatus * Julian + (1|Species) + (1|Site), 
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



# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- lmer(Fat ~ Biomass + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lmer(Fat ~ Diversity + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of diversity



# SPEI with IWLS ----
m <- lmer(Fat ~ SPEI + (1 | Species), data = birds.cs)

birds.cs <- birds.cs %>% 
  group_by(Species) %>% 
  mutate(weight = 1 / (var(Fat) + 0.0001)) %>%  # Avoid division by zero
  ungroup()


m.wls <- lmer(Fat ~ SPEI + (1|Species),
              data = birds.cs,
              REML = FALSE,
              weights = weight)

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


# do neonics explain any further variation of fat than event * time? ----
# informative covariates: Julian * Migratory Status, Event * Time (no informed null)


# summary statistics----
table(birds$PlasmaDetection) # n: 110, y: 60 (n = 170)
table(birds$WaterNeonicDetection) # n: 150, y: 25 (n = 175)
table(birds$AnyDetection) # n: 53, y: 124 (n = 177)
table(birds$WaterOrInvertDetection) # n: 85, y: 92 (n = 177)
table(birds$InvertPesticideDetection) # n: 56, y: 67 (n = 123)

mean(birds$OverallNeonic, na.rm = TRUE) # 8.76 ug/L
sd(birds$OverallNeonic, na.rm = TRUE) # 79.4 ug/L

# water neonic detection --> neonics not informative
birds.clean.water <- birds.cs[!is.na(birds.cs$WaterNeonicDetection), ] #n = 175

m1 <- lmer(Fat ~ Event * seconds_since_midnight + (1|Species), 
           data = birds.clean.water,
           REML = FALSE)

m2 <- lmer(Fat ~ Julian * MigStatus + (1|Species), data = birds.clean.water,
           REML = FALSE)

m3 <- lmer(Fat ~ Event * seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.water, REML = FALSE)

m4 <- lmer(Fat ~ Event * seconds_since_midnight + WaterNeonicDetection +
             (1|Species), data = birds.clean.water, REML = FALSE)

m5 <- lmer(Fat ~ Julian * MigStatus + WaterNeonicDetection +
             (1|Species), data = birds.clean.water, REML = FALSE)

m6 <- lmer(Fat ~ Event * seconds_since_midnight + 
             Julian * MigStatus + WaterNeonicDetection +
             (1|Species), data = birds.clean.water, REML = FALSE)

m7 <- lmer(Fat ~ WaterNeonicDetection + (1|Species), data = birds.clean.water,
           REML = FALSE)

# invertebrate pesticide detection --> neonics not informative
birds.clean.invert <- birds.cs[!is.na(birds.cs$InvertPesticideDetection), ] #n = 123

m1 <- lmer(Fat ~ Event * seconds_since_midnight + (1|Species), 
           data = birds.clean.invert,
           REML = FALSE)

m2 <- lmer(Fat ~ Julian * MigStatus + (1|Species), data = birds.clean.invert,
           REML = FALSE)

m3 <- lmer(Fat ~ Event * seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.invert, REML = FALSE)

m4 <- lmer(Fat ~ Event * seconds_since_midnight + InvertPesticideDetection +
             (1|Species), data = birds.clean.invert, REML = FALSE)

m5 <- lmer(Fat ~ Julian * MigStatus + InvertPesticideDetection +
             (1|Species), data = birds.clean.invert, REML = FALSE)

m6 <- lmer(Fat ~ Event * seconds_since_midnight + 
             Julian * MigStatus + InvertPesticideDetection +
             (1|Species), data = birds.clean.invert, REML = FALSE)

m7 <- lmer(Fat ~ InvertPesticideDetection + (1|Species), data = birds.clean.invert,
           REML = FALSE)

# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
birds.clean.waterorinvert <- birds.cs[!is.na(birds.cs$WaterOrInvertDetection), ] #n = 177

m1 <- lmer(Fat ~ Event * seconds_since_midnight + (1|Species), 
           data = birds.clean.waterorinvert,
           REML = FALSE)

m2 <- lmer(Fat ~ Julian * MigStatus + (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m3 <- lmer(Fat ~ Event * seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.waterorinvert, REML = FALSE)

m4 <- lmer(Fat ~ Event * seconds_since_midnight + WaterOrInvertDetection +
             (1|Species), data = birds.clean.waterorinvert, REML = FALSE)

m5 <- lmer(Fat ~ Julian * MigStatus + WaterOrInvertDetection +
             (1|Species), data = birds.clean.waterorinvert, REML = FALSE)

m6 <- lmer(Fat ~ Event * seconds_since_midnight + 
             Julian * MigStatus + WaterOrInvertDetection +
             (1|Species), data = birds.clean.waterorinvert, REML = FALSE)

m7 <- lmer(Fat ~ WaterOrInvertDetection + (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

# shorebird plasma detection --> neonics not informative
birds.clean.plasma <- birds.cs[!is.na(birds.cs$PlasmaDetection), ] #n = 170

m1 <- lmer(Fat ~ Event * seconds_since_midnight + (1|Species), 
           data = birds.clean.plasma,
           REML = FALSE)

m2 <- lmer(Fat ~ Julian * MigStatus + (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m3 <- lmer(Fat ~ Event * seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.plasma, REML = FALSE)

m4 <- lmer(Fat ~ Event * seconds_since_midnight + PlasmaDetection +
             (1|Species), data = birds.clean.plasma, REML = FALSE)

m5 <- lmer(Fat ~ Julian * MigStatus + PlasmaDetection +
             (1|Species), data = birds.clean.plasma, REML = FALSE)

m6 <- lmer(Fat ~ Event * seconds_since_midnight + 
             Julian * MigStatus + PlasmaDetection +
             (1|Species), data = birds.clean.plasma, REML = FALSE)

m7 <- lmer(Fat ~ PlasmaDetection + (1|Species), data = birds.clean.plasma,
           REML = FALSE)

# any detection (plasma or environmental) --> neonics not informative
birds.clean.any <- birds.cs[!is.na(birds.cs$AnyDetection), ] #n = 177


m1 <- lmer(Fat ~ Event * seconds_since_midnight + (1|Species), 
           data = birds.clean.any,
           REML = FALSE)

m2 <- lmer(Fat ~ Julian * MigStatus + (1|Species), data = birds.clean.any,
           REML = FALSE)

m3 <- lmer(Fat ~ Event * seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.any, REML = FALSE)

m4 <- lmer(Fat ~ Event * seconds_since_midnight + AnyDetection +
             (1|Species), data = birds.clean.any, REML = FALSE)

m5 <- lmer(Fat ~ Julian * MigStatus + AnyDetection +
             (1|Species), data = birds.clean.any, REML = FALSE)

m6 <- lmer(Fat ~ Event * seconds_since_midnight + 
             Julian * MigStatus + AnyDetection +
             (1|Species), data = birds.clean.any, REML = FALSE)

m7 <- lmer(Fat ~ AnyDetection + (1|Species), data = birds.clean.any,
           REML = FALSE)


### ...AIC 
models <- list(m1, m2, m3, m4, m5, m6, m7)
model.sel(models)

# model summaries:
summary(m6)
confint(m6)


