#----------------------------------------------#
#    All Species Uric Acid Habitat Analysis    #
#            Created 2025-05-05                #
#           Modified 2025-05-30                #
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


# filter birds that only contain metabolite information (n = 85)
birds <- birds %>% 
  filter(!is.na(Uric))


# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("Uric"), scale))

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
# % ag and dominant crop
# % veg and % shoreline
# event and julian
# SPEI and julian
# invert pesticide detection & site

# temporal
m1 <- lmer(Uric ~ Julian + MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(Uric ~ Julian * MigStatus + (1|Species), data = birds.cs, 
           REML = FALSE)
m3 <- lmer(Uric ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Uric ~ Julian + (1|Species), data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model with interaction is more informative

# transformation for time needed? visually no ----------------------------------
plot(birds$seconds_since_midnight, birds$Uric)

# interactions between season and capture time? no
m1 <- lmer(Uric ~ Event + seconds_since_midnight + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(Uric ~ Event * seconds_since_midnight + (1|Species), data = birds.cs, 
           REML = FALSE)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# interaction needed between SPEI and Ag?
m1 <- lmer(Uric ~ PercentAg + (1|Species), data = birds.cs, 
           REML = FALSE)
m2 <- lmer(Uric ~ SPEI + (1|Species), data = birds.cs, 
           REML = FALSE)
m3 <- lmer(Uric ~ PercentAg * SPEI + (1|Species), data = birds.cs, 
           REML = FALSE)

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model without interaction performs better

# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lmer(Uric ~ PercentAg + (1|Species), data = birds.cs, REML = FALSE)

# vegetation
m2 <- lmer(Uric ~ Percent_Total_Veg + (1|Species), data = birds.cs, REML = FALSE)

# habitat
m3 <- lmer(Uric ~ Permanence + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Uric ~ Percent_Exposed_Shoreline + (1|Species), 
           data = birds.cs, REML = FALSE)
m5 <- lmer(Uric ~ Dist_Closest_Wetland_m + (1|Species),  data = birds.cs, 
           REML = FALSE)

# weather
m6 <- lmer(Uric ~ SPEI + (1|Species), 
           data = birds.cs, REML = FALSE)

# life history
m7 <- lmer(Uric ~ Sex + (1|Species), 
           data = birds.cs, REML = FALSE)
m8 <- lmer(Uric ~ MigStatus + (1|Species), 
           data = birds.cs, REML = FALSE)

# temporal
m9 <- lmer(Uric ~ Julian + (1|Species), 
           data = birds.cs, REML = FALSE)
m10 <- lmer(Uric ~ Event + (1|Species), data = birds.cs,
            REML = FALSE)
m11 <- lmer(Uric ~ seconds_since_midnight + (1 | Species), 
            data = birds.cs, REML = FALSE)

# flock
m12 <- lmer(Uric ~ Max_Flock_Size + (1|Species), data = birds.cs, REML = FALSE)

# null
m13 <- lmer(Uric ~ 1 + (1 | Species), data = birds.cs, REML = FALSE)

# interactions
m14 <- lmer(Uric ~ Julian * MigStatus + (1 | Species), 
            data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:14)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m14)
confint(m14) # julian * migratory status very significant

summary(m11)
confint(m11) # capture time significant

summary(m8)
confint(m8) # migrants have higher levels of uric acid (makes sense)

summary(m4)
confint(m4) # % exposed shoreline informative; more shoreline, lower uric acid

summary(m2)
confint(m2) # more vegetation, higher uric acid levels

summary(m1)
confint(m1) # marginally significant? but in the opposite direction you would expect..

plot(birds$Percent_Exposed_Shoreline, birds$Uric)

# is random effect of site needed? ----
# exclude sites with only 1 bird
birds.m <- birds.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()


m1 <- lmer(Uric ~ PercentAg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# vegetation
m2 <- lmer(Uric ~ Percent_Total_Veg  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# habitat
m3 <- lmer(Uric ~ Permanence  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)
m4 <- lmer(Uric ~ Percent_Exposed_Shoreline  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m5 <- lmer(Uric ~ Dist_Closest_Wetland_m  + (1|Species) + (1|Site),  data = birds.m, 
           REML = FALSE)

# weather
m6 <- lmer(Uric ~ SPEI  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# life history
m7 <- lmer(Uric ~ Sex  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m8 <- lmer(Uric ~ MigStatus  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)

# temporal
m9 <- lmer(Uric ~ Julian  + (1|Species) + (1|Site), 
           data = birds.m, REML = FALSE)
m10 <- lmer(Uric ~ Event  + (1|Species) + (1|Site), data = birds.m,
            REML = FALSE)
m11 <- lmer(Uric ~ seconds_since_midnight + (1 | Species) + (1|Site), 
            data = birds.m, REML = FALSE)

# flock
m12 <- lmer(Uric ~ Max_Flock_Size  + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# null
m13 <- lmer(Uric ~ 1 + (1|Species) + (1|Site), data = birds.m, REML = FALSE)

# interactions
m14 <- lmer(Uric ~ Julian * MigStatus + (1 | Species) + (1|Site), 
            data = birds.m, REML = FALSE)

model_names <- paste0("m", 1:14)

models <- mget(model_names)

aictab(models, modnames = model_names) # results do not change


# Is there enough support to include the random effect? No ----
m1 <- lmer(Uric ~ seconds_since_midnight + (1 | Species), 
          data = birds.m, REML = FALSE)
m2 <- lmer(Uric ~ seconds_since_midnight + (1 | Species) + (1|Site), 
            data = birds.m, REML = FALSE)

m1 <- lmer(Uric ~ PercentAg  + (1|Species),  data = birds.m, 
           REML = FALSE)
m2 <- lmer(Uric ~ PercentAg  + (1|Species) + (1|Site),  data = birds.m, 
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



# Model Selection (Stage 2) ----------------------------------------------------
m1 <- lmer(Uric ~ Julian + (1|Species), data = birds.cs, REML = FALSE)
m2 <- lmer(Uric ~ MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m3 <- lmer(Uric ~ Julian * MigStatus + (1|Species), data = birds.cs, REML = FALSE)
m4 <- lmer(Uric ~ seconds_since_midnight + (1 | Species), data = birds.cs, 
          REML = FALSE)
m5 <- lmer(Uric ~ Percent_Exposed_Shoreline + (1|Species), 
            data = birds.cs, REML = FALSE)

m6 <- lmer(Uric ~ Percent_Exposed_Shoreline + seconds_since_midnight + (1|Species), 
            data = birds.cs, REML = FALSE)
m7 <- lmer(Uric ~ Percent_Exposed_Shoreline + Julian * MigStatus + (1|Species), 
            data = birds.cs, REML = FALSE)
m8 <- lmer(Uric ~ Percent_Exposed_Shoreline + Julian + (1|Species), 
            data = birds.cs, REML = FALSE)
m9 <- lmer(Uric ~ Percent_Exposed_Shoreline + MigStatus + (1|Species), 
            data = birds.cs, REML = FALSE)


m10 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus + (1|Species), 
            data = birds.cs, REML = FALSE)
m11 <- lmer(Uric ~ seconds_since_midnight + Julian + (1|Species), 
            data = birds.cs, REML = FALSE)
m12 <- lmer(Uric ~ seconds_since_midnight + MigStatus + (1|Species), 
            data = birds.cs, REML = FALSE)

m13 <- lmer(Uric ~ Percent_Exposed_Shoreline + seconds_since_midnight + 
             Julian * MigStatus + (1|Species), 
           data = birds.cs, REML = FALSE)

m14 <- lmer(Uric ~ Percent_Exposed_Shoreline + Julian + seconds_since_midnight +
             (1|Species), 
           data = birds.cs, REML = FALSE)

m15 <- lmer(Uric ~ Percent_Exposed_Shoreline + MigStatus + 
             seconds_since_midnight + (1|Species), 
           data = birds.cs, REML = FALSE)

model_names <- paste0("m", 1:15)

models <- mget(model_names)

aictab(models, modnames = model_names)

# plot exposed shoreline -------------------------------------------------------
m <- lmer(Uric ~ Percent_Exposed_Shoreline + seconds_since_midnight + 
            Julian * MigStatus + (1 | Species), 
          data = birds, REML = FALSE)

d <- expand.grid(Percent_Exposed_Shoreline = seq(min(birds$Percent_Exposed_Shoreline), 
                                                 max(birds$Percent_Exposed_Shoreline), 
                                                 length = 1000),
                 Julian = mean(birds$Julian),
                 MigStatus = unique(birds$MigStatus),
                 seconds_since_midnight = mean(birds$seconds_since_midnight),
                 Species = unique(birds$Species)) 

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = Percent_Exposed_Shoreline, y = fit)) +
  geom_line(aes(color = MigStatus), size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = MigStatus), 
              alpha = 0.25, color = NA, show.legend = TRUE) +
  theme_classic() +
  labs(x = "% Exposed Shoreline", 
       y = expression("Uric Acid Levels (" *mu*"mol/L)"),
       color = "Migratory Status",
       fill = "Migratory Status") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = birds, 
             aes(x = Percent_Exposed_Shoreline, y = Uric, color = MigStatus), size = 2.5) +
  scale_color_viridis_d(alpha = 1, begin = 0, end = 0.8) +
  scale_fill_viridis_d(alpha = 1, begin = 0, end = 0.8) +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 2250, by = 500)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(0, 100, by = 10)) +
  coord_cartesian(ylim = c(-15,2250),
                  xlim = c(5,105))




summary(m)
confint(m)

# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- lmer(Uric ~ Biomass + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lmer(Uric ~ Diversity + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of diversity


# do neonics explain any variation in uric acid? ----
# informative covariates: time, Julian * MigStatus, exposed shoreline

# summary statistics----
table(birds$PlasmaDetection) # n: 57, y: 28 (n = 85)
table(birds$WaterNeonicDetection) # n: 81, y: 4 (n = 85)
table(birds$AnyDetection) # n: 22, y: 63 (n = 85)
table(birds$WaterOrInvertDetection) # n: 37, y: 48 (n = 85)
table(birds$InvertPesticideDetection) # n: 23, y: 44 (n = 67)

mean(birds$OverallNeonic, na.rm = TRUE) # 14.2 ug/L
sd(birds$OverallNeonic, na.rm = TRUE) # 111 ug/L


# water neonic detection --> neonics not informative
birds.clean.water <- birds.cs[!is.na(birds.cs$WaterNeonicDetection), ] #n = 85

m1 <- lmer(Uric ~ seconds_since_midnight + (1|Species), 
           data = birds.clean.water,
           REML = FALSE)

m2 <- lmer(Uric ~ Julian * MigStatus +
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m3 <- lmer(Uric ~ Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m4 <- lmer(Uric ~ WaterNeonicDetection + 
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m5 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m6 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m7 <- lmer(Uric ~ seconds_since_midnight + WaterNeonicDetection +
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m8 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m9 <- lmer(Uric ~ Julian * MigStatus + WaterNeonicDetection + 
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m10 <- lmer(Uric ~ Percent_Exposed_Shoreline + WaterNeonicDetection + 
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m11 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
             Percent_Exposed_Shoreline +
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m12 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
             WaterNeonicDetection +
             (1|Species), data = birds.clean.water,
           REML = FALSE)

m13 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
             WaterNeonicDetection +
           (1|Species), data = birds.clean.water,
           REML = FALSE)


m14 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
             WaterNeonicDetection + seconds_since_midnight +
           (1|Species), data = birds.clean.water,
           REML = FALSE)



# invertebrate pesticide detection --> invert detections correspond with lower uric acid levels...
birds.clean.invert <- birds.cs[!is.na(birds.cs$InvertPesticideDetection), ] #n = 67

m1 <- lmer(Uric ~ seconds_since_midnight + (1|Species), 
           data = birds.clean.invert,
           REML = FALSE)

m2 <- lmer(Uric ~ Julian * MigStatus +
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m3 <- lmer(Uric ~ Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m4 <- lmer(Uric ~ InvertPesticideDetection + 
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m5 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m6 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m7 <- lmer(Uric ~ seconds_since_midnight + InvertPesticideDetection +
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m8 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m9 <- lmer(Uric ~ Julian * MigStatus + InvertPesticideDetection + 
             (1|Species), data = birds.clean.invert,
           REML = FALSE)

m10 <- lmer(Uric ~ Percent_Exposed_Shoreline + InvertPesticideDetection + 
              (1|Species), data = birds.clean.invert,
            REML = FALSE)

m11 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
              Percent_Exposed_Shoreline +
              (1|Species), data = birds.clean.invert,
            REML = FALSE)

m12 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
              InvertPesticideDetection +
              (1|Species), data = birds.clean.invert,
            REML = FALSE)

m13 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              InvertPesticideDetection +
              (1|Species), data = birds.clean.invert,
            REML = FALSE)


m14 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              InvertPesticideDetection + seconds_since_midnight +
              (1|Species), data = birds.clean.invert,
            REML = FALSE)

# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
birds.clean.waterorinvert <- birds.cs[!is.na(birds.cs$WaterOrInvertDetection), ] #n = 85

m1 <- lmer(Uric ~ seconds_since_midnight + (1|Species), 
           data = birds.clean.waterorinvert,
           REML = FALSE)

m2 <- lmer(Uric ~ Julian * MigStatus +
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m3 <- lmer(Uric ~ Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m4 <- lmer(Uric ~ WaterOrInvertDetection + 
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m5 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m6 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m7 <- lmer(Uric ~ seconds_since_midnight + WaterOrInvertDetection +
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m8 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m9 <- lmer(Uric ~ Julian * MigStatus + WaterOrInvertDetection + 
             (1|Species), data = birds.clean.waterorinvert,
           REML = FALSE)

m10 <- lmer(Uric ~ Percent_Exposed_Shoreline + WaterOrInvertDetection + 
              (1|Species), data = birds.clean.waterorinvert,
            REML = FALSE)

m11 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
              Percent_Exposed_Shoreline +
              (1|Species), data = birds.clean.waterorinvert,
            REML = FALSE)

m12 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
              WaterOrInvertDetection +
              (1|Species), data = birds.clean.waterorinvert,
            REML = FALSE)

m13 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              WaterOrInvertDetection +
              (1|Species), data = birds.clean.waterorinvert,
            REML = FALSE)


m14 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              WaterOrInvertDetection + seconds_since_midnight +
              (1|Species), data = birds.clean.waterorinvert,
            REML = FALSE)
# shorebird plasma detection --> neonics not informative
birds.clean.plasma <- birds.cs[!is.na(birds.cs$PlasmaDetection), ] #n = 85


m1 <- lmer(Uric ~ seconds_since_midnight + (1|Species), 
           data = birds.clean.plasma,
           REML = FALSE)

m2 <- lmer(Uric ~ Julian * MigStatus +
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m3 <- lmer(Uric ~ Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m4 <- lmer(Uric ~ PlasmaDetection + 
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m5 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m6 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m7 <- lmer(Uric ~ seconds_since_midnight + PlasmaDetection +
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m8 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m9 <- lmer(Uric ~ Julian * MigStatus + PlasmaDetection + 
             (1|Species), data = birds.clean.plasma,
           REML = FALSE)

m10 <- lmer(Uric ~ Percent_Exposed_Shoreline + PlasmaDetection + 
              (1|Species), data = birds.clean.plasma,
            REML = FALSE)

m11 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
              Percent_Exposed_Shoreline +
              (1|Species), data = birds.clean.plasma,
            REML = FALSE)

m12 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
              PlasmaDetection +
              (1|Species), data = birds.clean.plasma,
            REML = FALSE)

m13 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              PlasmaDetection +
              (1|Species), data = birds.clean.plasma,
            REML = FALSE)


m14 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              PlasmaDetection + seconds_since_midnight +
              (1|Species), data = birds.clean.plasma,
            REML = FALSE)

# any detection (plasma or environmental) --> neonics not informative
birds.clean.any <- birds.cs[!is.na(birds.cs$AnyDetection), ] #n = 85


m1 <- lmer(Uric ~ seconds_since_midnight + (1|Species), 
           data = birds.clean.any,
           REML = FALSE)

m2 <- lmer(Uric ~ Julian * MigStatus +
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m3 <- lmer(Uric ~ Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m4 <- lmer(Uric ~ AnyDetection + 
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m5 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m6 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m7 <- lmer(Uric ~ seconds_since_midnight + AnyDetection +
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m8 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline + 
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m9 <- lmer(Uric ~ Julian * MigStatus + AnyDetection + 
             (1|Species), data = birds.clean.any,
           REML = FALSE)

m10 <- lmer(Uric ~ Percent_Exposed_Shoreline + AnyDetection + 
              (1|Species), data = birds.clean.any,
            REML = FALSE)

m11 <- lmer(Uric ~ seconds_since_midnight + Julian * MigStatus +
              Percent_Exposed_Shoreline +
              (1|Species), data = birds.clean.any,
            REML = FALSE)

m12 <- lmer(Uric ~ seconds_since_midnight + Percent_Exposed_Shoreline +
              AnyDetection +
              (1|Species), data = birds.clean.any,
            REML = FALSE)

m13 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              AnyDetection +
              (1|Species), data = birds.clean.any,
            REML = FALSE)


m14 <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              AnyDetection + seconds_since_midnight +
              (1|Species), data = birds.clean.any,
            REML = FALSE)

### ...AIC 
models <- list(m1, m2, m3,m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
model.sel(models)

# model summaries:
summary(m14)
confint(m14)


plot(m14)

# remove random effect for invert detection and uric acid result ----
m <- lmer(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
              InvertPesticideDetection + seconds_since_midnight +
              (1|Species), data = birds,
            REML = FALSE)

m <- lm(Uric ~ Julian * MigStatus + Percent_Exposed_Shoreline +
            InvertPesticideDetection + seconds_since_midnight, data = birds.cs)
