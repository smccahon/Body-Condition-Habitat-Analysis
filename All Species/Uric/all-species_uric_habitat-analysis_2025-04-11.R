#----------------------------------------------#
#    All Species Uric Acid Habitat Analysis    #
#            Created 2025-04-11                #
#           Modified 2025-04-11                #
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
  mutate(across(where(is.numeric), scale))

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
# % ag and dominant crop
# % veg and % shoreline
# event and julian
# SPEI and julian

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

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m7)
confint(m10)

# plot exposed shoreline -------------------------------------------------------
m <- lmer(Uric ~ Percent_Exposed_Shoreline + Julian * MigStatus + (1 | Species), 
          data = birds)

d <- expand.grid(Percent_Exposed_Shoreline = seq(min(birds$Percent_Exposed_Shoreline), 
                                                 max(birds$Percent_Exposed_Shoreline), 
                                                 length = 1000),
                 Julian = mean(birds$Julian),
                 MigStatus = unique(birds$MigStatus),
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
       y = expression("All Species Uric Acid Levels (" *mu*"mol/L)"),
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
             aes(x = Percent_Exposed_Shoreline, y = Uric, color = MigStatus), size = 2) +
  scale_color_viridis_d(alpha = 1, begin = 0, end = 0.8) +
  scale_fill_viridis_d(alpha = 1, begin = 0, end = 0.8)







# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(birds.cs, !is.na(Biomass))

m <- lmer(Uric ~ Biomass + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lmer(Uric ~ Diversity + (1 | Species), data = birds.sub)

summary(m)
confint(m) # no effect of diversity








