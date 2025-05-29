#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#               linear regression                    #
#               Created 2025-04-24                   #
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


# filter birds that only contain triglycerides (n = 31)
leye <- leye %>% 
  filter(!is.na(Tri))

# standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("Tri"), scale))

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
# nearest crop distance and ag category
# % veg and dominant crop
# permanence & % exposed shoreline
# permanence and SPEI
# SPEI and % exposed shoreline
# julian and detection (0.60 exactly)
# dominant crop and max flock size

# is cyclical time transformation needed? no
m1 <- lm(Tri ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Tri ~ seconds_since_midnight + sin + cos, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# spei and ag?
m1 <- lm(Beta ~ PercentAg + seconds_since_midnight, data = leye.cs)
m2 <- lm(Beta ~ SPEI + seconds_since_midnight, data = leye.cs)
m3 <- lm(Beta ~ PercentAg + SPEI + seconds_since_midnight, data = leye.cs)
m4 <- lm(Beta ~ PercentAg * SPEI + seconds_since_midnight, data = leye.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Stage 1 Model Selection ----
m1 <- lm(Tri ~ PercentAg + seconds_since_midnight, data = leye.cs)

# vegetation
m2 <- lm(Tri ~ Percent_Total_Veg + seconds_since_midnight, data = leye.cs)

# habitat
m3 <- lm(Tri ~ Permanence + seconds_since_midnight, data = leye.cs)
m4 <- lm(Tri ~ Percent_Exposed_Shoreline + seconds_since_midnight, data = leye.cs)
m5 <- lm(Tri ~ Dist_Closest_Wetland_m + seconds_since_midnight, data = leye.cs)

# weather
m6 <- lm(Tri ~ SPEI + seconds_since_midnight, data = leye.cs)

# life history
m7 <- lm(Tri ~ Age + seconds_since_midnight, data = leye.cs)
m8 <- lm(Tri ~ Sex + seconds_since_midnight, data = leye.cs)

# temporal
m9 <- lm(Tri ~ Julian + seconds_since_midnight, data = leye.cs)
m10 <- lm(Tri ~ Event + seconds_since_midnight, data = leye.cs)
m11 <- lm(Tri ~ seconds_since_midnight, data = leye.cs)

# flock
m12 <- lm(Tri ~ Max_Flock_Size + seconds_since_midnight, data = leye.cs)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)
summary(m1) # percent ag not significant

plot(leye$PercentAg, leye$Tri)



# BETA ----

# filter birds that only contain beta (n = 29)
leye <- leye %>% 
  filter(!is.na(Beta))

# standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

# is cyclical time transformation needed? YES
m1 <- lm(Beta ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Stage 1 Model Selection ----
m1 <- lm(Beta ~ PercentAg + seconds_since_midnight + sin + cos, data = leye.cs)

# vegetation
m2 <- lm(Beta ~ Percent_Total_Veg + seconds_since_midnight + sin + cos, data = leye.cs)

# habitat
m3 <- lm(Beta ~ Permanence + seconds_since_midnight + sin + cos, data = leye.cs)
m4 <- lm(Beta ~ Percent_Exposed_Shoreline + seconds_since_midnight + sin + cos, data = leye.cs)
m5 <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + sin + cos, data = leye.cs)

# weather
m6 <- lm(Beta ~ SPEI + seconds_since_midnight + sin + cos, data = leye.cs)

# life history
m7 <- lm(Beta ~ Age + seconds_since_midnight + sin + cos, data = leye.cs)
m8 <- lm(Beta ~ Sex + seconds_since_midnight + sin + cos, data = leye.cs)

# temporal
m9 <- lm(Beta ~ Julian + seconds_since_midnight + sin + cos, data = leye.cs)
m10 <- lm(Beta ~ Event + seconds_since_midnight + sin + cos, data = leye.cs)
m11 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.cs)

# flock
m12 <- lm(Beta ~ Max_Flock_Size + seconds_since_midnight + sin + cos, data = leye.cs)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# informative parameters: dist to wetland & time

confint(m5)
summary(m5) # lose fat the closer you are from another wetland

summary(m2)
confint(m2) # no effect of vegetation

plot(leye$Dist_Closest_Wetland_m, leye$Beta)

m <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye)

d <- expand.grid(Dist_Closest_Wetland_m = seq(min(leye$Dist_Closest_Wetland_m), 
                                 max(leye$Dist_Closest_Wetland_m), 
                                 length.out = 1000),
                 seconds_since_midnight = mean(leye$seconds_since_midnight))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = Dist_Closest_Wetland_m, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Dist. to Nearest Wetland (m)", 
       y = "LEYE Fat Metabolism (BHOB mmol/L)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = Dist_Closest_Wetland_m, y = Beta), size = 3)


cor(leye$Tri, leye$Beta)


# assess heteroscedasticity
plot(predict(m), rstudent(m)) # very mild heteroscedasticity

# weighted least squares --> didn't change anything...
leye %>% summarize(variance = var(Beta), weight = 1/var(Beta))

# put results back into data frame
leye <- leye %>% mutate(w = 1/var(Beta))

m.wls <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye,
        weights = w)

# no difference for some reason
plot(predict(m.wls),rstudent(m.wls))
plot(predict(m),rstudent(m))


# iteratively weighted least squares
m.ols <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + 
              sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
              cos(2 * pi * seconds_since_midnight / (24 * 3600)),
            data = leye)

leye$w <- 1

for (i in 1:10) {
  m.wls <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + 
                sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
                cos(2 * pi * seconds_since_midnight / (24 * 3600)),
              data = leye,
              weights = w)
  print(coef(m.wls))
  leye$w <- 1 / predict(m.wls)
}

plot(predict(m.wls), rstudent(m.wls))
confint(m.wls)
summary(m.wls)

# plot dist to wetland & time with WLS ----
d <- expand.grid(Dist_Closest_Wetland_m = seq(min(leye$Dist_Closest_Wetland_m), 
                                              max(leye$Dist_Closest_Wetland_m), 
                                              length.out = 1000),
                 seconds_since_midnight = mean(leye$seconds_since_midnight))

predictions <- predict(m.wls, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m.wls, newdata = d)

d$se <- predict(m.wls, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = Dist_Closest_Wetland_m, y = yhat)) +
  geom_line(size = 1, col = "goldenrod3") +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "goldenrod3", alpha = 0.2) +
  labs(x = "Distance to Nearest Wetland (m)", 
       y = "Lesser Yellowlegs BHB Levels (mmol/L)") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = Dist_Closest_Wetland_m, y = Beta), size = 2.5,
             col = "goldenrod3") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 1.60, by = 0.25)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(0, 375, by = 50)) +
  coord_cartesian(ylim = c(0,1.6),
                  xlim = c(-0.2, 375))

  # geom_text(data = leye, aes(x = Dist_Closest_Wetland_m, y = Beta, 
  #                            label = seconds_since_midnight),
  #           size = 3.5, hjust = -0.1, vjust = 0.5, color = "black", inherit.aes = FALSE)
  

# do neonics explain any further variation of triglycerides than time (cyclical)? ----
# informative covariates: time (cyclical)

# summary statistics----
table(leye$PlasmaDetection) # n: 23, y: 8 (n = 31)
table(leye$WaterNeonicDetection) # n: 27, y: 4 (n = 31)
table(leye$AnyDetection) # n: 11, y: 20 (n = 31)
table(leye$WaterOrInvertDetection) # n: 12, y: 19 (n = 31)
table(leye$InvertPesticideDetection) # n: 11, y: 15 (n = 26)

mean(leye$OverallNeonic, na.rm = TRUE) # 0.681 ug/L
sd(leye$OverallNeonic, na.rm = TRUE) # 1.9 ug/L

# water neonic detection --> neonics not informative
leye.clean.water <- leye.cs[!is.na(leye.cs$WaterNeonicDetection), ] #n = 31

m1 <- lm(Tri ~ seconds_since_midnight + sin + cos, data = leye.clean.water)
m2 <- lm(Tri ~ seconds_since_midnight + sin + cos + WaterNeonicDetection, data = leye.clean.water)

# invertebrate pesticide detection --> neonics not informative
leye.clean.invert <- leye.cs[!is.na(leye.cs$InvertPesticideDetection), ] #n = 26

m1 <- lm(Tri ~ seconds_since_midnight + sin + cos, data = leye.clean.invert)
m2 <- lm(Tri ~ seconds_since_midnight + sin + cos + InvertPesticideDetection, data = leye.clean.invert)

# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
leye.clean.waterorinvert <- leye.cs[!is.na(leye.cs$WaterOrInvertDetection), ] #n = 31
m1 <- lm(Tri ~ seconds_since_midnight + sin + cos, data = leye.clean.waterorinvert)
m2 <- lm(Tri ~ seconds_since_midnight + sin + cos + 
           WaterOrInvertDetection, data = leye.clean.waterorinvert)

# shorebird plasma detection --> neonics not informative
m1 <- lm(Tri ~ seconds_since_midnight + sin + cos, data = leye.cs)
m2 <- lm(Tri ~ seconds_since_midnight + sin + cos + PlasmaDetection, data = leye.cs)

# any detection (plasma or environmental) --> neonics not informative
m1 <- lm(Tri ~ seconds_since_midnight + sin + cos, data = leye.cs)
m2 <- lm(Tri ~ seconds_since_midnight + sin + cos + AnyDetection, data = leye.cs)

### ...AIC 
models <- list(m1, m2)
model.sel(models)

# model summaries:
summary(m2)
confint(m2)



# do neonics explain any further variation of beta than time (cyclical)? ----
# informative covariates: time (cyclical), dist. to wetland

# filter birds that only contain beta (n = 29)
leye <- leye %>% 
  filter(!is.na(Beta))

# standardize data
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("Beta"), scale))

# summary statistics----
table(leye$PlasmaDetection) # n: 21, y: 8 (n = 29)
table(leye$WaterNeonicDetection) # n: 25, y: 4 (n = 29)
table(leye$AnyDetection) # n: 10, y: 19 (n = 29)
table(leye$WaterOrInvertDetection) # n: 11, y: 18 (n = 29)
table(leye$InvertPesticideDetection) # n: 10, y: 14 (n = 24)

mean(leye$OverallNeonic, na.rm = TRUE) # 0.728 ug/L
sd(leye$OverallNeonic, na.rm = TRUE) # 1.96 ug/L

# water neonic detection --> neonics not informative
leye.clean.water <- leye.cs[!is.na(leye.cs$WaterNeonicDetection), ] #n = 29

m1 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.clean.water)
m2 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m, data = leye.clean.water)
m3 <- lm(Beta ~ seconds_since_midnight + sin + cos + WaterNeonicDetection, data = leye.clean.water)
m4 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m +
           WaterNeonicDetection, data = leye.clean.water)

# invertebrate pesticide detection --> neonics not informative
leye.clean.invert <- leye.cs[!is.na(leye.cs$InvertPesticideDetection), ] #n = 24

m1 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.clean.invert)
m2 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m, data = leye.clean.invert)
m3 <- lm(Beta ~ seconds_since_midnight + sin + cos + InvertPesticideDetection, data = leye.clean.invert)
m4 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m +
           InvertPesticideDetection, data = leye.clean.invert)

# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
leye.clean.waterorinvert <- leye.cs[!is.na(leye.cs$WaterOrInvertDetection), ] #n = 29

m1 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.clean.waterorinvert)
m2 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m, data = leye.clean.waterorinvert)
m3 <- lm(Beta ~ seconds_since_midnight + sin + cos + WaterOrInvertDetection, data = leye.clean.waterorinvert)
m4 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m +
           WaterOrInvertDetection, data = leye.clean.waterorinvert)


# shorebird plasma detection --> neonics not informative #n = 29
m1 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.cs)
m2 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m, data = leye.cs)
m3 <- lm(Beta ~ seconds_since_midnight + sin + cos + PlasmaDetection, data = leye.cs)
m4 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m +
           PlasmaDetection, data = leye.cs)

# any detection (plasma or environmental) --> neonics not informative #n = 29
m1 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.cs)
m2 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m, data = leye.cs)
m3 <- lm(Beta ~ seconds_since_midnight + sin + cos + AnyDetection, data = leye.cs)
m4 <- lm(Beta ~ seconds_since_midnight + sin + cos + Dist_Closest_Wetland_m +
           AnyDetection, data = leye.cs)

### ...AIC 
models <- list(m1, m2, m3, m4)
model.sel(models)

# model summaries:
summary(m4)
confint(m4)



