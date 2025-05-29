#----------------------------------------------#
# Lesser Yellowlegs Uric Acid Habitat Analysis #
#             linear regression                #
#             Created 2025-04-10               #
#             Modified 2025-05-29             #
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

# categorize factors
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

# subset birds that have uric acid levels
leye <- leye %>% 
  filter(!is.na(Uric))

# ...standardize data except for response ----
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("Uric"), scale))

# Test for Correlations --------------------------------------------------------
# subset data
sample <- leye[, c("PercentAg",
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
                   "Age",
                   "DominantCrop",
                   "NearestCropDistance_m",
                   "Dist_Closest_Wetland_m",
                   "Max_Flock_Size",
                   "InvertPesticideDetection",
                   "WaterNeonicDetection",
                   "WaterOrInvertDetection",
                   "AnyDetection"
)]


# convert categorical to numeric for correlation matrix
sample$Event <- as.numeric(sample$Event)
sample$Permanence <- as.numeric(sample$Permanence)
sample$AgCategory <- as.numeric(sample$AgCategory)
sample$DominantCrop <- as.numeric(sample$DominantCrop)
sample$Sex <- as.numeric(sample$Sex)
sample$PlasmaDetection <- as.numeric(sample$PlasmaDetection)
sample$Site <- as.numeric(sample$Site)
sample$Age <- as.numeric(sample$Age)
sample$InvertPesticideDetection <- as.numeric(sample$InvertPesticideDetection)
sample$WaterNeonicDetection <- as.numeric(sample$WaterNeonicDetection)
sample$WaterOrInvertDetection <- as.numeric(sample$WaterOrInvertDetection)
sample$AnyDetection <- as.numeric(sample$AnyDetection)

cor(sample)
# correlations > 0.6:

cor
# PercentVeg + Dominant Crop (-0.71)
# PercentAg & AgCategory (0.80)
# PercentAg + Dominant Crop (0.62)
# NearestCropDistance + PercentAg (-0.656)
# permanence and exposed shoreline (-0.62)
# Ag Category and permanence (-0.94)
# SPEI & %exposed shoreline (-0.63)
# Ag category & nearest crop distance (-0.73)
# SPEI and permanence (-0.94)
# Max flock size and dominant crop (0.69)
# Water or Invert Detection & permanence (-0.7487)
# Any detection & permanence (-0.64)
# Water neonic detection & Event
# Julian & water neonic detection


# not really appropriate to look at drought or Event (small sample size)

# is cyclical time transformation needed? no
m1 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Uric ~ seconds_since_midnight + sin + cos, data = leye.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# ag and SPEI?
m1 <- lm(Uric ~ SPEI, data = leye.cs)
m2 <- lm(Uric ~ PercentAg, data = leye.cs)
m3 <- lm(Uric ~ SPEI * PercentAg, data = leye.cs)
m4 <- lm(Uric ~ SPEI + PercentAg, data = leye.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)


# Stage 1 Model Selection ----
m1 <- lm(Uric ~ PercentAg + seconds_since_midnight, data = leye.cs)

# vegetation
m2 <- lm(Uric ~ Percent_Total_Veg + seconds_since_midnight, data = leye.cs)

# habitat
m3 <- lm(Uric ~ Permanence + seconds_since_midnight, data = leye.cs)
m4 <- lm(Uric ~ Percent_Exposed_Shoreline + seconds_since_midnight, data = leye.cs)
m5 <- lm(Uric ~ Dist_Closest_Wetland_m + seconds_since_midnight, data = leye.cs)

# weather
m6 <- lm(Uric ~ SPEI + seconds_since_midnight, data = leye.cs)

# life history
m7 <- lm(Uric ~ Age + seconds_since_midnight, data = leye.cs)
m8 <- lm(Uric ~ Sex + seconds_since_midnight, data = leye.cs)

# temporal
m9 <- lm(Uric ~ Julian + seconds_since_midnight, data = leye.cs)
m10 <- lm(Uric ~ Event + seconds_since_midnight, data = leye.cs)
m11 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)

# flock
m12 <- lm(Uric ~ Max_Flock_Size + seconds_since_midnight, data = leye.cs)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m6)
confint(m6)

# plot top model for thesis ----
m <- lm(Uric ~ SPEI + seconds_since_midnight,
        data = leye)

d <- expand.grid(SPEI = seq(min(leye$SPEI), 
                            max(leye$SPEI), 
                            length.out = 1000),
                 seconds_since_midnight = mean(leye$seconds_since_midnight))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = SPEI, y = yhat)) +
  geom_line(size = 1, col = "goldenrod3") +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,  fill = "goldenrod3") +
  labs(x = "Drought Index (SPEI)", 
       y = expression("Lesser Yellowlegs Uric Acid Levels ("*mu*"mol/L)")) +
  theme(axis.title.x = element_text(size = 20,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = SPEI, y = Uric), size = 2.5,
             col = "goldenrod3") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, 2250, by = 500)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(-1.5, 1.5, by = 0.5)) +
  coord_cartesian(ylim = c(-15,2250),
                  xlim = c(-1.75,1.1)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1, color = "gray")



# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(leye.cs, !is.na(Biomass))

m <- lm(Uric ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lm(Uric ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of diversity


# do neonics explain any further variation of uric  than event * time? ----
# informative covariates: event * time

# summary statistics----
table(leye$PlasmaDetection) # n: 21, y: 8 (n = 29)
table(leye$WaterNeonicDetection) # n: 25, y: 4 (n = 29)
table(leye$AnyDetection) # n: 10, y: 19 (n = 29)
table(leye$WaterOrInvertDetection) # n: 11, y: 18 (n = 29)
table(leye$InvertPesticideDetection) # n: 10, y: 14 (n = 24)

mean(leye$OverallNeonic, na.rm = TRUE) # 0.728 ug/L
sd(leye$OverallNeonic, na.rm = TRUE) # 1.96 ug/L

# water neonic detection --> neonics not informative
m1 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Uric ~ seconds_since_midnight + SPEI, data = leye.cs)
m3 <- lm(Uric ~ seconds_since_midnight + WaterNeonicDetection, data = leye.cs)
m4 <- lm(Uric ~ seconds_since_midnight + SPEI +
           WaterNeonicDetection, data = leye.cs)

# invertebrate pesticide detection --> neonics not informative
leye.clean.invert <- leye.cs[!is.na(leye.cs$InvertPesticideDetection), ] #n = 24

m1 <- lm(Uric ~ seconds_since_midnight, data = leye.clean.invert)
m2 <- lm(Uric ~ seconds_since_midnight + SPEI, data = leye.clean.invert)
m3 <- lm(Uric ~ seconds_since_midnight + InvertPesticideDetection, data = leye.clean.invert)
m4 <- lm(Uric ~ seconds_since_midnight + SPEI +
           InvertPesticideDetection, data = leye.clean.invert)


# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
m1 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Uric ~ seconds_since_midnight + SPEI, data = leye.cs)
m3 <- lm(Uric ~ seconds_since_midnight + WaterOrInvertDetection, data = leye.cs)
m4 <- lm(Uric ~ seconds_since_midnight + SPEI +
           WaterOrInvertDetection, data = leye.cs)


# shorebird plasma detection --> neonics not informative
m1 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Uric ~ seconds_since_midnight + SPEI, data = leye.cs)
m3 <- lm(Uric ~ seconds_since_midnight + PlasmaDetection, data = leye.cs)
m4 <- lm(Uric ~ seconds_since_midnight + SPEI +
           PlasmaDetection, data = leye.cs)

# any detection (plasma or environmental) --> neonics not informative
m1 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Uric ~ seconds_since_midnight + SPEI, data = leye.cs)
m3 <- lm(Uric ~ seconds_since_midnight + AnyDetection, data = leye.cs)
m4 <- lm(Uric ~ seconds_since_midnight + SPEI +
           AnyDetection, data = leye.cs)

### ...AIC 
models <- list(m1, m2, m3, m4)
model.sel(models)

# model summaries:
summary(m4)
confint(m4)













