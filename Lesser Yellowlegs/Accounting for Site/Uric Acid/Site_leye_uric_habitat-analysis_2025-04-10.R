#----------------------------------------------#
# Lesser Yellowlegs Uric Acid Habitat Analysis #
#             linear regression                #
#             Created 2025-04-10               #
#             Modified 2025-05-05             #
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

# not really appropriate to look at drought or Event (small sample size)

# is cyclical time transformation needed? no
m1 <- lm(Uric ~ seconds_since_midnight, data = leye.cs)
m2 <- lm(Uric ~ seconds_since_midnight + sin + cos, data = leye.cs)

model_names <- paste0("m", 1:2)

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
       y = expression("Lesser Yellowlegs Uric Acid Levels ( "*mu*"mol/L)")) +
  theme(axis.title.x = element_text(size = 20,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 20,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_point(data = leye, aes(x = SPEI, y = Uric), size = 2,
             col = "goldenrod3")




# are random effects necessary? ----
leye.m <- leye.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()

m1 <- lmer(Uric ~ PercentAg + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)

# vegetation
m2 <- lmer(Uric ~ Percent_Total_Veg + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)

# habitat
m3 <- lmer(Uric ~ Permanence + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)
m4 <- lmer(Uric ~ Percent_Exposed_Shoreline + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)
m5 <- lmer(Uric ~ Dist_Closest_Wetland_m + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)

# weather
m6 <- lmer(Uric ~ SPEI + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)

# life history
m7 <- lmer(Uric ~ Age + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)
m8 <- lmer(Uric ~ Sex + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)

# temporal
m9 <- lmer(Uric ~ Julian + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)
m10 <- lmer(Uric ~ Event + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)
m11 <- lmer(Uric ~ seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)

# flock
m12 <- lmer(Uric ~ Max_Flock_Size + seconds_since_midnight+(1|Site), REML = FALSE, data = leye.m)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names) # results changed

# run reduced dataset without site ----
# Stage 1 Model Selection ----
m1 <- lm(Uric ~ PercentAg + seconds_since_midnight, data = leye.m)

# vegetation
m2 <- lm(Uric ~ Percent_Total_Veg + seconds_since_midnight, data = leye.m)

# habitat
m3 <- lm(Uric ~ Permanence + seconds_since_midnight, data = leye.m)
m4 <- lm(Uric ~ Percent_Exposed_Shoreline + seconds_since_midnight, data = leye.m)
m5 <- lm(Uric ~ Dist_Closest_Wetland_m + seconds_since_midnight, data = leye.m)

# weather
m6 <- lm(Uric ~ SPEI + seconds_since_midnight, data = leye.m)

# life history
m7 <- lm(Uric ~ Age + seconds_since_midnight, data = leye.m)
m8 <- lm(Uric ~ Sex + seconds_since_midnight, data = leye.m)

# temporal
m9 <- lm(Uric ~ Julian + seconds_since_midnight, data = leye.m)
m10 <- lm(Uric ~ Event + seconds_since_midnight, data = leye.m)
m11 <- lm(Uric ~ seconds_since_midnight, data = leye.m)

# flock
m12 <- lm(Uric ~ Max_Flock_Size + seconds_since_midnight, data = leye.m)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)




# Is there enough support to include the random effect? No ----
m1 <- lm(Uric ~ SPEI + seconds_since_midnight, data = leye.m)
m2 <- lmer(Uric ~ SPEI + seconds_since_midnight + (1|Site),
           data = leye.m, REML = FALSE)


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
birds.sub <- subset(leye.cs, !is.na(Biomass))

m <- lm(Uric ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lm(Uric ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of diversity



# assess heteroscedasticity - none; good model
m <- lm(Uric ~ SPEI + seconds_since_midnight, data = leye.cs)
plot(predict(m), rstudent(m))












