#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#               linear regression                    #
#               Created 2025-04-24                   #
#              Modified 2025-05-05                   #
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


# filter birds that only contain Tri (n = 29)
leye <- leye %>% 
  filter(!is.na(Tri))

# ...standardize data except for response ----
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

# ...standardize data except for response ----
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("Beta"), scale))


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
summary(m5) # burn more fat the closer you are to another wetland

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
       y = "LEYE Fat Metabolism (BHB mmol/L)") +
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


# are random effects necessary? ----
leye.m <- leye.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()

m1 <- lmer(Beta ~ PercentAg + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)

# vegetation
m2 <- lmer(Beta ~ Percent_Total_Veg + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)

# habitat
m3 <- lmer(Beta ~ Permanence + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)
m4 <- lmer(Beta ~ Percent_Exposed_Shoreline + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)
m5 <- lmer(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)

# weather
m6 <- lmer(Beta ~ SPEI + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)

# life history
m7 <- lmer(Beta ~ Age + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)
m8 <- lmer(Beta ~ Sex + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)

# temporal
m9 <- lmer(Beta ~ Julian + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)
m10 <- lmer(Beta ~ Event + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)
m11 <- lmer(Beta ~ seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)

# flock
m12 <- lmer(Beta ~ Max_Flock_Size + seconds_since_midnight + sin + cos + (1|Site), REML = FALSE, data = leye.m)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names) # dist. to closest wetland still informative after incorporating site

# run reduced dataset with site as a random effect ----
# Stage 1 Model Selection ----
m1 <- lm(Beta ~ PercentAg + seconds_since_midnight + sin + cos, data = leye.m)

# vegetation
m2 <- lm(Beta ~ Percent_Total_Veg + seconds_since_midnight + sin + cos, data = leye.m)

# habitat
m3 <- lm(Beta ~ Permanence + seconds_since_midnight + sin + cos, data = leye.m)
m4 <- lm(Beta ~ Percent_Exposed_Shoreline + seconds_since_midnight + sin + cos, data = leye.m)
m5 <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + sin + cos, data = leye.m)

# weather
m6 <- lm(Beta ~ SPEI + seconds_since_midnight + sin + cos, data = leye.m)

# life history
m7 <- lm(Beta ~ Age + seconds_since_midnight + sin + cos, data = leye.m)
m8 <- lm(Beta ~ Sex + seconds_since_midnight + sin + cos, data = leye.m)

# temporal
m9 <- lm(Beta ~ Julian + seconds_since_midnight + sin + cos, data = leye.m)
m10 <- lm(Beta ~ Event + seconds_since_midnight + sin + cos, data = leye.m)
m11 <- lm(Beta ~ seconds_since_midnight + sin + cos, data = leye.m)

# flock
m12 <- lm(Beta ~ Max_Flock_Size + seconds_since_midnight + sin + cos, data = leye.m)


model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)



# Is there enough support to include the random effect? No ----
m1 <- lm(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + sin + cos, data = leye.m)
m2 <- lmer(Beta ~ Dist_Closest_Wetland_m + seconds_since_midnight + sin + cos + (1|Site),
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




# -------------------------------- #



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
  geom_point(data = leye, aes(x = Dist_Closest_Wetland_m, y = Beta), size = 2,
             col = "goldenrod3")


