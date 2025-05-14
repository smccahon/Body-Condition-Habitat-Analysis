#----------------------------------------------------#
# Lesser Yellowlegs Pectoral Muscle Habitat Analysis #
#                linear regression                   #
#               Created 2025-04-30                   #
#              Modified 2025-04-30                   #
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
                      levels = c("Spring 2022",
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
                     levels = c("Spring 2022", "Fall 2023"),
                     labels = c("Spring 2022", "Fall 2023"))

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

# standardize data except for response
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("PecSizeBest"), scale))


# Average pectoral size per site
# leye <- leye %>% 
#   group_by(Site) %>% 
#   summarize(Avg_PecSizeBest = mean(PecSizeBest), PercentAg = first(PercentAg))
# 
# m <- lm(Avg_PecSizeBest ~ PercentAg, data = leye)
# summary(m)
# confint(m)
# 
# d <- expand.grid(PercentAg = seq(min(leye$PercentAg),
#                                  max(leye$PercentAg),
#                                  length = 1000)) 
# 
# d <- cbind(d, predict(m, newdata = d, interval = "confidence"))
# 
# # days into season
# ggplot(d, aes(x = PercentAg, y = fit)) +
#   geom_line(size = 0.8) + 
#   geom_ribbon(aes(ymin = lwr, ymax = upr), 
#               alpha = 0.25, color = NA, show.legend = FALSE) +
#   theme_classic() +
#   labs(x = "% Surrounding Agriculture within 500 m", 
#        y = expression("Lesser Yellowlegs Average Pectoral Muscle Size" ~~~ (mm[score]))) +
#   theme(axis.title.x = element_text(size =18,
#                                     margin = margin(t = 12)),
#         axis.title.y = element_text(size = 18,
#                                     margin = margin(r = 12)),
#         axis.text.x = element_text(size = 18),
#         axis.text.y = element_text(size = 18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.position = "none") +
#   geom_point(data = leye, aes(x = PercentAg, y = Avg_PecSizeBest), size = 2) +
#   scale_color_viridis(alpha = 1, begin = 0, end = 1)

# Model Selection --------------------------------------------------------------

# remove sites with only one observation ----
leye.m <- leye.cs %>% 
  group_by(Site) %>% 
  filter(n() > 1) %>% 
  ungroup()

# agriculture
m1 <- lmer(PecSizeBest ~ PercentAg + (1|Site), data = leye.m, REML = FALSE)

# vegetation
m2 <- lmer(PecSizeBest ~ Percent_Total_Veg + (1|Site), data = leye.m, REML = FALSE)

# habitat
m3 <- lmer(PecSizeBest ~ Permanence + (1|Site), data = leye.m, REML = FALSE)
m4 <- lmer(PecSizeBest ~ Percent_Exposed_Shoreline + (1|Site), data = leye.m, REML = FALSE)
m5 <- lmer(PecSizeBest ~ Dist_Closest_Wetland_m + (1|Site), data = leye.m, REML = FALSE)

# weather
m6 <- lmer(PecSizeBest ~ SPEI + (1|Site), data = leye.m, REML = FALSE)

# life history
m7 <- lmer(PecSizeBest ~ Age + (1|Site), data = leye.m, REML = FALSE)
m8 <- lmer(PecSizeBest ~ Sex + (1|Site), data = leye.m, REML = FALSE)

# temporal
m9 <- lmer(PecSizeBest ~ Julian + (1|Site), data = leye.m, REML = FALSE)
m10 <- lmer(PecSizeBest ~ seconds_since_midnight + (1|Site), data = leye.m, REML = FALSE)
m11 <- lmer(PecSizeBest ~ Event + (1|Site), data = leye.m, REML = FALSE)

# flock
m12 <- lmer(PecSizeBest ~ Max_Flock_Size + (1|Site), data = leye.m, REML = FALSE)

# null
m13 <- lmer(PecSizeBest ~ 1 + (1|Site), data = leye.m, REML = FALSE)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m7)
confint(m7) # age is important

summary(m1)
confint(m1) # ag is important

summary(m9)
confint(m9) # date is important

# do these results hold without random effect on reduced dataset? yes ----

# agriculture
m1 <- lm(PecSizeBest ~ PercentAg, data = leye.m)

# vegetation
m2 <- lm(PecSizeBest ~ Percent_Total_Veg, data = leye.m)

# habitat
m3 <- lm(PecSizeBest ~ Permanence, data = leye.m)
m4 <- lm(PecSizeBest ~ Percent_Exposed_Shoreline, data = leye.m)
m5 <- lm(PecSizeBest ~ Dist_Closest_Wetland_m, data = leye.m)

# weather
m6 <- lm(PecSizeBest ~ SPEI, data = leye.m)

# life history
m7 <- lm(PecSizeBest ~ Age, data = leye.m)
m8 <- lm(PecSizeBest ~ Sex, data = leye.m)

# temporal
m9 <- lm(PecSizeBest ~ Julian, data = leye.m)
m10 <- lm(PecSizeBest ~ seconds_since_midnight, data = leye.m)
m11 <- lm(PecSizeBest ~ Event, data = leye.m)

# flock
m12 <- lm(PecSizeBest ~ Max_Flock_Size, data = leye.m)

# null
m13 <- lm(PecSizeBest ~ 1, data = leye.m)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Stage 2: Multiple combinations with informative parameters--------------------
# excluding variables with correlations
m1 <- lmer(PecSizeBest ~ PercentAg + (1|Site), data = leye.m, REML = FALSE)
m2 <- lmer(PecSizeBest ~ Age + (1|Site), data = leye.m, REML = FALSE)
m3 <- lmer(PecSizeBest ~ Julian + (1|Site), data = leye.m, REML = FALSE)

m4 <- lmer(PecSizeBest ~ PercentAg + Age + (1|Site), data = leye.m, REML = FALSE)
m5 <- lmer(PecSizeBest ~ PercentAg + Julian + (1|Site), data = leye.m, REML = FALSE)
m6 <- lmer(PecSizeBest ~ Age + Julian + (1|Site), data = leye.m, REML = FALSE)

m7 <- lmer(PecSizeBest ~ PercentAg + Julian + Age + (1|Site), data = leye.m, REML = FALSE)

model_names <- paste0("m", 1:7)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m4)
confint(m4) # age and percent ag top model

summary(m2)
confint(m2)

summary(m1)
confint(m1)

summary(m7)
confint(m7)

# second stage without random effect
m1 <- lm(PecSizeBest ~ PercentAg, data = leye.m)
m2 <- lm(PecSizeBest ~ Age, data = leye.m)
m3 <- lm(PecSizeBest ~ Julian, data = leye.m)

m4 <- lm(PecSizeBest ~ PercentAg + Age, data = leye.m)
m5 <- lm(PecSizeBest ~ PercentAg + Julian, data = leye.m)
m6 <- lm(PecSizeBest ~ Age + Julian, data = leye.m)

m7 <- lm(PecSizeBest ~ PercentAg + Julian + Age, data = leye.m)

model_names <- paste0("m", 1:7)

models <- mget(model_names)

aictab(models, modnames = model_names)


# Is there enough support to include the random effect (stage 2)? No ----
m1 <- lm(PecSizeBest ~ PercentAg + Age, data = leye.m)
m2 <- lmer(PecSizeBest ~ PercentAg + Age + (1|Site), data = leye.m, REML = FALSE)


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

# plotting
m <- lmer(PecSizeBest ~ PercentAg + Age + (1|Site), data = leye, REML = FALSE)

d <- expand.grid(Age = c("Juvenile", "Adult"),
                 Site = unique(leye$Site),
                 PercentAg = seq(min(leye$PercentAg),
                                 max(leye$PercentAg),
                                 length = 1000)) 

predictions <- predict(m, newdata = d, se.fit = TRUE, re.form = NA)

d$fit <- predictions$fit

d$lwr <- d$fit - 1.96 * predictions$se.fit
d$upr <- d$fit + 1.96 * predictions$se.fit

ggplot(d, aes(x = (PercentAg), y = fit, color = Age)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Age), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = expression("Lesser Yellowlegs Pectoral Muscle Size" ~~~ (mm[score])),
       color = "Age") +
  theme(axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        axis.title.y = element_text(size = 21, margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = "top") +
  geom_point(data = leye, aes(x = PercentAg, y = PecSizeBest, 
                              color = Age), size = 2) +
  scale_color_manual(values = c("Juvenile" = "#009E73", 
                                "Adult" = "#CC79A7")) +
  scale_fill_manual(values = c("Juvenile" = "#009E73",  
                               "Adult" = "#CC79A7")) 
