#-----------------------------------------#
#  Lesser Yellowlegs Fat Habitat Analysis #
#          Linear regression              #
#          Created 2025-04-23             #
#         Modified 2025-05-29             #
#-----------------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(trtools)
library(lme4)
library(car)
library(viridis)
library(lubridate)
library(nlme)

options(digits = 3)

# read data
birds <- read.csv("Body_Condition_Habitat_Analysis_2025-05-29.csv")

# ...make new columns ----

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

# standardize data except for response
leye.cs <- leye %>%
  mutate(across(where(is.numeric) & !matches("Fat"), scale))

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
# % veg and dominant crop
# permanence & SPEI
# julian and detection (0.60 exactly)

# interaction between capture time and event? yes, seems reasonable
m1 <- lm(Fat ~ Event + seconds_since_midnight, data = leye)
m2 <- lm(Fat ~ Event * seconds_since_midnight, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# interaction between drought and ag?
m1 <- lm(Fat ~ SPEI + Julian, data = leye)
m2 <- lm(Fat ~ PercentAg + Julian, data = leye)
m3 <- lm(Fat ~ PercentAg * SPEI + Julian, data = leye)
m4 <- lm(Fat ~ PercentAg + SPEI + Julian, data = leye)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# interaction between julian and event? no
m1 <- lm(Fat ~ Event + Julian, data = leye)
m2 <- lm(Fat ~ Event * Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)


# Model Selection (Stage 1) ----------------------------------------------------

# agriculture
m1 <- lm(Fat ~ PercentAg, data = leye.cs)

# vegetation
m2 <- lm(Fat ~ Percent_Total_Veg, data = leye.cs)

# habitat
m3 <- lm(Fat ~ Permanence, data = leye.cs)
m4 <- lm(Fat ~ Percent_Exposed_Shoreline, data = leye.cs)
m5 <- lm(Fat ~ Dist_Closest_Wetland_m, data = leye.cs)

# weather
m6 <- lm(Fat ~ SPEI, data = leye.cs)

# life history
m7 <- lm(Fat ~ Age, data = leye.cs)
m8 <- lm(Fat ~ Sex, data = leye.cs)

# temporal
m9 <- lm(Fat ~ Julian, data = leye.cs)
m10 <- lm(Fat ~ seconds_since_midnight, data = leye.cs)
m11 <- lm(Fat ~ Event, data = leye.cs)
m12 <- lm(Fat ~ Event * seconds_since_midnight,
           data = leye.cs)

# flock
m13 <- lm(Fat ~ Max_Flock_Size, data = leye.cs)

# null
m14 <- lm(Fat ~ 1, data = leye.cs)


model_names <- paste0("m", 1:14)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m9)
confint(m9) # fat increases as date progresses

summary(m11)
confint(m11) # fat lower in fall 2021

summary(m12)
confint(m12) # interaction between capture time and season is significant

## Date is extremely important: include it in every model as an informed null

# informative parameters -------------------------------------------------------
# julian, event, event * time

# Plot Julian ----
m <- lm(Fat ~ Julian, data = leye)

d <- expand.grid(Julian = seq(min(leye$Julian),
                            max(leye$Julian),
                            length = 1000)) 

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = Julian, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Julian Date", 
       y = "Lesser Yellowlegs Fat Score") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) +
  geom_point(data = leye, aes(x = Julian, y = Fat), 
             size = 2, width = 0.1, height = 0.05)



# Plot SPEI: No longer significant when treated as numeric ----
m <- lm(Fat ~ SPEI, data = leye)

d <- expand.grid(SPEI = seq(min(leye$SPEI),
                            max(leye$SPEI),
                            length = 1000)) 

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Drought Condition Index", 
       y = "Lesser Yellowlegs Fat Score") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) +
  geom_point(data = leye, aes(x = SPEI, y = Fat), 
              size = 2, width = 0.1, height = 0.05)

# Plot Event ----
m <- lm(Fat ~ Event, data = leye)

d <- expand.grid(Event = unique(leye$Event))

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = Event, y = fit)) +
  geom_point(size = 5, col = "black") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1,
                col = "black",
                size = 1) +
  theme_classic() +
  labs(x = "Sampling Event", 
       y = "Lesser Yellowlegs Fat Score") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top")


# Add Julian as informed null (Report in Thesis as Stage 1) ----

# agriculture
m1 <- lm(Fat ~ PercentAg + Julian, data = leye.cs)

# vegetation
m2 <- lm(Fat ~ Percent_Total_Veg + Julian, data = leye.cs)

# habitat
m3 <- lm(Fat ~ Permanence + Julian, data = leye.cs)
m4 <- lm(Fat ~ Percent_Exposed_Shoreline + Julian, data = leye.cs)
m5 <- lm(Fat ~ Dist_Closest_Wetland_m + Julian, data = leye.cs)

# weather
m6 <- lm(Fat ~ SPEI + Julian, data = leye.cs)

# life history
m7 <- lm(Fat ~ Age + Julian, data = leye.cs)
m8 <- lm(Fat ~ Sex + Julian, data = leye.cs)

# temporal
m9 <- lm(Fat ~ Julian, data = leye.cs)
m10 <- lm(Fat ~ seconds_since_midnight + Julian, data = leye.cs)
m11 <- lm(Fat ~ Event + Julian, data = leye.cs)
m12 <- lm(Fat ~ Event * seconds_since_midnight + Julian,
          data = leye.cs)

# flock
m13 <- lm(Fat ~ Max_Flock_Size + Julian, data = leye.cs)

model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m6)
confint(m6)

# Plot Top Model with Julian as Informed Null ----

m <- lm(Fat ~ SPEI + Julian, data = leye)

d <- expand.grid(SPEI = seq(min(leye$SPEI),
                            max(leye$SPEI),
                            length = 1000),
                 Julian = mean(leye$Julian)) 

predictions <- predict(m, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(size = 0.8) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, color = NA, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Drought Condition Index", 
       y = "Lesser Yellowlegs Fat Score") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) +
  geom_point(data = leye, aes(x = SPEI, y = Fat), 
             size = 2, width = 0.1, height = 0.05)
  

#--------------#

## Identify optimal variance structure ####
m <- lm(Fat ~ SPEI + Julian, data = leye)
plot(predict(m),rstudent(m))

leye$yhat <- predict(m)
leye$rese <- residuals(m)
leye$rest <- rstudent(m)

ggplot(leye, aes(x = yhat, y = rest)) + geom_segment(aes(x = yhat, 
                                                         xend = yhat,
                                                         y = 0,
                                                         yend = rest)) +
  geom_point() + theme_classic() + labs(x = "Predicted Value",
                                        y = "Studentized Residual") +
  geom_hline(yintercept = c(-2,0,2), linetype = 2)

fitted_vals <- fitted(m)  
residuals_vals <- resid(m)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

### Variance structures ####

# VarFixed (Fixed variance): 
m.varFixed <- gls(Fat ~ SPEI + Julian,
                  weights = varFixed(~ Julian),
                  data = leye)

plot(m.varFixed)

# Did it fix our issue? Nope
fitted_vals <- fitted(m.varFixed)  
residuals_vals <- resid(m.varFixed)  
plot(fitted_vals, residuals_vals, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# weighted least squares
leye %>% summarize(variance = var(Fat), weight = 1/var(Fat))

# put results back into data frame
leye <- leye %>% mutate(w = 1/var(Fat))

m.wls <- lm(Fat ~ Julian + SPEI, weights = w, data = leye)

# did not fix 
plot(predict(m.wls),rstudent(m.wls))
plot(predict(m),rstudent(m))

confint(m)
confint(m.wls)

# iteratively weighted least squares
# calculate weighted last squares
# m.ols <- lm(Fat ~ SPEI + Julian, data = leye.cs)
# 
# leye.cs$w <- 1/predict(m.ols)
# m.wls <- lm(Fat ~ SPEI + Julian, weights = w, data = leye.cs)
# plot(predict(m.wls),rstudent(m.wls)) # hint of a megaphone
# 
# leye.cs$w <- 1/predict(m.ols)^2
# m.wls <- lm(Fat ~ SPEI + Julian, weights = w, data = leye.cs)
# plot(predict(m.wls),rstudent(m.wls)) # worse
# 
# # repeat until variance stops changing
# m.ols <- lm(Fat ~ SPEI + Julian, data = leye.cs)
# leye.cs$w <- 1
# 
# for (i in 1:10) {
#   m.wls <- lm(Fat ~ SPEI + Julian, weights = w, data = leye.cs)
#   print(coef(m.wls))
#   leye.cs$w <- 1 / predict(m.wls)
# }
# 
# plot(predict(m.wls), rstudent(m.wls))
# confint(m.wls)
# summary(m.wls)
# 
# leye$yhat <- predict(m.wls)
# leye$rest <- rstudent(m.wls)
# 
# ggplot(leye, aes(x = yhat, y = rest)) + geom_point() + theme_classic() +
#   ylim(-3,3)
# 
# abline(h = 0, col = "red")
# 
# hist(fitted(m.wls), breaks = 30, main = "Histogram of Fitted Values")
# 
# # assess model fit
# ggplot(leye, aes(x = yhat, y = Fat)) +
#   geom_point() + 
#   geom_abline(slope = 1, intercept = 0, color = "blue") +
#   theme_classic() +
#   labs(x = "Fitted Fat", y = "Observed Fat")

# Plot WLS Model ----
d <- expand.grid(SPEI = seq(min(leye$SPEI),
                            max(leye$SPEI),
                            length = 1000),
                 Julian = mean(leye$Julian)) 

predictions <- predict(m.wls, newdata = d, type = "response", se.fit = TRUE) 

d$fit <- predictions$fit

d$lwr <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI
d$upr <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI

ggplot(d, aes(x = SPEI, y = fit)) +
  geom_line(size = 0.8, col = "seagreen") + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.25, fill = "seagreen", show.legend = FALSE) +
  theme_classic() +
  labs(x = "Drought Index (SPEI)", 
       y = "Lesser Yellowlegs Fat Score") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) +
  geom_point(data = leye, aes(x = SPEI, y = Fat), 
             size = 2.5, col = "seagreen") +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,5, by = 1)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(-2, 2, by = 0.5)) +
  coord_cartesian(ylim = c(-0.05,5.25),
                  xlim = c(-1.9,2.25))
  # geom_vline(xintercept = 0, linetype = "dashed", size = 1, color = "gray")



summary(m.wls)
confint(m.wls)


# add macroinvertebrate diversity and biomass data for 2023 --------------------
birds.sub <- subset(leye.cs, !is.na(Biomass))

m <- lm(Fat ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lm(Fat ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of diversity



# do neonics explain further variation in fat? ----
# informative covariates: SPEI + Julian


# summary statistics----
table(leye$PlasmaDetection) # n: 31, y: 23 (n = 54)
table(leye$WaterNeonicDetection) # n: 39, y: 14 (n = 53)
table(leye$AnyDetection) # n: 17, y: 37 (n = 54)
table(leye$WaterOrInvertDetection) # n: 25, y: 29 (n = 54)
table(leye$InvertPesticideDetection) # n: 15, y: 15 (n = 30)

mean(leye$OverallNeonic, na.rm = TRUE) # 2.39 ug/L
sd(leye$OverallNeonic, na.rm = TRUE) # 10.2 ug/L

# water neonic detection --> neonics not informative
leye.clean.water <- leye.cs[!is.na(leye.cs$WaterNeonicDetection), ] #n = 53

m1 <- lm(Fat ~ SPEI + Julian, data = leye.clean.water)
m2 <- lm(Fat ~ SPEI + Julian + WaterNeonicDetection, data = leye.clean.water)

# invertebrate pesticide detection --> neonics not informative
leye.clean.invert <- leye.cs[!is.na(leye.cs$InvertPesticideDetection), ] #n = 30

m1 <- lm(Fat ~ SPEI + Julian, data = leye.clean.invert)
m2 <- lm(Fat ~ SPEI + Julian + InvertPesticideDetection, data = leye.clean.invert)

# invertebrate or water pesticide detection (environmental detection) --> neonics not informative
leye.clean.waterorinvert <- leye.cs[!is.na(leye.cs$WaterOrInvertDetection), ] #n = 54
m1 <- lm(Fat ~ SPEI + Julian, data = leye.clean.waterorinvert)
m2 <- lm(Fat ~ SPEI + Julian + WaterOrInvertDetection, data = leye.clean.waterorinvert)

# shorebird plasma detection --> neonics not informative
m1 <- lm(Fat ~ SPEI + Julian, data = leye.cs)
m2 <- lm(Fat ~ SPEI + Julian + PlasmaDetection, data = leye.cs)

# any detection (plasma or environmental) --> neonics not informative
m1 <- lm(Fat ~ SPEI + Julian, data = leye.cs)
m2 <- lm(Fat ~ SPEI + Julian + AnyDetection, data = leye.cs)

### ...AIC 
models <- list(m1, m2)
model.sel(models)

# model summaries:
summary(m2)
confint(m2)

