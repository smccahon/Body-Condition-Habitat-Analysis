#----------------------------------------------------#
# Lesser Yellowlegs Refueling Rates Habitat Analysis #
#             Non-linear regression                  #
#               Created 2025-04-01                   #
#              Modified 2025-04-02                   #
#----------------------------------------------------#

# load packages
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(trtools)
library(lme4)
library(car)
library(viridis)


# read data
birds <- read.csv("Body_Condition_Habitat_Analysis_2025-03-31.csv")

# Transform Data ---------------------------------------------------------------

# ...make new columns ----
# neonicotinoid detection column
birds$Detection <- ifelse(birds$OverallNeonic > 0, 
                          "Detection", "Non-detection")

# convert capture time to number of minutes after sunrise
birds <- birds %>%
  mutate(
    hour_of_day = as.numeric(format(strptime(Time, "%H:%M"), "%H")),  # Extract hour (0-23)
    minute_of_day = as.numeric(format(strptime(Time, "%H:%M"), "%M")),  # Extract minute (0-59)
    time_in_minutes = hour_of_day * 60 + minute_of_day  # Total time in minutes
  ) %>% 
  mutate(
    hour_of_day_s = as.numeric(format(strptime(Sunrise, "%H:%M"), "%H")),  # Extract hour (0-23)
    minute_of_day_s = as.numeric(format(strptime(Sunrise, "%H:%M"), "%M")),  # Extract minute (0-59)
    time_in_minutes_s = hour_of_day_s * 60 + minute_of_day_s  # Total time in minutes
  ) %>% 
  mutate(
    ts.sunrise = time_in_minutes - time_in_minutes_s)

# ...reorder and manipulate relevant factor variables ----
birds$Sex <- factor(birds$Sex,
                    levels = c("M", "F"),
                    labels = c("Male", "Female"))

birds$Buffered <- factor(birds$Buffered,
                    levels = c("Y", "N"),
                    labels = c("Buffered", "Not Buffered"))

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

# filter birds that only contain metabolite information (n = 29)
leye <- leye %>% 
  filter(!is.na(Tri) & !is.na(Beta))

# is there enough variation in buffer presence?
table(leye$Buffered) # nope
table(birds$Buffered) # also nope

# Perform PCA ------------------------------------------------------------------

# Subset the dataset to only include 'Tri' and 'Beta'
leye_subset <- leye[, c("Tri", "Beta")]

# Remove rows with NAs for PCA
leye_subset_clean <- leye_subset[complete.cases(leye_subset), ]

# Run PCA on the cleaned data
pca_result <- prcomp(leye_subset_clean, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye$PC1 <- NA  # Initialize with NA values
leye$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]

# View the updated dataset with PCA scores
print(leye)

# Merge PCA scores back to the original LEYE unstandardized dataset
# Create a data frame to store the PCA scores for the rows with no missing data
leye$PC1 <- NA  # Initialize with NA values
leye$PC2 <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
leye[complete.cases(leye_subset), "PC1"] <- pca_scores[, 1]
leye[complete.cases(leye_subset), "PC2"] <- pca_scores[, 2]


# which correlated covariates should I drop? -----------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg, data = leye)
m2 <- lm(PC1 ~ DominantCrop, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# % ag is a much better predictor

# weather and time
m1 <- lm(PC1 ~ SPEI, data = leye)
m2 <- lm(PC1 ~ Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# if included in the same model, Julian better predictor

# habitat and vegetation
m1 <- lm(PC1 ~ Percent_Total_Veg, data = leye)
m2 <- lm(PC1 ~ Percent_Exposed_Shoreline, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# if included in the same model, % exposed shoreline better predictor

# which data variable to use for LEYE? -----------------------------------------
m1 <- lm(PC1 ~ DaysIntoSeason, data = leye)
m2 <- lm(PC1 ~ Julian, data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Julian is better

# aside from time, what else is important? -------------------------------------

# agriculture
m1 <- lm(PC1 ~ PercentAg, data = leye)

# vegetation
m2 <- lm(PC1 ~ Percent_Total_Veg, data = leye)

# habitat
m3 <- lm(PC1 ~ Permanence, data = leye)
m4 <- lm(PC1 ~ Percent_Exposed_Shoreline, data = leye)
m5 <- lm(PC1 ~ Dist_Closest_Wetland_m, data = leye)

# weather
m6 <- lm(PC1 ~ SPEI, data = leye)

# life history
m7 <- lm(PC1 ~ Age, data = leye)
m8 <- lm(PC1 ~ Sex, data = leye)

# temporal
m9 <- lm(PC1 ~ Julian, data = leye)
m10 <- lm(PC1 ~ Event, data = leye)
m11 <- lm(PC1 ~ ts.sunrise, data = leye)

# flock
m12 <- lm(PC1 ~ Max_Flock_Size, data = leye)

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

## wetland permanence is better than linear transformation of time
## percent ag is also significant (not within 2 delta AICc)

# are any variables important after transforming time into a nonlinear variable? ----

# ...cos(time) ----
# agriculture
m1 <- nls(PC1 ~ b0 + b1*PercentAg, data = leye,
          start = list(b0 = 0, b1 = 0))

# vegetation
m2 <- nls(PC1 ~ b0 + b1*Percent_Total_Veg, data = leye,
         start = list(b0 = 0, b1 = 0))

# habitat
leye$Permanence <- relevel(leye$Permanence, ref = "Temporary/Seasonal")
m3 <- nls(PC1 ~ b0 + b1*(Permanence == "Semipermanent") + 
            b2*(Permanence == "Permanent"), data = leye,
         start = list(b0 = 0, b1 = 0, b2 = 0))

m4 <- nls(PC1 ~ b0 + b1*Percent_Exposed_Shoreline, data = leye,
         start = list(b0 = 0, b1 = 0))

m5 <- nls(PC1 ~ b0 + b1*Dist_Closest_Wetland_m, data = leye,
         start = list(b0 = 0, b1 = 0))

# weather
m6 <- nls(PC1 ~ b0 + b1*SPEI, data = leye,
         start = list(b0 = 0, b1 = 0))

# life history
leye$Age <- relevel(leye$Age, ref = "Adult")
m7 <- nls(PC1 ~ b0 + b1*(Age == "Juvenile"), data = leye,
         start = list(b0 = 0, b1 = 0))

leye$Sex <- relevel(leye$Sex, ref = "Male")
m8 <- nls(PC1 ~ b0 + b1*(Sex == "Female"), data = leye,
         start = list(b0 = 0, b1 = 0))

# temporal
m9 <- nls(PC1 ~ b0 + b1*Julian, data = leye,
          start = list(b0 = 0, b1 = 0))

leye$Event <- relevel(leye$Event, ref = "Fall 2021")
m10 <- nls(PC1 ~ b0 + b1*(Event == "Spring 2022") + 
             b2*(Event == "Fall 2023"), data = leye,
          start = list(b0 = 0, b1 = 0, b2 = 0))

m11 <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n), 
          data = leye,
          start = list(beta0 = 0, beta1 = 2, a = 113, n = -1.8),
          control = nls.control(maxiter = 500))

# flock
m12 <- nls(PC1 ~ b0 + b1*Max_Flock_Size, data = leye,
          start = list(b0 = 0, b1 = 0))

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# no variables have ANY effect beyond cos(time)

# ...time^3 ----
# agriculture
m1 <- nls(PC1 ~ b0 + b1*PercentAg, data = leye,
          start = list(b0 = 0, b1 = 0))

# vegetation
m2 <- nls(PC1 ~ b0 + b1*Percent_Total_Veg, data = leye,
          start = list(b0 = 0, b1 = 0))

# habitat
leye$Permanence <- relevel(leye$Permanence, ref = "Temporary/Seasonal")
m3 <- nls(PC1 ~ b0 + b1*(Permanence == "Semipermanent") + 
            b2*(Permanence == "Permanent"), data = leye,
          start = list(b0 = 0, b1 = 0, b2 = 0))

m4 <- nls(PC1 ~ b0 + b1*Percent_Exposed_Shoreline, data = leye,
          start = list(b0 = 0, b1 = 0))
m5 <- nls(PC1 ~ b0 + b1*Dist_Closest_Wetland_m, data = leye,
          start = list(b0 = 0, b1 = 0))

# weather
m6 <- nls(PC1 ~ b0 + b1*SPEI, data = leye,
          start = list(b0 = 0, b1 = 0))

# life history
leye$Age <- relevel(leye$Age, ref = "Adult")
m7 <- nls(PC1 ~ b0 + b1*(Age == "Juvenile"), data = leye,
          start = list(b0 = 0, b1 = 0))

leye$Sex <- relevel(leye$Sex, ref = "Male")
m8 <- nls(PC1 ~ b0 + b1*(Sex == "Female"), data = leye,
          start = list(b0 = 0, b1 = 0))

# temporal
m9 <- nls(PC1 ~ b0 + b1*Julian, data = leye,
           start = list(b0 = 0, b1 = 0))

leye$Event <- relevel(leye$Event, ref = "Fall 2021")
m10 <- nls(PC1 ~ b0 + b1*(Event == "Spring 2022") + 
             b2*(Event == "Fall 2023"), data = leye,
           start = list(b0 = 0, b1 = 0, b2 = 0))

m11 <- nls(PC1 ~ beta0 + beta1*ts.sunrise + beta2*I(ts.sunrise^2) +
            beta3*I(ts.sunrise^3), 
          data = leye,
          start = list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0))

# flock
m12 <- nls(PC1 ~ b0 + b1*Max_Flock_Size, data = leye,
           start = list(b0 = 0, b1 = 0))

model_names <- paste0("m", 1:12)

models <- mget(model_names)

aictab(models, modnames = model_names)

# no variables have any effect beyond time^3

# should time be linear or non-linear? -----------------------------------------
# plot on desmos to find starting values
# leye_sub <- leye[, c("ts.sunrise", "PC1")]

m1 <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n), 
          data = leye,
          start = list(beta0 = 0, beta1 = 2, a = 113, n = -1.8),
          control = nls.control(maxiter = 500))

m2 <- nls(PC1 ~ beta0 + beta1*ts.sunrise, 
          data = leye,
          start = list(beta0 = 0, beta1 = 0))

m3 <- nls(PC1 ~ beta0 + beta1*ts.sunrise + beta2*I(ts.sunrise^2) +
            beta3*I(ts.sunrise^3), 
          data = leye,
          start = list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0))

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

summary(m3)

# ...non-linear cosine relationship by far the most informative-----------------

# plot on desmos
# leye_sub <- leye[, c("ts.sunrise", "PC1")]

# plot model of time with cosine transformation---------------------------------
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n), data = leye,
         start = list(beta0 = 0, beta1 = 2, a = 113, n = -1.8))

summary(m)
confint(m)

d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise), 
                                  max(leye$ts.sunrise), 
                                  length.out = 1000),  
                 PercentAg = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d <- cbind(d, nlsint(m, newdata = d))

ggplot(d, aes(x = ts.sunrise, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "Time of Capture Since Sunrise (min)", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2, 
           label = "Sunrise", angle = 90, 
           vjust = -1, hjust = 1,
           size = 5) +
  geom_point(data = leye, aes(x = ts.sunrise, y = PC1), size = 3)

# plot model of time with quadratic transformation------------------------------
m <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2) + I(ts.sunrise^3), 
        data = leye)

d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise), 
                                  max(leye$ts.sunrise), 
                                  length.out = 1000))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = ts.sunrise, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Time of Capture Since Sunrise (min)", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2, 
           label = "Sunrise", angle = 90, 
           vjust = -1, hjust = 1,
           size = 5) +
  geom_point(data = leye, aes(x = ts.sunrise, y = PC1), size = 3)

# plot model of time with linear transformation---------------------------------
m <- lm(PC1 ~ ts.sunrise,  data = leye)


d <- expand.grid(ts.sunrise = seq(min(leye$ts.sunrise), 
                                  max(leye$ts.sunrise), 
                                  length.out = 1000))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = ts.sunrise, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Time of Capture Since Sunrise (min)", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 0, y = 2, 
           label = "Sunrise", angle = 90, 
           vjust = -1, hjust = 1,
           size = 5) +
  geom_point(data = leye, aes(x = ts.sunrise, y = PC1), size = 3)

# should % ag be linear or nonlinear? ------------------------------------------
m1 <- nls(PC1 ~ beta0 + beta1 * -exp(PercentAg/a), 
          data = leye,
          start = list(beta0 = 2, beta1 = 2, a = 40),
          control = nls.control(maxiter = 500))

m2 <- nls(PC1 ~ beta0 + beta1*PercentAg, 
          data = leye,
          start = list(beta0 = 0, beta1 = 0))

m3 <- nls(PC1 ~ beta0 + beta1*PercentAg + beta2*I(PercentAg^2), 
          data = leye,
          start = list(beta0 = 0, beta1 = 0, beta2 = 0))

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

# quadratic relationship is best

summary(m1) # not significant
summary(m2) # significant
summary(m3) # significant   

confint(m3)
confint(m2)
trtools::lincon(m1, fcov=vcov)

# ...plot model of ag with exponential decay transformation---------------------
m <- nls(PC1 ~ beta0 + beta1 * -exp(PercentAg/a), 
         data = leye,
         start = list(beta0 = 2, beta1 = 2, a = 40),
         control = nls.control(maxiter = 1000))

summary(m)
nlsint(m)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length.out = 1000),  
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d <- cbind(d, nlsint(m, newdata = d))

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), size = 3)


# ...plot model of ag with linear transformation---------------------
m <- lm(PC1 ~ PercentAg, data = leye)

summary(m)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length.out = 1000))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), size = 3)


# ...plot model of ag with quadratic transformation-----------------------------
m <- lm(PC1 ~ PercentAg + I(PercentAg^2), data = leye)

summary(m)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length.out = 1000))

predictions <- predict(m, newdata = d)

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), size = 3)

# ...model % ag with different transformations of time--------------------------

# cos(time) + % surrounding ag (exponential)
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n) + 
           beta2 * -exp(PercentAg/b), 
         data = leye,
         start = list(beta0 = 5, beta1 = 4, a = 113, n = -1.8,
                      beta2 = 3, b = 55),
         control = nls.control(maxiter = 1000))

# model too complex and won't converge

# cos(time) + % surrounding ag (linear)
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n) + 
           beta2 * PercentAg, 
         data = leye,
         start = list(beta0 = 5, beta1 = 4, a = 113, n = -1.8,
                      beta2 = 0),
         control = nls.control(maxiter = 1000))

summary(m)

# % ag not significant

# cos(time) + % ag (quadratic)
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n) + 
           beta2 * PercentAg + beta3 * I(PercentAg^2), 
         data = leye,
         start = list(beta0 = 5, beta1 = 4, a = 113, n = -1.8,
                      beta2 = 0, beta3 = 0),
         control = nls.control(maxiter = 1000))

summary(m)
confint(m)

# time and % ag significant


# does cos(time) + "variable" plot properly? answer is no ----------------------
# add vegetation
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n) + 
           beta2*Percent_Total_Veg, 
         data = leye,
         start = list(beta0 = 0, beta1 = 2, a = 113, n = -1.8,
                      beta2 = 0),
         control = nls.control(maxiter = 1000))

summary(m)
confint(m)

d <- expand.grid(Percent_Total_Veg = seq(min(leye$Percent_Total_Veg), 
                                         max(leye$Percent_Total_Veg), 
                                         length.out = 1000),  
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d <- cbind(d, nlsint(m, newdata = d))

ggplot(d, aes(x = Percent_Total_Veg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) 

# add % ag (quadratic)
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n) + 
           beta2*PercentAg + beta3*I(PercentAg^2), 
         data = leye,
         start = list(beta0 = 0, beta1 = 2, a = 113, n = -1.8,
                      beta2 = 0, beta3 = 0),
         control = nls.control(maxiter = 1000))

summary(m)
confint(m)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                  max(leye$PercentAg), 
                                  length.out = 1000),  
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d <- cbind(d, nlsint(m, newdata = d))

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1))

# add % ag (linear)
m <- nls(PC1 ~ beta0 + beta1*cos((ts.sunrise/a) + n) + 
           beta2*PercentAg, 
         data = leye,
         start = list(beta0 = 0, beta1 = 2, a = 113, n = -1.8,
                      beta2 = 0),
         control = nls.control(maxiter = 1000))

summary(m)
confint(m)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length.out = 1000),  
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d <- cbind(d, nlsint(m, newdata = d))

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1))



# does time^3 + "variable" plot properly? answer is yes ----------------------
# add ag (linear)
m <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2) + I(ts.sunrise^3) + PercentAg, 
        data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                  max(leye$PercentAg), 
                                  length.out = 1000),
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), size = 3)

# add ag (quadratic)
m <- lm(PC1 ~ ts.sunrise + I(ts.sunrise^2) + I(ts.sunrise^3) + PercentAg +
          I(PercentAg^2), 
        data = leye)

d <- expand.grid(PercentAg = seq(min(leye$PercentAg), 
                                 max(leye$PercentAg), 
                                 length.out = 1000),
                 ts.sunrise = mean(leye$ts.sunrise))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = PercentAg, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "% Surrounding Agriculture within 500 m", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = PercentAg, y = PC1), size = 3)



# do any other variables need nonlinear transformations? -----------------------

ggplot(data = leye, aes(x = Percent_Total_Veg, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = Dist_Closest_Wetland_m, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = Percent_Exposed_Shoreline, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = DaysIntoSeason, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = Julian, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = MaxBufferWidth_m, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = SPEI, y = PC1)) + geom_point() # no
ggplot(data = leye, aes(x = Max_Flock_Size, y = PC1)) + geom_point() # maybe exponential
ggplot(data = leye, aes(x = PercentAg, y = PC1)) + geom_point() # no


# max flock size transformations
m1 <- nls(PC1 ~ beta0 + beta1 * -exp(Max_Flock_Size/a), 
         data = leye,
         start = list(beta0 = 4, beta1 = 1, a = 30),
         control = nls.control(maxiter = 500))

m2 <- nls(PC1 ~ beta0 + beta1 * Max_Flock_Size, 
         data = leye,
         start = list(beta0 = 0, beta1 = 0),
         control = nls.control(maxiter = 500))

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# linear is best


# standardize time to something more simple ------------------------------------
leye$Time <- strptime(leye$Time, format = "%H:%M")
leye$Time <- as.POSIXct(leye$Time, tz = "America/Chicago")
attributes(leye$Time)$tzone

# Calculate seconds since midnight (start of the day)
leye$seconds_since_midnight <- as.numeric(difftime(leye$Time, 
                                                   floor_date(leye$Time, "day"), 
                                                   units = "secs"))

# Check the first few values to confirm
head(leye$seconds_since_midnight)

# cosine and sine
m <- lm(PC1 ~ seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye)

# just sin
m <- lm(PC1 ~ seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye)

summary(m)


d <- expand.grid(seconds_since_midnight = seq(min(leye$seconds_since_midnight), 
                                  max(leye$seconds_since_midnight), 
                                  length.out = 1000))

predictions <- predict(m, newdata = d)

d$yhat <- predictions

d$yhat <- predict(m, newdata = d)

d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
d$lower <- d$yhat - 2*d$se
d$upper <- d$yhat + 2*d$se

ggplot(d, aes(x = seconds_since_midnight, y = yhat)) +
  geom_line(size = 1) +  
  theme_classic() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Time of Capture", 
       y = "Lesser Yellowlegs Fattening Index") +
  theme(axis.title.x = element_text(size = 21,
                                    margin = margin(t = 12)),
        axis.title.y = element_text(size = 21,
                                    margin = margin(r = 12)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.position = "top") +
  geom_hline(yintercept = 0, linetype = "twodash", color = "red",
             size = 1) +
  geom_point(data = leye, aes(x = seconds_since_midnight, y = PC1), size = 3) +
  scale_x_time(labels = scales::time_format("%H:%M"),
               breaks = seq(0, 86400, by = 7200))

plot(predict(m), rstudent(m))

# do I need both sine and cosine?
m1 <- lm(PC1 ~ seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)) +  
          cos(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye)

# just sin
m2 <- lm(PC1 ~ seconds_since_midnight + 
          sin(2 * pi * seconds_since_midnight / (24 * 3600)),
        data = leye)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# you need both sine and cosine in model

     