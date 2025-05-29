#-------------------------------------------------------#
# All Species Size-Corrected Body Mass Habitat Analysis #
#                Created 2025-04-11                     #
#               Modified 2025-05-13                     #
#-------------------------------------------------------#

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

# Logarithmic transformation of mass
birds <- birds %>% 
  mutate(LogMass = log(Mass))

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

# Only include species with at least three individuals that also have mass, tarsus, wing length, and bill length

birds <- birds %>%
  filter(!is.na(Mass))

birds <- birds %>%
  filter(!is.na(Culmen))

birds <- birds %>%
  filter(!is.na(DiagTarsus))

birds <- birds %>%
  filter(!is.na(Wing))

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# ---------------------------------------------------------------------------- #

# Subset data for Lesser Yellowlegs
leye <- subset(birds, Species %in% c("Lesser Yellowlegs"))

leye.cs <- leye %>%
  mutate(across(where(is.numeric), scale))

leye_female <- subset(leye.cs, Sex == "Female")
leye_male <- subset(leye.cs, Sex == "Male")

leye_female <- leye_female[, c("Wing", "Culmen", "DiagTarsus")]
leye_male <- leye_male[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Least Sandpiper
lesa <- subset(birds, Species %in% c("Least Sandpiper"))

lesa.cs <- lesa %>%
  mutate(across(where(is.numeric), scale))

lesa.cs.subset <- lesa.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Semipalmated Sandpiper
sesa <- subset(birds, Species %in% c("Semipalmated Sandpiper"))

sesa.cs <- sesa %>%
  mutate(across(where(is.numeric), scale))

sesa.cs.subset <- sesa.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Pectoral Sandpiper
pesa <- subset(birds, Species %in% c("Pectoral Sandpiper"))

pesa.cs <- pesa %>%
  mutate(across(where(is.numeric), scale))

pesa.cs.subset <- pesa.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Wilson's Phalarope
wiph <- subset(birds, Species %in% c("Wilsons Phalarope"))

wiph.cs <- wiph %>%
  mutate(across(where(is.numeric), scale))

wiph_female <- subset(wiph.cs, Sex == "Female")
wiph_male <- subset(wiph.cs, Sex == "Male")

wiph_female <- wiph_female[, c("Wing", "Culmen", "DiagTarsus")]
wiph_male <- wiph_male[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for American Avocet
amav <- subset(birds, Species %in% c("American Avocet"))

amav.cs <- amav %>%
  mutate(across(where(is.numeric), scale))

amav.cs.subset <- amav.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Killdeer
kill <- subset(birds, Species %in% c("Killdeer"))

kill.cs <- kill %>%
  mutate(across(where(is.numeric), scale))

kill.cs.subset <- kill.cs[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Willet
will <- subset(birds, Species %in% c("Willet"))

will.cs <- will %>%
  mutate(across(where(is.numeric), scale))

will_female <- subset(will.cs, Sex == "Female")
will_male <- subset(will.cs, Sex == "Male")

will_female <- will_female[, c("Wing", "Culmen", "DiagTarsus")]
will_male <- will_male[, c("Wing", "Culmen", "DiagTarsus")]

# Subset data for Long-billed Dowitcher
lbdo <- subset(birds, Species %in% c("Longbilled Dowitcher"))

lbdo.cs <- lbdo %>%
  mutate(across(where(is.numeric), scale))

lbdo.cs.subset <- lbdo.cs[, c("Wing", "Culmen", "DiagTarsus")]

# ---------------------------------------------------------------------------- #

# Perform PCA ####
## Lesser Yellowlegs: Females ####

# Run PCA on the subsetted data
pca_result_f <- prcomp(leye_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_f)

# Extract residuals 
PC1_leye_female <- pca_result_f$x[, 1]

# Regress PC1 against body mass
m_leye_female <- lm(PC1_leye_female ~ Mass, data = subset(leye.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.leye.f <- resid(m_leye_female)

# Add size-corrected mass back into the full dataset
leye.cs$sc.mass <- NA
leye$sc.mass <- NA
leye.cs[leye.cs$Sex == "Female", "sc.mass"] <- sc.mass.leye.f
leye[leye$Sex == "Female", "sc.mass"] <- sc.mass.leye.f

## Lesser Yellowlegs: Males ####

# Run PCA on the subsetted data
pca_result_m <- prcomp(leye_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_m)

# Extract residuals 
PC1_leye_male <- pca_result_m$x[, 1]

# Regress PC1 against body mass
m_leye_male <- lm(PC1_leye_male ~ Mass, data = subset(leye.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.leye.m <- resid(m_leye_male)

# Add size-corrected mass back into the full datasest
leye.cs[leye.cs$Sex == "Male", "sc.mass"] <- sc.mass.leye.m
leye[leye$Sex == "Male", "sc.mass"] <- sc.mass.leye.m

# Add leye size corrected mass to full dataset

#---------#

## Wilson's Phalaropes: Females ####

# Run PCA on the subsetted data
pca_result_f <- prcomp(wiph_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_f)

# Extract residuals 
PC1_wiph_female <- pca_result_f$x[, 1]

# Regress PC1 against body mass
m_wiph_female <- lm(PC1_wiph_female ~ Mass, data = subset(wiph.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.wiph.f <- resid(m_wiph_female)

# Add size-corrected mass back into the full dataset
wiph.cs$sc.mass <- NA
wiph$sc.mass <- NA
wiph.cs[wiph.cs$Sex == "Female", "sc.mass"] <- sc.mass.wiph.f
wiph[wiph$Sex == "Female", "sc.mass"] <- sc.mass.wiph.f

## Wilson's Phalaropes: Males ####

# Run PCA on the subsetted data
pca_result_m <- prcomp(wiph_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_m)

# Extract residuals 
PC1_wiph_male <- pca_result_m$x[, 1]

# Regress PC1 against body mass
m_wiph_male <- lm(PC1_wiph_male ~ Mass, data = subset(wiph.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.wiph.m <- resid(m_wiph_male)

# Add size-corrected mass back into the full datasest
wiph.cs[wiph.cs$Sex == "Male", "sc.mass"] <- sc.mass.wiph.m
wiph[wiph$Sex == "Male", "sc.mass"] <- sc.mass.wiph.m

#---------#

## Willet: Females ####

# Run PCA on the subsetted data
pca_result_f <- prcomp(will_female, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_f)

# Extract residuals 
PC1_will_female <- pca_result_f$x[, 1]

# Regress PC1 against body mass
m_will_female <- lm(PC1_will_female ~ Mass, data = subset(will.cs, Sex == "Female"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.will.f <- resid(m_will_female)

# Add size-corrected mass back into the full dataset
will.cs$sc.mass <- NA
will$sc.mass <- NA
will.cs[will.cs$Sex == "Female", "sc.mass"] <- sc.mass.will.f
will[will$Sex == "Female", "sc.mass"] <- sc.mass.will.f

## Willet: Males ####

# Run PCA on the subsetted data
pca_result_m <- prcomp(will_male, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result_m)

# Extract residuals 
PC1_will_male <- pca_result_m$x[, 1]

# Regress PC1 against body mass
m_will_male <- lm(PC1_will_male ~ Mass, data = subset(will.cs, Sex == "Male"))

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.will.m <- resid(m_will_male)

# Add size-corrected mass back into the full datasest
will.cs[will.cs$Sex == "Male", "sc.mass"] <- sc.mass.will.m
will[will$Sex == "Male", "sc.mass"] <- sc.mass.will.m

#---------#

## Least Sandpiper ####

# Run PCA on the subsetted data
pca_result <- prcomp(lesa.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_lesa <- pca_result$x[, 1]

# Regress PC1 against body mass
m_lesa <- lm(PC1_lesa ~ Mass, data = lesa.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.lesa <- resid(m_lesa)

# Add size-corrected mass back into the full datasest
lesa.cs$sc.mass <- NA
lesa$sc.mass <- NA
lesa.cs["sc.mass"] <- sc.mass.lesa
lesa["sc.mass"] <- sc.mass.lesa

#---------#

## Semipalmated Sandpiper ####

# Run PCA on the subsetted data
pca_result <- prcomp(sesa.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_sesa <- pca_result$x[, 1]

# Regress PC1 against body mass
m_sesa <- lm(PC1_sesa ~ Mass, data = sesa.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.sesa <- resid(m_sesa)

# Add size-corrected mass back into the full datasest
sesa.cs$sc.mass <- NA
sesa$sc.mass <- NA
sesa.cs["sc.mass"] <- sc.mass.sesa
sesa["sc.mass"] <- sc.mass.sesa

#---------#

## Pectoral Sandpiper ####

# Run PCA on the subsetted data
pca_result <- prcomp(pesa.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_pesa <- pca_result$x[, 1]

# Regress PC1 against body mass
m_pesa <- lm(PC1_pesa ~ Mass, data = pesa.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.pesa <- resid(m_pesa)

# Add size-corrected mass back into the full datasest
pesa.cs$sc.mass <- NA
pesa$sc.mass <- NA
pesa.cs["sc.mass"] <- sc.mass.pesa
pesa["sc.mass"] <- sc.mass.pesa

#---------#

## Killdeer ####

# Run PCA on the subsetted data
pca_result <- prcomp(kill.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_kill <- pca_result$x[, 1]

# Regress PC1 against body mass
m_kill <- lm(PC1_kill ~ Mass, data = kill.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.kill <- resid(m_kill)

# Add size-corrected mass back into the full datasest
kill.cs$sc.mass <- NA
kill$sc.mass <- NA
kill.cs["sc.mass"] <- sc.mass.kill
kill["sc.mass"] <- sc.mass.kill

#---------#

## American Avocet ####

# Run PCA on the subsetted data
pca_result <- prcomp(amav.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_amav <- pca_result$x[, 1]

# Regress PC1 against body mass
m_amav <- lm(PC1_amav ~ Mass, data = amav.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.amav <- resid(m_amav)

# Add size-corrected mass back into the full datasest
amav.cs$sc.mass <- NA
amav$sc.mass <- NA
amav.cs["sc.mass"] <- sc.mass.amav
amav["sc.mass"] <- sc.mass.amav

#---------#

## Long-billed Dowitcher ####

# Run PCA on the subsetted data
pca_result <- prcomp(lbdo.cs.subset, center = TRUE, scale. = TRUE)

# View PCA summary to understand variance explained by principal components
summary(pca_result)

# Extract residuals 
PC1_lbdo <- pca_result$x[, 1]

# Regress PC1 against body mass
m_lbdo <- lm(PC1_lbdo ~ Mass, data = lbdo.cs)

# Save the residuals of the regression to obtain size-corrected mass
sc.mass.lbdo <- resid(m_lbdo)

# Add size-corrected mass back into the full datasest
lbdo.cs$sc.mass <- NA
lbdo$sc.mass <- NA
lbdo.cs["sc.mass"] <- sc.mass.lbdo
lbdo["sc.mass"] <- sc.mass.lbdo

# ---------------------------------------------------------------------------- #

# Combine all species together

birds.cs <- rbind(leye.cs, pesa.cs, lbdo.cs, amav.cs, kill.cs, lesa.cs, 
                  will.cs, sesa.cs, wiph.cs)

birds <- rbind(leye, pesa, lbdo, amav, kill, lesa, 
               will, sesa, wiph)

# standardize data except for response
birds.cs <- birds %>%
  mutate(across(where(is.numeric) & !matches("sc.mass"), scale))

# ---------------------------------------------------------------------------- #

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
# % ag and ag category
# % ag and dominant crop
# SPEI and julian
# ag category and dominant crop

# which correlated variables should I drop? ------------------------------------

# agriculture
m1 <- lm(sc.mass ~ PercentAg, data = birds.cs )
m2 <- lm(sc.mass ~ DominantCrop, data = birds.cs )
m3 <- lm(sc.mass ~ AgCategory, data = birds.cs )

model_names <- paste0("m", 1:3)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)

# % ag is best

# temporal
m1 <- lm(sc.mass ~ Julian + MigStatus, data = birds.cs)
m2 <- lm(sc.mass ~ Julian * MigStatus, data = birds.cs)
m3 <- lm(sc.mass ~ MigStatus, data = birds.cs )
m4 <- lm(sc.mass ~ Julian, data = birds.cs )

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

confint(m1)
confint(m2)
# model with interaction is more informative but not significant, drop i

# transformation for time needed? no
plot(birds$seconds_since_midnight, birds$sc.mass)

# interactions between season and capture time? no
m1 <- lm(sc.mass ~ Event + seconds_since_midnight, data = birds.cs)
m2 <- lm(sc.mass ~ Event * seconds_since_midnight, data = birds.cs)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)

# model without interaction is much better

# interaction between ag and SPEI needed? no
m1 <- lm(sc.mass ~ SPEI + PercentAg, data = birds.cs)
m2 <- lm(sc.mass ~ SPEI * PercentAg, data = birds.cs)
m3 <- lm(sc.mass ~ SPEI, data = birds.cs)
m4 <- lm(sc.mass ~ PercentAg, data = birds.cs)

model_names <- paste0("m", 1:4)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Model Selection with Informed Null (Stage 1) ---------------------------------

# agriculture
m1 <- lm(sc.mass ~ PercentAg, data = birds.cs)

# vegetation
m2 <- lm(sc.mass ~ Percent_Total_Veg, data = birds.cs)

# habitat
m3 <- lm(sc.mass ~ Permanence, data = birds.cs)
m4 <- lm(sc.mass ~ Percent_Exposed_Shoreline, data = birds.cs)
m5 <- lm(sc.mass ~ Dist_Closest_Wetland_m, data = birds.cs)

# weather
m6 <- lm(sc.mass ~ SPEI, data = birds.cs)

# life history
m7 <- lm(sc.mass ~ MigStatus, data = birds.cs)
m8 <- lm(sc.mass ~ Sex, data = birds.cs)

# temporal
m9 <- lm(sc.mass ~ Julian, data = birds.cs)
m10 <- lm(sc.mass ~ seconds_since_midnight, data = birds.cs)
m11 <- lm(sc.mass ~ Event, data = birds.cs)

# flock
m12 <- lm(sc.mass ~ Max_Flock_Size, data = birds.cs)

# null
m13 <- lm(sc.mass ~ 1, data = birds.cs)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# results
summary(m8)
confint(m8) # sex not significant

summary(m9)
confint(m9) # date not significant

# nothing significant

# is random effect of site helpful? no ----
# Only include sites with at least three individuals
birds.s <- birds %>% 
  group_by(Site) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# standardize data except for response
birds.s.cs <- birds.s %>%
  mutate(across(where(is.numeric) & !matches("sc.mass"), scale))

# agriculture
m1 <- lmer(sc.mass ~ PercentAg + (1|Site), REML = FALSE, data = birds.s.cs)

# vegetation
m2 <- lmer(sc.mass ~ Percent_Total_Veg + (1|Site), REML = FALSE, data = birds.s.cs)

# habitat
m3 <- lmer(sc.mass ~ Permanence + (1|Site), REML = FALSE, data = birds.s.cs)
m4 <- lmer(sc.mass ~ Percent_Exposed_Shoreline + (1|Site), REML = FALSE, 
         data = birds.s.cs)
m5 <- lmer(sc.mass ~ Dist_Closest_Wetland_m + (1|Site), REML = FALSE, data = birds.s.cs)

# weather
m6 <- lmer(sc.mass ~ SPEI + (1|Site), REML = FALSE, data = birds.s.cs)

# life history
m7 <- lmer(sc.mass ~ MigStatus + (1|Site), REML = FALSE, data = birds.s.cs)
m8 <- lmer(sc.mass ~ Sex + (1|Site), REML = FALSE, data = birds.s.cs)

# temporal
m9 <- lmer(sc.mass ~ Julian + (1|Site), REML = FALSE, data = birds.s.cs)
m10 <- lmer(sc.mass ~ seconds_since_midnight + (1|Site), REML = FALSE, data = birds.s.cs)
m11 <- lmer(sc.mass ~ Event + (1|Site), REML = FALSE, data = birds.s.cs)

# flock
m12 <- lmer(sc.mass ~ Max_Flock_Size + (1|Site), REML = FALSE, data = birds.s.cs)

# null
m13 <- lmer(sc.mass ~ 1 + (1|Site), REML = FALSE, data = birds.s.cs)


model_names <- paste0("m", 1:13)

models <- mget(model_names)

aictab(models, modnames = model_names)

# Is there enough support to include the random effect? No ----
m2 <- lmer(sc.mass ~ Sex + (1|Site), 
           data = birds.s.cs, REML = FALSE)

m1 <- lm(sc.mass ~ Sex, data = birds.s.cs)

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

m <- lm(sc.mass ~ Biomass, data = birds.sub)

summary(m)
confint(m) # no effect of biomass

m <- lm(sc.mass ~ Diversity, data = birds.sub)

summary(m)
confint(m) # no effect of diversity






