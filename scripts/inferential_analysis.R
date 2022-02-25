devtools::install_github("glmmTMB/glmmTMB/glmmTMB") # Load the developer version of the package to 
# run 0 and 1 inflated beta regression

library(ggplot2)
library(corrplot)
source("http://www.sthda.com/upload/rquery_cormat.r") # to use corrmat
library(glmmTMB)
library(devtools)

library(readxl)
library(tidyverse)
library(dplyr)
library(epiDisplay) # for frequency graphs (tab1)

library(broom.mixed)
library(dotwhisker)
library(sjPlot) # to run plot_model

# 1 - Prepare data
##################################################

df <- read_csv("consumables_df.csv", col_names = TRUE)

# Sort month
month_order = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September',
                'October', 'November', 'December')

# Keep only LMIS data
df <- df %>% 
  filter(stringr::str_detect(data_source, 'lmis')) %>% 
  arrange(fct_relevel(month, month_order))
df$dist_todh <- df$dist_todh/1000 # convert meters to kilometers
df$dist_torms <- df$dist_torms/1000 # convert meters to kilometers
df$lndist_todh <- log(df$dist_todh + 1)
df$lndist_torms <- log(df$dist_torms + 1)

# Convert closing balance figure to opening balance figure
# for each fac, item and month, opening bal is the closing balance of the previous month
items_list <- unique(df$item_code)
facs_list <- unique(df$fac_name)
df$opening_bal = NaN
for (fac in facs_list){
  for (item in items_list){
    for (i in seq(2, 12, by=1)){
      cond_currmonth = df$item_code == item & df$fac_name == fac & df$month == month_order[i]
      cond_prevmonth = df$item_code == item & df$fac_name == fac & df$month == month_order[i-1]
      df$opening_bal[cond_currmonth] = df$closing_bal[cond_prevmonth]
    }
  }
}

# Add number of days in a month column for weights entry in the glmer function
df$mthdays = NaN
cond_31 = df$month %in% c("January", "March", "May", "July", "August", "October", "December")
cond_30 = df$month %in% c("April", "June", "September", "November")
cond_28 = df$month %in% "February"
df$mthdays[cond_28] = 28
df$mthdays[cond_30] = 30
df$mthdays[cond_31] = 31

# Keep only relevant columns for regression
cols <- c('district', 'fac_type_tlo', 'fac_name', 'category', 'item_code', 'month', 
          'available_prop', 'mthdays',
          'amc', 'dispensed', 'received', 'opening_bal',
          'dist_todh', 'dist_torms', 'drivetime_todh', 'drivetime_torms',
          'lndist_todh', 'lndist_torms')
regdf <- df[cols]

# Broad summaries to check for data inconsistencies
summary_by_mth_datasource <- df %>%
  group_by(month, data_source) %>%
  summarize(mean_size = mean(available_prop, na.rm = TRUE)) %>%
  arrange(fct_relevel(month, month_order))

summary_by_datasource <- df %>%
  group_by(data_source) %>%
  summarize(mean_size = mean(available_prop, na.rm = TRUE))

# 2 - Descriptive analysis
################################################
# Descriptive graphs
#********************
# Monthly availability
ggplot(aes(x = month, y = available_prop), data = regdf) + stat_summary(fun = "mean", geom = "bar")

# Availability by level of care
ggplot(aes(x = fac_type_tlo, y = available_prop), data = regdf) + stat_summary(fun = "mean", geom = "bar")

# Availability by program
ggplot(aes(x = category, y = available_prop), data = regdf) + stat_summary(fun = "mean", geom = "bar")

# Availability by district
ggplot(aes(x = district, y = available_prop), data = regdf) + stat_summary(fun = "mean", geom = "bar")

# Availability by distance from DHO
plot(regdf$dist_todh, regdf$available_prop, main = "Distance from DHO and Consumable availability",
     xlab = "Distance from District Health Office (kilometers)", ylab = "Availability (%)",
     pch = 19, frame = FALSE)
abline(lm(regdf$available_prop ~ regdf$dist_todh, data = mtcars), col = "blue")

# Availability by datasource

# Availability by distance from DHO
plot(regdf$dist_torms, regdf$available_prop, main = "Distance from RMS and Consumable availability",
     xlab = "Distance from Regional Medical Store (kilometers)", ylab = "Availability (%)",
     pch = 19, frame = FALSE)
abline(lm(regdf$available_prop ~ regdf$dist_torms, data = mtcars), col = "blue")

# 3 - Regression analysis
##################################################

# --- 3.1 Check for collinearity --- #
#------------------------------------#
numcols <- c('available_prop', 
             'amc', 'dispensed', 'received', 'opening_bal',
             'dist_todh', 'dist_torms', 'drivetime_todh', 'drivetime_torms')
cormat<-rquery.cormat(regdf[numcols], graphType="heatmap")
# highly correlated pairs - distance to DHO and distance to RMS; dispensed and amc; 
# drivetime to dist; ~dispensed and opening balance

# --- 3.2 Define a function to generate residual plots for models --- #
#---------------------------------------------------------------------#
gen_res_plots <- function(model, name){
  # One figure in row 1 and two figures in row 2
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  
  pred = predict(model, type = 'response')
  res = regdf$available_prop - pred
  
  # Create a density plot
  plot(density(res))
  
  # Create Q-Q plot for residuals
  qqnorm(res)
  qqline(res) #add a straight diagonal line to the plot
  
  # Create a scatter plot showing the relationship between residuals and fitted values
  plot(pred, res, 
       xlab="Fitted values", ylab="Residuals", 
       main="") 
  abline(0, 0)
  
  dev.copy(jpeg,filename=name);
  dev.off ()
  
}

# --- 3.3 Test different models --- #
#-----------------------------------#
# I. Basic Logit model with interactions
# (The results are the same as a fractional model)
model_logit <- logit(available_prop ~ amc + fac_type_tlo + category * month + lndist_todh * district + category * lndist_torms, data = regdf)
summary(model_logit)
gen_res_plots(model_logit, "model_logit_diagnosis.jpeg")

# II. Multilevel binomial model with facility random effects
mlt_binom_fac <- glmmTMB(
  available_prop ~ closing_bal + amc + dispensed
  + fac_type_tlo + district
  + category + month + lndist_todh + lndist_torms
  + (1|fac_name),
  data = regdf,
  family = binomial,
  REML = TRUE
)
gen_res_plots(mlt_binom_fac, "mlt_binom_fac_diagnosis.jpeg")

# III. Multilevel binomial  model with item and facility random effects
mlt_binom_facitem <- glmmTMB(
  available_prop ~ closing_bal + amc + dispensed
  + fac_type_tlo + district
  + category + month + lndist_todh + lndist_torms
  + (1|fac_name/item_code),
  data = regdf,
  family = binomial,
  REML = TRUE
)
gen_res_plots(mlt_binom_facitem, "mlt_binom_facitem_diagnosis.jpeg")

# IV. Multilevel beta model with item and facility level random effects
regdf$available_prop[regdf$available_prop == 0] = 0.00000000000000000001
mlt_beta1 <- glmmTMB(available_prop ~ closing_bal + amc + dispensed
                                    + fac_type_tlo + district
                                    + category + month + lndist_todh + lndist_torms
                                    + (1|fac_name/item_code), 
                                    data = regdf, ziformula=~1,family=list(family="beta",link="logit"))
gen_res_plots(mlt_beta1, "mlt_beta1_diagnosis.jpeg")

# --- 3.4 Plot results from the chosen model --- #
#------------------------------------------------#
# Multilevel binomial  model with item and facility random effects
# 1. Full model representation
plot_model(mixedeff_binom_item)

# 2. Plot marginal effects
plot_model(mlt_binom_facitem, type = "pred", terms = "category") # Program
plot_model(mlt_binom_facitem, type = "pred", terms = "district") # District 
plot_model(mlt_binom_facitem, type = "pred", terms = "month") # Month
plot_model(mlt_binom_facitem, type = "pred", terms = "fac_type_tlo") # Level of care
plot_model(model_lm, type = "pred", terms = "lndist_todh") # Distance from DHO
