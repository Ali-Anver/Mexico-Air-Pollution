# Load required packages
library(haven)
library(dplyr)
library(lubridate)
library(sandwich)
library(lmtest)
library(stargazer)

# Function to compute the exact critical value following Armstrong & Kolesár:
# Find c such that:
#    Φ(c - t_star) + Φ(c + t_star) - 1 = 1 - alpha.
# Dynamically determine the upper bound, as the solution might be larger for large t_star.
critical_value <- function(t_star, alpha) {
  lower_bound <- 0
  upper_bound <- max(10, t_star + qnorm(1 - alpha) + 1)
  objective <- function(c) {
    pnorm(c - t_star) + pnorm(c + t_star) - 1 - (1 - alpha)
  }
  uniroot(objective, lower = lower_bound, upper = upper_bound)$root
}

# Set alpha and smoothness parameter
alpha <- 0.05
M <- 1   # Smoothness constant; choose based on domain knowledge

# Set a local window (in days)
h_window <- 2920

# Load data
data <- read_dta("temp5a.dta")

# Add "hour" if missing (default to 0)
if(!"hour" %in% names(data)) {
  data$hour <- 0
}

# Transform data, including original variables
data <- data %>%
  mutate(
    date = as.Date(date, origin = "1960-01-01"),  # Stata date
    dow = wday(date) - 1,                           # Sunday = 0
    month = month(date),
    year = year(date),
    Wtmp2 = Wtmp^2,
    Wtmp3 = Wtmp^3,
    Wtmp4 = Wtmp^4,
    Wwsp2 = Wwsp^2,
    Wwsp3 = Wwsp^3,
    Wwsp4 = Wwsp^4,
    Wrh2 = Wrh^2,
    Wrh3 = Wrh^3,
    Wrh4 = Wrh^4,
    hnc = ifelse(date >= as.Date("1989-11-20"), 1, 0),
    t = as.numeric(date - as.Date("1986-01-01")),
    weekday = ifelse(dow %in% c(0,6), 0, 1),
    hourweekday = hour * ifelse(dow %in% c(0,6), 0, 1),
    cvar = as.integer(as.numeric(date) / 35),
    timevariable = as.numeric(date) * 24 + hour
  ) %>% 
  filter(date >= as.Date("1986-01-01")) %>%
  filter(!(date > as.Date("1989-11-20") & date < as.Date("1989-11-28"))) %>%
  arrange(t)

# Determine cutoff in t (days)
t_cutoff <- as.numeric(as.Date("1989-11-20") - as.Date("1986-01-01"))
cat("Cutoff t =", t_cutoff, "\n")

# Restrict data to a local window around cutoff (± h_window days)
data_local <- data %>% filter(abs(t - t_cutoff) <= h_window)

# Scale t-variable so support is [–1, 1]
data_local <- data_local %>% mutate(t_scaled = (t - t_cutoff) / h_window)

# Create lagged variables in local dataset
data_local <- data_local %>%
  mutate(
    lag_Wtmp    = lag(Wtmp, 1),
    lag_Wtmp2   = lag(Wtmp2, 1),
    lag_Wtmp3   = lag(Wtmp3, 1),
    lag_Wtmp4   = lag(Wtmp4, 1),
    lag_Wwsp    = lag(Wwsp, 1),
    lag_Wwsp2   = lag(Wwsp2, 1),
    lag_Wwsp3   = lag(Wwsp3, 1),
    lag_Wwsp4   = lag(Wwsp4, 1),
    lag_Wrh     = lag(Wrh, 1),
    lag_Wrh2    = lag(Wrh2, 1),
    lag_Wrh3    = lag(Wrh3, 1),
    lag_Wrh4    = lag(Wrh4, 1)
  )

data_local <- as.data.frame(data_local)

# Check existence of hourweekday in local data; recalculate if missing
if(!"hourweekday" %in% names(data_local)) {
  data_local <- data_local %>% mutate(hourweekday = hour * ifelse(dow %in% c(0,6), 0, 1))
}

# List of pollutants
pollutants <- c("mco", "mno2", "mo3", "mnox", "mso2")

###############################
# Loop: GLM and exact bias-aware CIs (with scaled t)
###############################
for (x in pollutants) {
  
  # Add log transformation
  data_local[[paste0("log", x)]] <- log(data_local[[x]])
  
  # Construct regression formula with scaled t-variable and polynomial terms up to order 9
  formula <- as.formula(paste0(
    "log", x, " ~ hnc + t_scaled + I(t_scaled^2) + I(t_scaled^3) + I(t_scaled^4) + I(t_scaled^5) + ",
    "I(t_scaled^6) + I(t_scaled^7) + I(t_scaled^8) + I(t_scaled^9) + ",
    "Wtmp + Wtmp2 + Wtmp3 + Wtmp4 + lag(Wtmp, 1) + lag(Wtmp2, 1) + lag(Wtmp3, 1) + lag(Wtmp4, 1) + ",
    "Wwsp + Wwsp2 + Wwsp3 + Wwsp4 + lag(Wwsp, 1) + lag(Wwsp2, 1) + lag(Wwsp3, 1) + lag(Wwsp4, 1) + ",
    "Wrh + Wrh2 + Wrh3 + Wrh4 + lag(Wrh, 1) + lag(Wrh2, 1) + lag(Wrh3, 1) + lag(Wrh4, 1) + ",
    "factor(month) + factor(dow) + factor(hour) + factor(hourweekday)"
  ))
  
  # Estimate model using glm (gaussian)
  model <- glm(formula, data = data_local, family = gaussian())
  
  # Extract the model.frame to identify actual observations used
  model_data <- model.frame(model)
  idx <- as.numeric(rownames(model_data))
  
  # Obtain point estimate and standard error for hnc coefficient
  tau_hat <- coef(model)["hnc"]
  se <- sqrt(vcov(model)["hnc", "hnc"])
  
  # In scaled space, bandwidth is 1, so the worst-case bias bound is:
  B <- M / 2
  t_star <- B / se
  
  # Compute exact critical value using critical_value()
  cv_exact <- critical_value(t_star, alpha)
  
  # Compute bias-aware confidence interval
  lower_exact <- tau_hat - cv_exact * se
  upper_exact <- tau_hat + cv_exact * se
  
  cat("Exact bias-aware CI for pollutant:", x, "\n")
  cat("tau_hat =", tau_hat, ", se =", se, "\n")
  cat("Bias bound B =", B, ", t_star =", t_star, "\n")
  cat("Exact Critical value =", cv_exact, "\n")
  cat("CI: [", lower_exact, ",", upper_exact, "]\n")
  cat("---------------------------------------------\n")
}

###############################
# Summary table with hnc coefficients and exact bias-aware CIs
###############################
create_hnc_table <- function(data_local, pollutants, alpha = 0.05, M = 0.1) {
  results <- data.frame(Pollutant = character(), 
                        Coefficient = numeric(), 
                        CI_Lower = numeric(), 
                        CI_Upper = numeric(), 
                        stringsAsFactors = FALSE)
  
  for (x in pollutants) {
    formula <- as.formula(paste0(
      "log", x, " ~ hnc + t_scaled + I(t_scaled^2) + I(t_scaled^3) + I(t_scaled^4) + I(t_scaled^5) + ",
      "I(t_scaled^6) + I(t_scaled^7) + I(t_scaled^8) + I(t_scaled^9) + ",
      "Wtmp + Wtmp2 + Wtmp3 + Wtmp4 + lag(Wtmp, 1) + lag(Wtmp2, 1) + lag(Wtmp3, 1) + lag(Wtmp4, 1) + ",
      "Wwsp + Wwsp2 + Wwsp3 + Wwsp4 + lag(Wwsp, 1) + lag(Wwsp2, 1) + lag(Wwsp3, 1) + lag(Wwsp4, 1) + ",
      "Wrh + Wrh2 + Wrh3 + Wrh4 + lag(Wrh, 1) + lag(Wrh2, 1) + lag(Wrh3, 1) + lag(Wrh4, 1) + ",
      "factor(month) + factor(dow) + factor(hour) + factor(hourweekday)"
    ))
    
    model <- glm(formula, data = data_local, family = gaussian())
    tau_hat <- coef(model)["hnc"]
    se <- sqrt(vcov(model)["hnc", "hnc"])
    B <- M / 2
    t_star <- B / se
    cv_exact <- critical_value(t_star, alpha)
    
    lower_exact <- tau_hat - cv_exact * se
    upper_exact <- tau_hat + cv_exact * se
    
    results <- rbind(results, data.frame(Pollutant = x, 
                                         Coefficient = round(tau_hat, 2), 
                                         CI_Lower = round(lower_exact, 2), 
                                         CI_Upper = round(upper_exact, 2)))
  }
  
  print(stargazer(results, type = "text", summary = FALSE))
  return(results)
}

hnc_results <- create_hnc_table(data_local, pollutants, alpha = 0.05, M = 0.1)