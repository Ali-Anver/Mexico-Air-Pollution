# Load required packages
library(haven)
library(dplyr)
library(lubridate)
library(rdrobust)    # for rdbwselect (IK methodology)
library(sandwich)
library(lmtest)
library(stargazer)
library(fixest)      # in case you need fixest for other analyses

# Function to compute exact critical value according to Armstrong & Kolesár:
# Solve equation: Φ(c - t_star) + Φ(c + t_star) - 1 = 1 - alpha.
critical_value <- function(t_star, alpha) {
  lower_bound <- 0
  upper_bound <- max(10, t_star + qnorm(1 - alpha) + 1)
  objective <- function(c) {
    pnorm(c - t_star) + pnorm(c + t_star) - 1 - (1 - alpha)
  }
  uniroot(objective, lower = lower_bound, upper = upper_bound)$root
}

# Set alpha and smoothness constant
alpha <- 0.05
M <- 1  # Smoothness constant; adjust based on domain knowledge

# Choose a local window (in days) and a small jitter range
h_window <- 60
jitter_range <- 0.0001

# Load and transform data
data <- read_dta("temp5a.dta") %>%
  mutate(
    date = as.Date(date, origin = "1960-01-01"),
    dow = wday(date) - 1,
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
    t = as.numeric(date - as.Date("1986-01-01"))
  ) %>%
  filter(date >= as.Date("1986-01-01"))

# If "hour" is not present, add it (default to 0)
if(!"hour" %in% names(data)) {
  data$hour <- 0
}

# Add additional variables
data <- data %>%
  mutate(
    weekday = ifelse(dow %in% c(0,6), 0, 1),
    hourweekday = hour * ifelse(dow %in% c(0,6), 0, 1),
    cvar = as.integer(as.numeric(date) / 35),
    timevariable = as.numeric(date)*24 + hour
  )

# Determine cutoff: days between 1986-01-01 and 1989-11-20
cutoff_value <- as.numeric(as.Date("1989-11-20") - as.Date("1986-01-01"))
cat("Cutoff (in days) =", cutoff_value, "\n")

# Center the running variable and add small jitter
data <- data %>% 
  mutate(t_centered = t - cutoff_value,
         t_centered_jitter = t_centered + runif(n(), min = -jitter_range, max = jitter_range))

# Define the list of pollutants (outcome variables)
pollutants <- c("mco", "mno2", "mo3", "mnox", "mso2")

# Compute optimal bandwidth via IK-method (rdbwselect) for one outcome (example)
data <- data %>%
  mutate(across(all_of(pollutants), ~ log(.), .names = "log{.col}"))

bw_obj <- rdbwselect(y = data[["logmco"]], x = data$t_centered_jitter, c = 0, bwselect = "mserd")
bw_opt <- bw_obj$bws[1]
cat("Optimal bandwidth (IK) =", bw_opt, "\n")

# Restrict data to local window around cutoff
data_local <- data %>% 
  filter(!(date > as.Date("1989-11-20") & date < as.Date("1989-11-28")),
         abs(t_centered) <= bw_opt)

# Scale t-variable so its support becomes [–1, 1]
data_local <- data_local %>% mutate(t_scaled = t_centered / bw_opt)

# Create lagged variables in local data (based on original t-variable)
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

# Function to create summary table (with bias-aware CIs)
create_hnc_table <- function(data, pollutants, alpha = 0.05, M = 1) {
  results <- data.frame()
  
  for (x in pollutants) {
    logvar <- paste0("log", x)
    
    bw_obj <- rdbwselect(y = data[[logvar]], x = data$t_centered_jitter, c = 0, bwselect = "mserd")
    bw_opt <- bw_obj$bws[1]
    
    data_sub <- data %>%
      filter(!(date > as.Date("1989-11-20") & date < as.Date("1989-11-28")),
             abs(t_centered) <= bw_opt) %>%
      mutate(t_scaled = t_centered / bw_opt)
    
    model <- glm(as.formula(paste0(logvar, " ~ hnc + poly(t_scaled,9) + Wtmp + Wtmp2 + Wtmp3 + Wtmp4 + lag(Wtmp,1) + lag(Wtmp2,1) + lag(Wtmp3,1) + lag(Wtmp4,1) + Wwsp + Wwsp2 + Wwsp3 + Wwsp4 + lag(Wwsp,1) + lag(Wwsp2,1) + lag(Wwsp3,1) + lag(Wwsp4,1) + Wrh + Wrh2 + Wrh3 + Wrh4 + lag(Wrh,1) + lag(Wrh2,1) + lag(Wrh3,1) + lag(Wrh4,1) + factor(month) + factor(dow) + factor(hour) + factor(hourweekday)")),
                 data = data_sub)
    
    tau_hat <- coef(model)["hnc"]
    se <- sqrt(vcov(model)["hnc","hnc"])
    B <- M/2
    t_star <- B/se
    cv_exact <- critical_value(t_star, alpha)
    
    results <- rbind(results, data.frame(Pollutant = x, Bandwidth = round(bw_opt,2), Coef = round(tau_hat,3), SE = round(se,3), CI_L = round(tau_hat - cv_exact*se,3), CI_U = round(tau_hat + cv_exact*se,3)))
  }
  print(stargazer(results,type="text",summary=FALSE))
  return(results)
}

hnc_results <- create_hnc_table(data, pollutants, alpha, M)


