# Laad benodigde pakketten
library(haven)
library(dplyr)
library(lubridate)
library(sandwich)
library(lmtest)
library(stargazer)

# Load data
data <- read_dta("temp5a.dta")

# Exacte omzetting van Stata-code:
data <- data %>%
  mutate(
    date = as.Date(date, origin = "1960-01-01"),  # Stata datum
    dow = wday(date) - 1,                           # zondag = 0
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
    t2 = t^2,
    t3 = t^3,
    t4 = t^4,
    t5 = t^5,
    t6 = t^6,
    t7 = t^7,
    t8 = t^8,
    t9 = t^9,
    weekday = ifelse(dow %in% c(0,6), 0, 1),
    hourweekday = hour * weekday,
    cvar = as.integer(as.numeric(date)/35),
    timevariable = as.numeric(date) * 24 + hour
  ) %>% 
  filter(date >= as.Date("1986-01-01")) %>%
  # Donut trimming: verwijder de eerste 14 dagen na de cutoff
  filter(!(date > as.Date("1989-11-20") & date < as.Date("1989-12-12"))) %>%
  arrange(timevariable)

# Maak de benodigde lagvariabelen (voor Wtmp, Wwsp en Wrh)
data <- data %>%
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

# Zorg dat de data een data.frame is (geen speciale panelclass)
data <- as.data.frame(data)

# Lijst van pollutants
pollutants <- c("mco", "mno2", "mo3", "mnox", "mso2")

# Loop voor regressies per pollutant met glm
for (x in pollutants) {
  
  # Maak de log-transform van de pollutantvariabele
  data[[paste0("log", x)]] <- log(data[[x]])
  
  # Stel de regressieformule samen
  formula <- as.formula(paste0(
    "log", x, " ~ hnc + t + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + ",
    "Wtmp + Wtmp2 + Wtmp3 + Wtmp4 + lag_Wtmp + lag_Wtmp2 + lag_Wtmp3 + lag_Wtmp4 + ",
    "Wwsp + Wwsp2 + Wwsp3 + Wwsp4 + lag_Wwsp + lag_Wwsp2 + lag_Wwsp3 + lag_Wwsp4 + ",
    "Wrh + Wrh2 + Wrh3 + Wrh4 + lag_Wrh + lag_Wrh2 + lag_Wrh3 + lag_Wrh4 + ",
    "factor(month) + factor(dow) + factor(hour) + factor(hourweekday)"
  ))
  
  # Voer de regressie uit met glm (lineaire regressie, family = gaussian)
  model <- glm(formula, data = data, family = gaussian())
  
  # Omdat glm() automatisch NA's verwijdert, halen we het gebruikte model.frame op
  model_data <- model.frame(model)
  # Bepaal de corresponderende clusterwaarden (cvar) via de rijnamen
  idx <- as.numeric(rownames(model_data))
  clusters <- data$cvar[idx]
  
  # Bereken cluster-robuste standaardfouten met vcovCL() uit sandwich (type = "HC0")
  vcov_cluster <- vcovCL(model, cluster = clusters, type = "HC0")
  summary_model <- coeftest(model, vcov. = vcov_cluster)
  
  cat("Resultaten voor pollutant:", x, "\n")
  print(summary_model)
  cat("---------------------------------------------\n")
}

###### Overzichtstabel ######

create_hnc_table <- function(data, pollutants) {
  results <- data.frame(Pollutant = character(), Coefficient = numeric(), `Standard Error` = numeric(), stringsAsFactors = FALSE)
  
  for (x in pollutants) {
    formula <- as.formula(paste0(
      "log", x, " ~ hnc + t + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + ",
      "Wtmp + Wtmp2 + Wtmp3 + Wtmp4 + lag_Wtmp + lag_Wtmp2 + lag_Wtmp3 + lag_Wtmp4 + ",
      "Wwsp + Wwsp2 + Wwsp3 + Wwsp4 + lag_Wwsp + lag_Wwsp2 + lag_Wwsp3 + lag_Wwsp4 + ",
      "Wrh + Wrh2 + Wrh3 + Wrh4 + lag_Wrh + lag_Wrh2 + lag_Wrh3 + lag_Wrh4 + ",
      "factor(month) + factor(dow) + factor(hour) + factor(hourweekday)"
    ))
    
    model <- glm(formula, data = data, family = gaussian())
    model_data <- model.frame(model)
    idx <- as.numeric(rownames(model_data))
    clusters <- data$cvar[idx]
    vcov_cluster <- vcovCL(model, cluster = clusters, type = "HC0")
    summary_model <- coeftest(model, vcov. = vcov_cluster)
    
    # Extraheer coëfficiënt en standaardfout voor de variabele 'hnc'
    coef_hnc <- summary_model["hnc", "Estimate"]
    se_hnc   <- summary_model["hnc", "Std. Error"]
    
    results <- rbind(results, data.frame(Pollutant = x, Coefficient = round(coef_hnc, 2), 
                                         `Standard Error` = round(se_hnc, 2)))
  }
  
  print(stargazer(results, type = "text", summary = FALSE))
  return(results)
}

hnc_results <- create_hnc_table(data, pollutants)

