# Laad benodigde packages
library(haven)
library(dplyr)
library(lubridate)
library(rdrobust)    # voor rdbwselect en rdrobust
library(sandwich)
library(lmtest)
library(stargazer)
library(fixest)      # indien nodig voor andere analyses

# Functie om het exacte kritieke waarde te berekenen volgens Armstrong & Kolesár:
# Los de vergelijking op: Φ(c - t_star) + Φ(c + t_star) - 1 = 1 - alpha.
critical_value <- function(t_star, alpha) {
  lower_bound <- 0
  upper_bound <- max(10, t_star + qnorm(1 - alpha) + 1)
  objective <- function(c) {
    pnorm(c - t_star) + pnorm(c + t_star) - 1 - (1 - alpha)
  }
  uniroot(objective, lower = lower_bound, upper = upper_bound)$root
}

# Stel alpha en de smoothness constant in
alpha <- 0.05
M <- 1  # Smoothness constant; aan te passen indien nodig

# Kies een lokale window (in dagen) en een kleine jitter-range
h_window <- 60
jitter_range <- 0.0001

# Laad en transformeer de data
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

# Als "hour" niet aanwezig is, voeg het dan toe (standaard op 0)
if(!"hour" %in% names(data)) {
  data$hour <- 0
}

# Voeg extra variabelen toe
data <- data %>%
  mutate(
    weekday = ifelse(dow %in% c(0,6), 0, 1),
    hourweekday = hour * ifelse(dow %in% c(0,6), 0, 1),
    cvar = as.integer(as.numeric(date) / 35),
    timevariable = as.numeric(date)*24 + hour
  )

# Bepaal cutoff: aantal dagen tussen 1986-01-01 en 1989-11-20
cutoff_value <- as.numeric(as.Date("1989-11-20") - as.Date("1986-01-01"))
cat("Cutoff (in days) =", cutoff_value, "\n")

# Centreer de running variable en voeg een kleine jitter toe
data <- data %>% 
  mutate(t_centered = t - cutoff_value,
         t_centered_jitter = t_centered + runif(n(), min = -jitter_range, max = jitter_range))

# Definieer de lijst van vervuilers (uitkomstvariabelen)
pollutants <- c("mco", "mno2", "mo3", "mnox", "mso2")

# Zet de log-transformaties van de vervuiler variabelen op
data <- data %>%
  mutate(across(all_of(pollutants), ~ log(.), .names = "log{.col}"))

# Bepaal de optimale bandbreedte voor één outcome (voorbeeld)
bw_obj <- rdbwselect(y = data[["logmco"]], x = data$t_centered_jitter, c = 0, bwselect = "mserd")
bw_opt <- bw_obj$bws[1]
cat("Optimale bandbreedte (IK) =", bw_opt, "\n")

# Beperk de data tot de lokale window rond de cutoff
data_local <- data %>% 
  filter(!(date > as.Date("1989-11-20") & date < as.Date("1989-11-28")),
         abs(t_centered) <= bw_opt)

# Scale t-variable zodat de support [–1, 1] wordt
data_local <- data_local %>% mutate(t_scaled = t_centered / bw_opt)

# Maak lagged variabelen in de lokale data (op basis van de originele t-variable)
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

# Aangepaste functie om een overzichtstabel te maken met bias-aware betrouwbaarheidsintervallen
# Hierbij worden de covariaten meegenomen in de rdrobust-schatting
create_hnc_table <- function(data, pollutants, alpha = 0.05, M = 1) {
  results <- data.frame()
  
  for (x in pollutants) {
    logvar <- paste0("log", x)
    
    # Bereken de optimale bandbreedte via IK-methodologie
    bw_obj <- rdbwselect(y = data[[logvar]], 
                         x = data$t_centered_jitter, 
                         c = 0, 
                         bwselect = "mserd")
    bw_opt <- bw_obj$bws[1]
    
    # Selecteer de lokale data binnen de optimale bandbreedte, nu met filtering op t_centered_jitter
    data_sub <- data %>%
      filter(!(date > as.Date("1989-11-20") & date < as.Date("1989-11-28")),
             abs(t_centered_jitter) <= bw_opt) %>%
      mutate(t_scaled = t_centered / bw_opt)
    
    # Zorg dat er geen ontbrekende waarden zijn in de respons en running variable
    keep_rows <- complete.cases(data_sub[[logvar]], data_sub$t_centered_jitter)
    data_sub <- data_sub[keep_rows, ]
    
    # Stel eerst de covariaat-dataframe samen (inclusief vaste effecten als factorvariabelen)
    covariates <- data_sub %>% 
      mutate(
        month = as.factor(month),
        dow = as.factor(dow),
        hour = as.factor(hour),
        weekday = as.factor(weekday)
      ) %>%
      dplyr::select(Wtmp, Wtmp2, Wtmp3, Wtmp4, 
                    lag_Wtmp, lag_Wtmp2, lag_Wtmp3, lag_Wtmp4, 
                    Wwsp, Wwsp2, Wwsp3, Wwsp4, 
                    lag_Wwsp, lag_Wwsp2, lag_Wwsp3, lag_Wwsp4,
                    Wrh, Wrh2, Wrh3, Wrh4, 
                    lag_Wrh, lag_Wrh2, lag_Wrh3, lag_Wrh4,
                    month, dow, hour, weekday)
    
    # Zorg ervoor dat we dezelfde rijen gebruiken: verwijder alle observaties met NA's in covariaten
    complete_idx <- complete.cases(covariates)
    data_sub <- data_sub[complete_idx, ]
    covariates <- covariates[complete_idx, ]
    
    # Bouw de covariaatmatrix op basis van de complete covariaten
    cov_matrix <- model.matrix(~ . - 1, data = covariates)
    
    # Controleer of het aantal rijen overeenkomt
    if(nrow(data_sub) != nrow(cov_matrix)){
      stop("Aantal rijen in data_sub en cov_matrix komen niet overeen.")
    }
    
    # Voer de kernelregressie uit met de driehoekkernel (triangular) en neem de covariaten mee
    rd_obj <- rdrobust(y = data_sub[[logvar]], 
                       x = data_sub$t_centered_jitter, 
                       c = 0, 
                       kernel = "epanechnikov", 
                       h = bw_opt,
                       covs = cov_matrix)
    
    # Haal het geschatte effect en de standaardfout uit de rdrobust output
    tau_hat <- rd_obj$coef[1]
    se <- rd_obj$se[1]
    
    # Bereken de biasterm en het exacte kritieke waarde
    B <- M / 2
    t_star <- B / se
    cv_exact <- critical_value(t_star, alpha)
    
    # Sla de resultaten op in de tabel
    results <- rbind(results, data.frame(Pollutant = x, 
                                         Bandwidth = round(bw_opt, 2), 
                                         Coef = round(tau_hat, 3), 
                                         SE = round(se, 3), 
                                         CI_L = round(tau_hat - cv_exact * se, 3), 
                                         CI_U = round(tau_hat + cv_exact * se, 3)))
  }
  
  print(stargazer(results, type = "text", summary = FALSE))
  return(results)
}
# Roep de functie aan met data_local (die de lag-variabelen bevat)
hnc_results <- create_hnc_table(data_local, pollutants, alpha, M)

# Voordat je rdrobust uitvoert, inspecteer de covariaatmatrix:
qr_cov <- qr(cov_matrix)
# Vind de indices van de niet-redundante kolommen:
keep_idx <- sort(qr_cov$pivot[1:qr_cov$rank])
kept_vars <- colnames(cov_matrix)[keep_idx]

# En bepaal welke kolommen mogelijk redundant zijn:
dropped_vars <- setdiff(colnames(cov_matrix), kept_vars)
cat("Gehouden covariaten:\n")
print(kept_vars)
cat("Gedropte (redundante) covariaten:\n")
print(dropped_vars)
