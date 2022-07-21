source("load_tables.R")
library(aod)

data_muts_severe$month_imputed <- data_muts_severe$day_imputed / 30
data_muts_severe$age_imputed_10 <- data_muts_severe$age_imputed / 10

base_regression_severity <- '`Severe Disease` ~ age_imputed + Gender.num + month_imputed + Vaccination.Status..Processed. + '

plot_base <- function(base_reg, base_vars=c("age_imputed_10", "Gender.num", "month_imputed", "Vaccination.Status..Processed.Yes"),
                      nice_vars=c("Age (10 years)", "Gender (Female)", "Time (months)", "Vaccination status")){
  estimate <- summary(base_reg)$coefficients[base_vars, 1]
  standard_error <- summary(base_reg)$coefficients[base_vars, 2]
  ci_lower <- estimate - qnorm(0.975) * standard_error
  ci_upper <- estimate + qnorm(0.975) * standard_error
  base_df <- data.frame("OR" = exp(estimate),
                        "Lower CI" = exp(ci_lower),
                        "Upper CI" = exp(ci_upper),
                        "p-value" = summary(base_reg)$coefficients[base_vars, 4],
                        "Value"= nice_vars)
  pp <- base_df %>%
    ggplot(aes(x = Value, y = OR)) +
    geom_point() +
    geom_errorbar(aes(ymin = `Lower.CI`, ymax = `Upper.CI`)) +
    ggtitle("Base variables") +
    coord_flip() +
    geom_hline(yintercept = 1, color="red") +
    scale_y_continuous(trans=log2_trans()) +
    ylab("OR") + xlab("") + theme_minimal_grid()
  return(list(base_df, pp, summary(base_reg)))
}

generate_comorbidity_plots <- function(base_regression, data) {
  comorbidities_of_interest <- c("`Other.Comorbidities`", "Diabetes",
                                   "Hypertension", "Ever.Smoke")
  comorbidities_of_interest_no_ticks <-
      c("Other", "Diabetes",
        "Hypertension", "Smoke")
  out <- vector('list', length(comorbidities_of_interest))

  comorbidity_odds_ratios <- rep(NA, length(comorbidities_of_interest))
  comorbidity_ci_lower <- rep(NA, length(comorbidities_of_interest))
  comorbidity_ci_upper <- rep(NA, length(comorbidities_of_interest))
  comorbidity_p_value <- rep(NA, length(comorbidities_of_interest))

  for (i in seq_along(comorbidities_of_interest)) {
    print(comorbidities_of_interest[i])
    out[[i]] <- glm(paste(base_regression, comorbidities_of_interest[i]),
                    data,
                    family = "binomial")
  }

  for (i in seq_along(out)) {
    regression_output <- broom::tidy(out[[i]])
    comorb_column <- grepl(gsub("`", "", comorbidities_of_interest[i]), regression_output$term)
    estimate <- regression_output[comorb_column, 2]
    standard_error <- regression_output[comorb_column, 3]
    comorbidity_odds_ratios[i] <- exp(estimate)
    comorbidity_ci_lower[i] <- exp(estimate - qnorm(0.975) * standard_error)
    comorbidity_ci_upper[i] <- exp(estimate + qnorm(0.975) * standard_error)
    comorbidity_p_value[i] <- regression_output[comorb_column, 5]
  }

  comorbidity_odds_ratios <- unlist(comorbidity_odds_ratios)
  comorbidity_ci_lower <- unlist(comorbidity_ci_lower)
  comorbidity_ci_upper <- unlist(comorbidity_ci_upper)
  comorbidity_p_value <- unlist(comorbidity_p_value)

  comorbidity_odds_ratios_df <- data.frame("OR" = comorbidity_odds_ratios,
                                           "Lower CI" = comorbidity_ci_lower,
                                           "Upper CI" = comorbidity_ci_upper,
                                           "p-value"  = comorbidity_p_value,
                                           "Comorbidity"= comorbidities_of_interest_no_ticks)

  comorbidity_odds_ratios_df$`Lower.CI`[comorbidity_odds_ratios_df$`Lower.CI`==0] <- NA
  comorbidity_odds_ratios_df$OR[comorbidity_odds_ratios_df$`Upper.CI`==Inf] <- NA
  comorbidity_odds_ratios_df$`Upper.CI`[comorbidity_odds_ratios_df$`Upper.CI`==Inf] <- NA

  pp <- comorbidity_odds_ratios_df %>%
    ggplot(aes(x = `Comorbidity`, y = `OR`)) +
    geom_point() +
    geom_errorbar(aes(ymin = `Lower.CI`, ymax = `Upper.CI`)) +
    ggtitle("Comorbidities") +
    coord_flip() +
    geom_hline(yintercept = 1, color="red") +
    scale_y_continuous(trans=log2_trans()) +
    ylab("Odds Ratio") + theme_minimal_grid()
  return(list(comorbidity_odds_ratios_df, pp))
}

generate_region_plots <- function(base_regression, data) {
  regions_of_interest <-
    c("region.processedCapital", "region.processedMuharraq",
      "region.processedNorthern", "region.processedSouthern",
      "region.processedNone")

  regions_of_interest_no_ticks <-
    c("Capital", "Muharraq", "Northern", "Southern", "Missing")

  region_odds_ratios <- rep(NA, length(regions_of_interest))
  region_ci_lower <- rep(NA, length(regions_of_interest))
  region_ci_upper <- rep(NA, length(regions_of_interest))
  region_p_value <- rep(NA, length(regions_of_interest))

  glm.out <- glm(paste(base_regression, "region.processed"),
                  data = data,
                  family = "binomial")

  regression_output <- broom::tidy(glm.out)

  for (i in seq_along(regions_of_interest)) {
    region_column <- grepl(gsub("`", "", regions_of_interest[i]), regression_output$term)
    if(sum(region_column)==0){
      estimate <- NA
      standard_error <- NA
      region_p_value[i] <- NA
    } else{
      estimate <- regression_output[region_column, 2]
      standard_error <- regression_output[region_column, 3]
      region_p_value[i] <- regression_output[region_column, 5]
    }
    region_odds_ratios[i] <- exp(estimate)
    region_ci_lower[i] <- exp(estimate - qnorm(0.975) * standard_error)
    region_ci_upper[i] <- exp(estimate + qnorm(0.975) * standard_error)

  }

  region_odds_ratios <- unlist(region_odds_ratios)
  region_ci_lower <- unlist(region_ci_lower)
  region_ci_upper <- unlist(region_ci_upper)
  region_p_value <- unlist(region_p_value)
  print(region_p_value)

  max_or <- max(region_odds_ratios)

  region_odds_ratios_df <- data.frame("OR" = region_odds_ratios,
                                      "Lower CI" = region_ci_lower,
                                      "Upper CI" = region_ci_upper,
                                      "p-value" = region_p_value,
                                      "Region"= regions_of_interest_no_ticks)

  pp <- region_odds_ratios_df %>%
    ggplot(aes(x = `Region`, y = `OR`)) +
    geom_point() +
    geom_errorbar(aes(ymin = `Lower.CI`, ymax = `Upper.CI`)) +
    ggtitle("Region") +
    coord_flip() +
    geom_hline(yintercept = 1, color="red") +
    scale_y_continuous(trans=log2_trans()) +
    ylab("Odds Ratio") + xlab("") + theme_minimal_grid()
  return(list(region_odds_ratios_df, pp, glm.out))
}

generate_nationality_plots <- function(base_regression, data) {
  nationalities_of_interest <-
    c("Nationality_processedBahraini", "Nationality_processedS_Asia",
      "Nationality_processedOther", "Nationality_processedNone")

  nationalities_of_interest_no_ticks <-
    c("Bahraini", "South Asia", "Other", "None")

  nationality_odds_ratios <- rep(NA, length(nationalities_of_interest))
  nationality_ci_lower <- rep(NA, length(nationalities_of_interest))
  nationality_ci_upper <- rep(NA, length(nationalities_of_interest))
  nationality_p_value <- rep(NA, length(nationalities_of_interest))

  glm.out <- glm(paste(base_regression, "Nationality_processed"),
                  data = data,
                  family = "binomial")
  regression_output <- broom::tidy(glm.out)

  for (i in seq_along(nationalities_of_interest)) {
    nationality_column <- grepl(gsub("`", "", nationalities_of_interest[i]), regression_output$term)
    if(sum(nationality_column)==0){
      estimate <- NA
      standard_error <- NA
      nationality_p_value[i] <- NA
    } else{
      estimate <- regression_output[nationality_column, 2]
      standard_error <- regression_output[nationality_column, 3]
      nationality_p_value[i] <-  regression_output[nationality_column, 5]
    }
    nationality_odds_ratios[i] <- exp(estimate)
    nationality_ci_lower[i] <- exp(estimate - qnorm(0.975) * standard_error)
    nationality_ci_upper[i] <- exp(estimate + qnorm(0.975) * standard_error)
  }

  nationality_odds_ratios <- unlist(nationality_odds_ratios)
  nationality_ci_lower <- unlist(nationality_ci_lower)
  nationality_ci_upper <- unlist(nationality_ci_upper)
  nationality_p_value <- unlist(nationality_p_value)

  nationality_odds_ratios_df <- data.frame("OR" = nationality_odds_ratios,
                                           "Lower CI" = nationality_ci_lower,
                                           "Upper CI" = nationality_ci_upper,
                                           "p-value" = nationality_p_value,
                                           "Nationality"= nationalities_of_interest_no_ticks)

  pp <- nationality_odds_ratios_df %>%
    ggplot(aes(x = `Nationality`, y = `OR`)) +
    geom_point() +
    geom_errorbar(aes(ymin = `Lower.CI`, ymax = `Upper.CI`)) +
    ggtitle("Nationality") +
    coord_flip() +
    geom_hline(yintercept = 1, color="red") +
    scale_y_continuous(trans=log2_trans()) +
    ylab("Odds Ratio") + xlab("") + theme_minimal_grid()
  return(list(nationality_odds_ratios_df, pp, glm.out))
}

severity_base <- plot_base(glm('`Severe Disease` ~ age_imputed_10 + Gender.num + month_imputed + Vaccination.Status..Processed.',
                               data_muts_severe,
                               family = "binomial"))

severity_comorb <- generate_comorbidity_plots(base_regression_severity, data_muts_severe)

severity_region <- generate_region_plots(base_regression_severity, data_muts_severe)

aod::wald.test(b = coef(severity_region[[3]]), Sigma = vcov(severity_region[[3]]), Terms=6:9)

severity_nationality <- generate_nationality_plots(base_regression_severity, data_muts_severe)

aod::wald.test(b = coef(severity_nationality[[3]]), Sigma = vcov(severity_nationality[[3]]), Terms=6:8)

plot_full_noage <- function(base_reg, base_vars=c("age_imputed_10",
                                            "Gender.num",
                                            "month_imputed",
                                            "Vaccination.Status..Processed.Yes",
                                            "Other.ComorbiditiesYes",
                                            "DiabetesYes",
                                            "HypertensionYes",
                                            "Ever.SmokeYes",
                                            "region.processedMuharraq", "region.processedNone", "region.processedNorthern", "region.processedSouthern",
                                            "Nationality_processedNone", "Nationality_processedOther", "Nationality_processedS_Asia"),
                      nice_vars=c("Age (10 years)",
                                  "Gender (Female)",
                                  "Time (months)",
                                  "Vaccination Status",
                                  "Misc. comorb.",
                                  "Diabetes",
                                  "Hypertension",
                                  "Smoke",
                                  "Region (Muharraq)", "Region (Missing)", "Region (Northern)", "Region (Southern)",
                                  "Nationality (Missing)", "Nationality (Other)", "Nationality (S. Asia)")){
  estimate <- summary(base_reg)$coefficients[base_vars, 1]
  standard_error <- summary(base_reg)$coefficients[base_vars, 2]
  ci_lower <- estimate - qnorm(0.975) * standard_error
  ci_upper <- estimate + qnorm(0.975) * standard_error
  base_df <- data.frame("OR" = exp(estimate),
                        "Lower CI" = exp(ci_lower),
                        "Upper CI" = exp(ci_upper),
                        "p-value" = summary(base_reg)$coefficients[base_vars, 4],
                        "Value"= nice_vars)
  base_df$Value <- factor(base_df$Value, levels=rev(nice_vars))
  pp <- base_df %>%
    ggplot(aes(x = Value, y = OR)) +
    geom_point() +
    geom_errorbar(aes(ymin = `Lower.CI`, ymax = `Upper.CI`)) +
    coord_flip() +
    geom_hline(yintercept = 1, color="red") +
    scale_y_continuous(trans=log2_trans()) +
    ylab("OR") + xlab("") + theme_minimal_grid()
  return(list(base_df, pp))
}


fit.full.noage <- glm('`Severe Disease` ~ age_imputed_10 + Gender.num+ month_imputed + Vaccination.Status..Processed. +
`Other.Comorbidities`+ Diabetes+
    Hypertension+ Ever.Smoke + region.processed + Nationality_processed',
                      data_muts_severe,
                      family = "binomial")

severity_full_noage <- plot_full_noage(fit.full.noage)

