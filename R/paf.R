#' Estimate attributable fraction in the exposed (AFE) and population attributable fraction (PAF)
#'
#' @description Computes the attributable fraction in the exposed (AFE) and the population attributable fraction (PAF)
#' with 95% Confidence Intervals (95% CIs) for a binary exposure (coded 0/1) and either:
#' * a binary outcome (prevalence / risk) when `y_t` is `NULL`, or
#' * a time-to-event outcome (incidence) using a Cox proportional hazards model when `y_t` is provided.
#'
#' For the binary-outcome case, risks are estimated as group-wise means and the relative risk (RR) is
#' `risk_exposed / risk_unexposed`. Wald 95% confidence intervals are calculated on the log(RR) scale
#' using a delta-method standard error, then transformed to RR, AFE, and PAF.
#'
#' For the time-to-event case, a Cox model is fitted and the hazard ratio (HR) for exposure is used as the
#' effect measure. Wald 95% confidence intervals are taken from the Cox model, then transformed to AFE and PAF.
#' Crude incidence rates (per 1,000 person-years) are also returned by exposure group (assumes time is in years).
#'
#' @details
#' **Coding assumptions:**
#' 1. `x` should be coded 0/1 (unexposed/exposed).
#' 2. For the binary-outcome case, `y` should be 0/1.
#' 3. For the time-to-event case, `y` is the event indicator (0/1) and `y_t`
#'         is follow-up time (in years if you want rates per 1,000 person-years).
#'
#' Model note (time-to-event): the Cox model is not adjusted by default. User can provide covariae string.
#'
#' @return
#' A data frame containing AFE and PAF (as proportions and percentages) with 95% CIs and p-values,
#' plus supporting statistics (RR/HR, group risks or rates, counts, and exposure prevalence).
#'
#' @param d A data frame containing exposure, outcome, and (optionally) follow-up time and covariates.
#' @param x Character string. Column name of the exposure variable (expected 0/1).
#' @param y Character string. Column name of the outcome variable (0/1; event indicator for survival).
#' @param y_t Character string. Optional column name of follow-up time for time-to-event analyses.
#'   Ideally this is in years, as crude incidence rates (per 1,000 person-years) are also returned.
#'   If not provided then treats \code{y} as a binary outcome and estimates RR-based AFE/PAF.
#'        \code{default=NULL} (character)
#' @param z A string. The covariate formula (e.g., "+age+sex") for variables found in `d`. Only used for time-to-event. 
#'   If not provided then the survival model runs without adjustment.
#'        \code{default=""} (character)
#' @param verbose Logical. Be verbose,
#'        \code{default=FALSE}
#'
#' @examples
#' 
#' # Binary exposure
#' example_data$hypertension <- dplyr::if_else(example_data$sbp>=140, 1, 0)
#' 
#' # Binary outcome (prevalence / risk)
#' result_1 <- paf(
#'   d = example_data,
#'   x = "hypertension",
#'   y = "event"
#' )
#' cat(sprintf("PAF: %.2f%% (95%% CI: %.2f%% to %.2f%%)\n",
#'             result_1$paf_percent,
#'             result_1$paf_ci_lower_percent,
#'             result_1$paf_ci_upper_percent))
#' 
#' # Time-to-event outcome (incidence)
#' result_2 <- paf(
#'   d = example_data,
#'   x = "hypertension",
#'   y = "event",
#'   y_t = "time",
#'   z = "age+sex"
#' )
#' cat(sprintf("PAF: %.2f%% (95%% CI: %.2f%% to %.2f%%)\n",
#'             result_2$paf_percent,
#'             result_2$paf_ci_lower_percent,
#'             result_2$paf_ci_upper_percent))
#'
#' @author Luke Pilling
#'
#' @export
#'
paf <- function(
  d, 
  x, 
  y, 
  y_t = NULL,
  z = "",
  verbose = FALSE
)  {
  
  v <- packageVersion("yodr")
  if (verbose)  cli::cli_alert_info("yodr v{v}")
  if (verbose)  cli::cli_alert("Estimating attributable fraction in the exposed (AFE) and population attributable fraction (PAF) with 95% CIs")
  start_time <- Sys.time()
  
  # check inputs
  if (class(x) != "character")  stop("x needs to be a string or character vector")
  if (class(y) != "character")  stop("y needs to be a string or character vector")
  if (!is.null(y))  if (class(y) != "character")  stop("y_t needs to be a string or character vector")
  if (class(z) != "character")  stop("z needs to be a string")
  if (! any(class(d) %in% c("data.frame","tbl","tbl_df")))  stop("d needs to be a data.frame or tibble. Best to explicitly provide inputs using `phewas(x=\"BMI\", y=\"diabetes\", z=\"sex\", d=data)` or `data |> phewas(x=\"BMI\", y=\"diabetes\", z=\"sex\")`")
  
  # check variables are all in d
  if (! x  %in% colnames(d))  cat("!! Exposure variable `x` not in the provided data frame\n")
  if (! y  %in% colnames(d))  cat("!! Outcome variable `y` not in the provided data frame\n")

  z_vars = stringr::str_replace_all(z, " |as.factor\\(|haven::|\\)", "") |> 
           stringr::str_split_1("\\+|\\*") |> 
           purrr::keep(\(x) stringr::str_length(x)>0)
  if (any(! z_vars %in% colnames(d)))  {
    cat("!! Not all covariate variables are in the provided data\n")
    missing <- z_vars[!z_vars %in% colnames(d)]
    int(missing)
    stop()
  }
  if (verbose)  cat("All x, y and z variables are in d\n")
  
  # check z formula starts with a "+" - if not, add one (unless string is empty)
  if (z != "")  {
    z = stringr::str_replace_all(z, " ", "")
    if (stringr::str_sub(z,start=1,end=1) != "+")  z = paste0("+",z)
  }
  
  # for binary outcome (prevalence)
  if (is.null(y_t)) {
    
    if (verbose)  cat("Estimating Relative Risk (RR)\n")
    
    # exposed group
    p_exposed_df <- d |>
      dplyr::filter(!!rlang::sym(x) == 1) |>
      dplyr::summarise(
        prev = mean(!!rlang::sym(y), na.rm = TRUE),
        n = dplyr::n(),
        n_cases = sum(!!rlang::sym(y), na.rm = TRUE)
      )
    
    # unexposed group
    p_unexposed_df <- d |>
      dplyr::filter(!!rlang::sym(x) == 0) |>
      dplyr::summarise(
        prev = mean(!!rlang::sym(y), na.rm = TRUE),
        n = dplyr::n(),
        n_cases = sum(!!rlang::sym(y), na.rm = TRUE)
      )
    
    p_e <- p_exposed_df$prev
    p_u <- p_unexposed_df$prev
    n_e <- p_exposed_df$n
    n_u <- p_unexposed_df$n
    
    # total sample size and cases
    n_total <- n_e + n_u
    prop_exposed <- n_e / n_total
    
    # calculate relative risk
    rr <- p_e / p_u
    
    # se of log(rr) using delta method
    se_log_rr <- sqrt((1 - p_e) / (n_e * p_e) + (1 - p_u) / (n_u * p_u))
    
    # 95% ci for rr
    rr_ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
    rr_ci_upper <- exp(log(rr) + 1.96 * se_log_rr)
    
    if (verbose)  cat("Estimating AFE and PAF\n")
    
    # attributable fraction in exposed: afe = (rr - 1) / rr
    afe <- (rr - 1) / rr
    afe_ci_lower <- (rr_ci_lower - 1) / rr_ci_lower
    afe_ci_upper <- (rr_ci_upper - 1) / rr_ci_upper
    
    # se of afe using delta method: se(afe) = se(log(rr)) / rr
    se_afe <- se_log_rr / rr
    z_afe <- afe / se_afe
    p_value_afe <- 2 * (1 - pnorm(abs(z_afe)))
    
    # population attributable fraction: paf = p_e * (rr - 1) / (p_e * (rr - 1) + 1)
    paf <- prop_exposed * (rr - 1) / (prop_exposed * (rr - 1) + 1)
    paf_ci_lower <- prop_exposed * (rr_ci_lower - 1) / (prop_exposed * (rr_ci_lower - 1) + 1)
    paf_ci_upper <- prop_exposed * (rr_ci_upper - 1) / (prop_exposed * (rr_ci_upper - 1) + 1)
    
    # se of paf using delta method
    se_paf <- se_log_rr * prop_exposed / ((prop_exposed * (rr - 1) + 1) * rr)
    z_paf <- paf / se_paf
    p_value_paf <- 2 * (1 - pnorm(abs(z_paf)))
    
    # number of excess cases in exposed
    excess_cases_exposed <- p_exposed_df$n_cases - (n_e * p_u)
    
    res <- data.frame(
      # attributable fraction in exposed
      afe = afe,
      afe_ci_lower = afe_ci_lower,
      afe_ci_upper = afe_ci_upper,
      afe_percent = afe * 100,
      afe_ci_lower_percent = afe_ci_lower * 100,
      afe_ci_upper_percent = afe_ci_upper * 100,
      afe_p_value = p_value_afe,
      
      # population attributable fraction
      paf = paf,
      paf_ci_lower = paf_ci_lower,
      paf_ci_upper = paf_ci_upper,
      paf_percent = paf * 100,
      paf_ci_lower_percent = paf_ci_lower * 100,
      paf_ci_upper_percent = paf_ci_upper * 100,
      paf_p_value = p_value_paf,
      
      # supporting statistics
      relative_risk = rr,
      rr_ci_lower = rr_ci_lower,
      rr_ci_upper = rr_ci_upper,
      risk_exposed = p_e,
      risk_unexposed = p_u,
      excess_risk = p_e - p_u,
      excess_cases_exposed = excess_cases_exposed,
      n_exposed = n_e,
      n_unexposed = n_u,
      n_cases_exposed = p_exposed_df$n_cases,
      n_cases_unexposed = p_unexposed_df$n_cases,
      prop_exposed = prop_exposed,
      x = x,
      y = y
    )
  }
  
  # for time-to-event outcome (incidence)
  if (!is.null(y_t)) {
    
    if (verbose)  cat("Fitting Cox model\n")
    
    if (! y_t  %in% colnames(d))  cat("!! Outcome time variable `y_t` not in the provided data frame\n")
    
    # fit cox model
    surv_obj <- survival::Surv(
      time = d[[y_t]], 
      event = d[[y]]
    )
    
    cox_model <- survival::coxph(
      as.formula(paste0("surv_obj ~ factor(", x, ")",z)),
      data = d
    )
    
    # hazard ratio, ci, and p-value
    hr <- exp(coef(cox_model)[1])
    hr_ci <- exp(confint(cox_model)[1, ])
    hr_p_value <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
    
    if (verbose)  cat("Transforming estimates to AFE and PAF\n")
    
    # proportion exposed
    n_exposed <- sum(d[[x]] == 1)
    n_total <- nrow(d)
    prop_exposed <- n_exposed / n_total
    
    # attributable fraction in exposed: afe = (hr - 1) / hr
    afe <- (hr - 1) / hr
    afe_ci_lower <- (hr_ci[1] - 1) / hr_ci[1]
    afe_ci_upper <- (hr_ci[2] - 1) / hr_ci[2]
    
    # se of log(hr)
    se_log_hr <- sqrt(summary(cox_model)$coefficients[1, "se(coef)"]^2)
    
    # se of afe using delta method
    se_afe <- se_log_hr / hr
    z_afe <- afe / se_afe
    afe_p_value <- 2 * (1 - pnorm(abs(z_afe)))
    
    # population attributable fraction: paf = p_e * (hr - 1) / (p_e * (hr - 1) + 1)
    paf <- prop_exposed * (hr - 1) / (prop_exposed * (hr - 1) + 1)
    paf_ci_lower <- prop_exposed * (hr_ci[1] - 1) / (prop_exposed * (hr_ci[1] - 1) + 1)
    paf_ci_upper <- prop_exposed * (hr_ci[2] - 1) / (prop_exposed * (hr_ci[2] - 1) + 1)
    
    # se of paf using delta method
    se_paf <- se_log_hr * prop_exposed / ((prop_exposed * (hr - 1) + 1) * hr)
    z_paf <- paf / se_paf
    paf_p_value <- 2 * (1 - pnorm(abs(z_paf)))
    
    # incidence rates per 1000 person-years
    rates <- d |>
      dplyr::group_by(!!rlang::sym(x)) |>
      dplyr::summarise(
        n_events = sum(!!rlang::sym(y)),
        person_years = sum(!!rlang::sym(y_t)),
        rate_per_1000py = (sum(!!rlang::sym(y)) / sum(!!rlang::sym(y_t))) * 1000,
        .groups = "drop"
      )
    
    rate_exposed <- rates |>
      dplyr::filter(!!rlang::sym(x) == 1) |>
      dplyr::pull(rate_per_1000py)
    
    rate_unexposed <- rates |>
      dplyr::filter(!!rlang::sym(x) == 0) |>
      dplyr::pull(rate_per_1000py)
    
    res <- data.frame(
      # attributable fraction in exposed (main result)
      afe = afe,
      afe_ci_lower = afe_ci_lower,
      afe_ci_upper = afe_ci_upper,
      afe_percent = afe * 100,
      afe_ci_lower_percent = afe_ci_lower * 100,
      afe_ci_upper_percent = afe_ci_upper * 100,
      afe_p_value = afe_p_value,
      
      # population attributable fraction
      paf = paf,
      paf_ci_lower = paf_ci_lower,
      paf_ci_upper = paf_ci_upper,
      paf_percent = paf * 100,
      paf_ci_lower_percent = paf_ci_lower * 100,
      paf_ci_upper_percent = paf_ci_upper * 100,
      paf_p_value = paf_p_value,
      
      # supporting statistics
      hazard_ratio = hr,
      hr_ci_lower = hr_ci[1],
      hr_ci_upper = hr_ci[2],
      hr_p_value = hr_p_value,
      rate_exposed_per_1000py = rate_exposed,
      rate_unexposed_per_1000py = rate_unexposed,
      prop_exposed = prop_exposed,
      n_exposed = n_exposed,
      #incidence_summary = rates
      x = x,
      y = y,
      y_t = y_t,
      z = z
    )
    rownames(res) <- NULL
  }
  
  if (verbose)  cli::cli_alert_success(c("Finished. Time taken: ", "{prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), start_time, units=\"secs\")))}."))
  
  return(res)
}
