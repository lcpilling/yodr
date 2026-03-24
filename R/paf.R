#' Estimate attributable fraction in the exposed (AFE) and population attributable fraction (PAF)
#'
#' @description Computes the attributable fraction in the exposed (AFE) and the population attributable fraction (PAF)
#' with 95% Confidence Intervals (95% CIs) for a binary exposure (coded 0/1) and either:
#' * a binary outcome (prevalence / risk) when `y_t` is `NULL`, or
#' * a time-to-event outcome (incidence) using a Cox proportional hazards model when `y_t` is provided.
#'
#' Binary-outcome case:
#' * If `z == ""` then risks are estimated as group-wise means and the relative risk (RR) is
#'   `risk_exposed / risk_unexposed`. Wald 95% confidence intervals are calculated on the log(RR) scale
#'   using a delta-method standard error, then transformed to RR, AFE, and PAF.
#' * If `z != ""` (i.e., covariates provided) then a logistic regression is fitted and adjusted AFE/PAF are estimated 
#'   by standardisation (g-formula): mean predicted risks under observed exposure and counterfactual exposure assignments.
#'   If `n_boot > 0`, non-parametric bootstrap 95% CIs are returned for AFE/PAF.
#'
#' Time-to-event case:
#' A Cox model is fitted and the hazard ratio (HR) for exposure is used as the effect measure.
#' Wald 95% confidence intervals are taken from the Cox model, then transformed to AFE and PAF.
#' Crude incidence rates (per 1,000 person-years) are also returned by exposure group (assumes time is in years).
#'
#' @details
#' **Coding assumptions:**
#' 1. `x` should be coded 0/1 (unexposed/exposed).
#' 2. For the binary-outcome case, `y` should be 0/1.
#' 3. For the time-to-event case, `y` is the event indicator (0/1) and `y_t`
#'         is follow-up time (in years if you want rates per 1,000 person-years).
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
#'   If not provided then treats `y` as a binary outcome and estimates AFE/PAF.
#'        `default=NULL` (character)
#' @param z A string. Covariate formula additions (e.g., "+age+sex") for variables found in `d`.
#'   If `y_t` is provided, covariates are included in the Cox model.
#'   If `y_t` is `NULL` and `z != ""`, covariates are included in a logistic regression and AFE/PAF are
#'   estimated by standardisation (g-formula).
#'        `default=""` (character)
#' @param n_boot Integer. Number of bootstrap replicates for adjusted (logistic) binary-outcome AFE/PAF.
#'   Only used when `y_t` is `NULL` and `z != ""`. If 0 then no bootstrap CIs are computed.
#'        `default=0` (integer)
#' @param verbose Logical. Be verbose, `default=FALSE`.
#'
#' @examples
#'
#' # Binary exposure
#' example_data$hypertension <- dplyr::if_else(example_data$sbp>=140, 1, 0)
#'
#' # Binary outcome (prevalence / risk) crude
#' res1 <- paf(
#'   d = example_data,
#'   x = "hypertension",
#'   y = "event"
#' )
#'
#' # Binary outcome adjusted with bootstrap CIs
#' res1_adj <- paf(
#'   d = example_data,
#'   x = "hypertension",
#'   y = "event",
#'   z = "age+sex",
#'   n_boot = 200
#' )
#'
#' # Time-to-event outcome (incidence)
#' res2 <- paf(
#'   d = example_data,
#'   x = "hypertension",
#'   y = "event",
#'   y_t = "time",
#'   z = "age+sex"
#' )
#'
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
  n_boot = 0L,
  verbose = FALSE
)  {

  v <- packageVersion("yodr")
  if (verbose) cli::cli_alert_info("yodr v{v}")
  if (verbose) cli::cli_alert("Estimating attributable fraction in the exposed (AFE) and population attributable fraction (PAF) with 95% CIs")
  start_time <- Sys.time()

  # check inputs
  if (class(x) != "character") stop("x needs to be a string or character vector")
  if (class(y) != "character") stop("y needs to be a string or character vector")
  if (!is.null(y)) if (class(y) != "character") stop("y_t needs to be a string or character vector")
  if (class(z) != "character") stop("z needs to be a string")
  if (!is.numeric(n_boot) || length(n_boot) != 1) stop("n_boot needs to be a single numeric/integer value")
  if (n_boot < 0) stop("n_boot must be >= 0")
  n_boot <- as.integer(n_boot)

  if (!any(class(d) %in% c("data.frame", "tbl", "tbl_df"))) {
    stop("d needs to be a data.frame or tibble.")
  }

  # check variables are all in d
  if (!x %in% colnames(d)) cat("!! Exposure variable `x` not in the provided data frame\n")
  if (!y %in% colnames(d)) cat("!! Outcome variable `y` not in the provided data frame\n")

  z_vars <- stringr::str_replace_all(z, " |as.factor\\(|haven::|\\)", "") |>
    stringr::str_split_1("\\+|\\*") |>
    purrr::keep(\(x) stringr::str_length(x) > 0)

  if (any(!z_vars %in% colnames(d))) {
    cat("!! Not all covariate variables are in the provided data\n")
    missing <- z_vars[!z_vars %in% colnames(d)]
    print(missing)
    stop()
  }
  if (verbose) cat("All x, y and z variables are in d\n")

  # check z formula starts with a "+" - if not, add one (unless string is empty)
  if (z != "") {
    z <- stringr::str_replace_all(z, " ", "")
    if (stringr::str_sub(z, start = 1, end = 1) != "+") z <- paste0("+", z)
  }
  
  # subset data frame to complete cases 
  varlist <- c(x, y)
  if (z != "")  varlist <- c(varlist, z_vars)
  if (!is.null(y_t))  varlist <- c(varlist, y_t)
  d <- d |>
    dplyr::select(dplyr::all_of(c(varlist))) |>
    tidyr::drop_na()

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
  n_cases_exposed <- p_exposed_df$n_cases
  n_cases_unexposed <- p_unexposed_df$n_cases

  # total sample size and cases
  n_total <- n_e + n_u
  prop_exposed <- n_e / n_total

  # for binary outcome (prevalence)
  if (is.null(y_t)) {

    # adjusted (logistic + standardisation) if covariates provided
    if (z != "") {

      if (verbose) cli::cli_alert("Fitting logistic regression for adjusted AFE/PAF (standardisation)")

      # coding checks on complete-case dataset
      if (!all(d[[x]] %in% c(0, 1))) stop("x must be coded 0/1 for adjusted AFE/PAF")
      if (!all(d[[y]] %in% c(0, 1))) stop("y must be coded 0/1 for adjusted AFE/PAF")

      # helper: compute adjusted risks + fractions from a fitted model and data
      calc_adj_from_glm <- function(glm_fit, d_fit, x, y) {
        p_obs <- stats::predict(glm_fit, type = "response")

        d_cf0 <- d_fit
        d_cf0[[x]] <- 0
        p_cf0 <- stats::predict(glm_fit, newdata = d_cf0, type = "response")

        d_cf1 <- d_fit
        d_cf1[[x]] <- 1
        p_cf1 <- stats::predict(glm_fit, newdata = d_cf1, type = "response")

        risk_obs <- mean(p_obs)
        risk_cf0 <- mean(p_cf0)
        risk_cf1 <- mean(p_cf1)

        paf <- 1 - (risk_cf0 / risk_obs)
        afe <- 1 - (risk_cf0 / risk_cf1)

        list(
          risk_observed = risk_obs,
          risk_counterfactual_unexposed = risk_cf0,
          risk_counterfactual_exposed = risk_cf1,
          paf = paf,
          afe = afe
        )
      }

      glm_formula <- stats::as.formula(paste0(y, " ~ ", x, z))
      glm_fit <- stats::glm(glm_formula, data = d, family = stats::binomial())

      adj <- calc_adj_from_glm(glm_fit = glm_fit, d_fit = d, x = x, y = y)

      paf <- adj$paf
      afe <- adj$afe
      excess_cases_exposed <- n_e * (adj$risk_counterfactual_exposed - adj$risk_counterfactual_unexposed)

      # bootstrap CIs (complete-case dataset)
      paf_ci_lower <- NA_real_
      paf_ci_upper <- NA_real_
      afe_ci_lower <- NA_real_
      afe_ci_upper <- NA_real_
      z_paf <- NA_real_
      z_afe <- NA_real_
      paf_p_value <- NA_real_
      afe_p_value <- NA_real_
      excess_cases_exposed_ci_lower <- NA_real_
      excess_cases_exposed_ci_upper <- NA_real_

      if (n_boot > 0L) {

        cli::cli_alert("Bootstrapping adjusted AFE/PAF CIs with {n_boot} iterations (can take a while)")

        n <- nrow(d)

        boot_res <- purrr::map_dfr(seq_len(n_boot), \(b) {
          boot_idx <- sample.int(n = n, size = n, replace = TRUE)
          d_boot <- d[boot_idx, , drop = FALSE]

          glm_boot <- stats::glm(glm_formula, data = d_boot, family = stats::binomial())
          adj_boot <- calc_adj_from_glm(glm_fit = glm_boot, d_fit = d_boot, x = x, y = y)
          n_exposed_boot <- sum(d_boot[[x]] == 1)

          data.frame(
            paf = adj_boot$paf,
            afe = adj_boot$afe,
            excess_cases_exposed_adj = n_exposed_boot * (adj_boot$risk_counterfactual_exposed - adj_boot$risk_counterfactual_unexposed)
          )
        })

        # bootstrap se and z (normal approximation)
        se_paf <- stats::sd(boot_res$paf, na.rm = TRUE)
        se_afe <- stats::sd(boot_res$afe, na.rm = TRUE)
        
        paf_ci_lower <- paf - (1.959964 * se_paf)
        paf_ci_upper <- paf + (1.959964 * se_paf)
        
        afe_ci_lower <- afe - (1.959964 * se_afe)
        afe_ci_upper <- afe + (1.959964 * se_afe)
        
        z_paf <- paf / se_paf
        z_afe <- afe / se_afe
        
        paf_p_value <- 2 * stats::pnorm(-abs(z_paf))
        afe_p_value <- 2 * stats::pnorm(-abs(z_afe))

        excess_cases_exposed_ci_lower <- stats::quantile(
          boot_res$excess_cases_exposed_adj, probs = 0.025, na.rm = TRUE, names = FALSE
        )
        excess_cases_exposed_ci_upper <- stats::quantile(
          boot_res$excess_cases_exposed_adj, probs = 0.975, na.rm = TRUE, names = FALSE
        )
      }

      res <- data.frame(
        # attributable fraction in exposed
        afe = afe,
        afe_ci_lower = afe_ci_lower,
        afe_ci_upper = afe_ci_upper,
        afe_percent = afe * 100,
        afe_ci_lower_percent = afe_ci_lower * 100,
        afe_ci_upper_percent = afe_ci_upper * 100,
        afe_z = z_afe,
        afe_p_value = afe_p_value,

        # population attributable fraction
        paf = paf,
        paf_ci_lower = paf_ci_lower,
        paf_ci_upper = paf_ci_upper,
        paf_percent = paf * 100,
        paf_ci_lower_percent = paf_ci_lower * 100,
        paf_ci_upper_percent = paf_ci_upper * 100,
        paf_z = z_paf,
        paf_p_value = paf_p_value,

        # supporting statistics
        n_total = n_total,
        n_exposed = n_e,
        n_unexposed = n_u,
        n_cases_exposed = n_cases_exposed,
        n_cases_unexposed = n_cases_unexposed,
        prop_exposed = prop_exposed,
        risk_observed = adj$risk_observed,
        risk_counterfactual_unexposed = adj$risk_counterfactual_unexposed,
        risk_counterfactual_exposed = adj$risk_counterfactual_exposed,
        excess_cases_exposed = excess_cases_exposed,
        excess_cases_exposed_ci_lower = excess_cases_exposed_ci_lower,
        excess_cases_exposed_ci_upper = excess_cases_exposed_ci_upper,
        x = x,
        y = y,
        z = z,
        n_boot = n_boot
      )

    } else {

      # crude RR method
      if (verbose) cli::cli_alert("Estimating Relative Risk (RR)\n")

      # calculate relative risk
      rr <- p_e / p_u

      # se of log(rr) using delta method
      se_log_rr <- sqrt((1 - p_e) / (n_e * p_e) + (1 - p_u) / (n_u * p_u))

      # 95% ci for rr
      rr_ci_lower <- exp(log(rr) - 1.959964 * se_log_rr)
      rr_ci_upper <- exp(log(rr) + 1.959964 * se_log_rr)

      if (verbose) cat("Estimating AFE and PAF\n")

      # attributable fraction in exposed: afe = (rr - 1) / rr
      afe <- (rr - 1) / rr
      afe_ci_lower <- (rr_ci_lower - 1) / rr_ci_lower
      afe_ci_upper <- (rr_ci_upper - 1) / rr_ci_upper

      # se of afe using delta method: se(afe) = se(log(rr)) / rr
      se_afe <- se_log_rr / rr
      z_afe <- afe / se_afe
      p_value_afe <- 2 * stats::pnorm(-abs(z_afe))

      # population attributable fraction: paf = p_e * (rr - 1) / (p_e * (rr - 1) + 1)
      paf <- prop_exposed * (rr - 1) / (prop_exposed * (rr - 1) + 1)
      paf_ci_lower <- prop_exposed * (rr_ci_lower - 1) / (prop_exposed * (rr_ci_lower - 1) + 1)
      paf_ci_upper <- prop_exposed * (rr_ci_upper - 1) / (prop_exposed * (rr_ci_upper - 1) + 1)

      # se of paf using delta method
      se_paf <- se_log_rr * prop_exposed / ((prop_exposed * (rr - 1) + 1) * rr)
      z_paf <- paf / se_paf
      p_value_paf <- 2 * stats::pnorm(-abs(z_paf))

      # number of excess cases in exposed
      excess_cases_exposed <- afe * n_cases_exposed
      excess_cases_exposed_ci_lower <- afe_ci_lower * n_cases_exposed
      excess_cases_exposed_ci_upper <- afe_ci_upper * n_cases_exposed

      res <- data.frame(
        # attributable fraction in exposed
        afe = afe,
        afe_ci_lower = afe_ci_lower,
        afe_ci_upper = afe_ci_upper,
        afe_percent = afe * 100,
        afe_ci_lower_percent = afe_ci_lower * 100,
        afe_ci_upper_percent = afe_ci_upper * 100,
        afe_z = z_afe,
        afe_p_value = p_value_afe,

        # population attributable fraction
        paf = paf,
        paf_ci_lower = paf_ci_lower,
        paf_ci_upper = paf_ci_upper,
        paf_percent = paf * 100,
        paf_ci_lower_percent = paf_ci_lower * 100,
        paf_ci_upper_percent = paf_ci_upper * 100,
        paf_z = z_paf,
        paf_p_value = p_value_paf,

        # supporting statistics
        n_total = n_total,
        n_exposed = n_e,
        n_unexposed = n_u,
        n_cases_exposed = n_cases_exposed,
        n_cases_unexposed = n_cases_unexposed,
        prop_exposed = prop_exposed,
        relative_risk = rr,
        rr_ci_lower = rr_ci_lower,
        rr_ci_upper = rr_ci_upper,
        risk_exposed = p_e,
        risk_unexposed = p_u,
        excess_risk = p_e - p_u,
        excess_cases_exposed = excess_cases_exposed,
        excess_cases_exposed_ci_lower = excess_cases_exposed_ci_lower,
        excess_cases_exposed_ci_upper = excess_cases_exposed_ci_upper,
        x = x,
        y = y
      )
    }
  }

  # for time-to-event outcome (incidence)
  if (!is.null(y_t)) {

    if (verbose) cli::cli_alert("Fitting Cox model\n")

    if (!y_t %in% colnames(d)) cat("!! Outcome time variable `y_t` not in the provided data frame\n")

    # fit cox model
    surv_obj <- survival::Surv(
      time = d[[y_t]],
      event = d[[y]]
    )

    cox_model <- survival::coxph(
      stats::as.formula(paste0("surv_obj ~ factor(", x, ")", z)),
      data = d
    )

    # hazard ratio, ci, and p-value
    hr <- exp(coef(cox_model)[1])
    hr_ci <- exp(confint(cox_model)[1, ])
    hr_p_value <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]

    if (verbose) cat("Transforming estimates to AFE and PAF\n")

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
    afe_p_value <- 2 * pnorm(-abs(z_afe))

    # population attributable fraction: paf = p_e * (hr - 1) / (p_e * (hr - 1) + 1)
    paf <- prop_exposed * (hr - 1) / (prop_exposed * (hr - 1) + 1)
    paf_ci_lower <- prop_exposed * (hr_ci[1] - 1) / (prop_exposed * (hr_ci[1] - 1) + 1)
    paf_ci_upper <- prop_exposed * (hr_ci[2] - 1) / (prop_exposed * (hr_ci[2] - 1) + 1)

    # se of paf using delta method
    se_paf <- se_log_hr * prop_exposed / ((prop_exposed * (hr - 1) + 1) * hr)
    z_paf <- paf / se_paf
    paf_p_value <- 2 * pnorm(-abs(z_paf))

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
      afe_z = z_afe,
      afe_p_value = afe_p_value,

      # population attributable fraction
      paf = paf,
      paf_ci_lower = paf_ci_lower,
      paf_ci_upper = paf_ci_upper,
      paf_percent = paf * 100,
      paf_ci_lower_percent = paf_ci_lower * 100,
      paf_ci_upper_percent = paf_ci_upper * 100,
      paf_z = z_paf,
      paf_p_value = paf_p_value,

      # supporting statistics
      n_total = n_total,
      n_exposed = n_e,
      n_unexposed = n_u,
      n_cases_exposed = n_cases_exposed,
      n_cases_unexposed = n_cases_unexposed,
      prop_exposed = prop_exposed,
      hazard_ratio = hr,
      hr_ci_lower = hr_ci[1],
      hr_ci_upper = hr_ci[2],
      hr_p_value = hr_p_value,
      rate_exposed_per_1000py = rate_exposed,
      rate_unexposed_per_1000py = rate_unexposed,
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
