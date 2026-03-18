#' Estimate haplotype frequencies and LD statistics from two biallelic variants
#'
#' @description Genotype data in R is usually loaded as genotype counts of biallelic variants. Accurately estimating haplotype frequencies without phase information is challenging. The key ambiguity is for double heterozygotes (1/1): you can't tell if the haplotypes are A1-B1/A2-B2 or A1-B2/A2-B1. 
#' 
#' To estimate haplotype frequency without the underlying haplotype phase we can use the Expectation-Maximisation (EM) algorithm by Excoffier & Slatkin (1995) (https://doi.org/10.1093/oxfordjournals.molbev.a040269). This function estimates the haplotype frequencies and returns the R^2 and D' statistics.
#' 
#' Estimates highly similar to `plink2 --r2-phased cols=+dprimeabs` for several examples. This function is the combination of some of my older scripts plus input from Copilot.
#' 
#' @return A data.frame with haplotype freqs, allele freqs, D, D', and R²
#'
#' @author Luke Pilling
#'
#' @name estimate_ld
#'
#' @param df A data frame of exactly two cols - both must be integer vector (0, 1, 2) — minor allele counts for variant A and B
#' @param tol Convergence tolerance for EM
#' @param max_iter Maximum EM iterations
#'
#' @examples
#'
#' # Estimate LD between HFE p.C282Y and HFE p.S65C
#' ukb |> 
#'   dplyr::select(rs1800562_g_a, rs1800730_a_t) |> 
#'   estimate_ld()
#' #>         P_A1B1     P_A1B2     P_A2B1    P_A2B2       p_A1      p_A2       p_B1     p_B2            D       D_max    D_prime unsigned_D_prime   R_squared em_iterations
#' #> 1 3.203526e-05 0.07335235 0.01563897 0.9109766 0.07338439 0.9266156 0.01567101 0.984329 -0.001117972 0.001150007 -0.9721434        0.9721434 0.001191575           110
#'
#' @export
#'
estimate_ld <- function(df, tol = 1e-8, max_iter = 1000) {

  # check provided data frame 
  if (! any(class(df) %in% c("data.frame","tbl","tbl_df")))  {
    cli::cli_abort(c(
      "{.var df} must be a data.frame (or tibble)",
      "x" = "You've supplied a {.cls {class(df)}}."
    ))
  }
  if (ncol(df)!=2)  {
    cli::cli_abort(c(
	  "Provided data frame must have exactly two cols. Both must be integer vectors (0, 1, 2) — minor allele counts for variant A and B",
      "x" = "{.var df} has {ncol(df)} cols."
	))
  }
  
  # drop missing values
  df <- as.data.frame(na.omit(df))
  
  # assign vectors
  geno_a <- as.vector(df[,1])
  geno_b <- as.vector(df[,2])
  
  # --- EM algorithm to estimate haplotype frequencies ---

  # Count the 9 genotype combinations (0/0 through 2/2)
  counts <- table(factor(geno_a, levels = 0:2),
                  factor(geno_b, levels = 0:2))
  n <- sum(counts)

  # Initialise assuming linkage equilibrium
  p_a <- mean(geno_a) / 2
  p_b <- mean(geno_b) / 2
  h <- c(p_a * p_b, p_a * (1 - p_b), (1 - p_a) * p_b, (1 - p_a) * (1 - p_b))

  for (iter in seq_len(max_iter)) {
    h_old <- h

    # E-step: only double heterozygotes (1,1) have ambiguous phase
    denom  <- h[1] * h[4] + h[2] * h[3]
    p_coup <- h[1] * h[4] / denom  # prob of coupling arrangement
    n_11   <- counts[2, 2]         # count of (1,1) genotypes

    # Tally haplotype contributions from all 9 genotype classes
    e_h <- c(
      2 * counts[3, 3] + counts[3, 2] + counts[2, 3] + n_11 * p_coup,
      2 * counts[3, 1] + counts[3, 2] + counts[2, 1] + n_11 * (1 - p_coup),
      2 * counts[1, 3] + counts[1, 2] + counts[2, 3] + n_11 * (1 - p_coup),
      2 * counts[1, 1] + counts[1, 2] + counts[2, 1] + n_11 * p_coup
    )

    # M-step
    h <- e_h / (2 * n)

    if (max(abs(h - h_old)) < tol) break
  }

  # --- LD statistics from estimated haplotype frequencies ---

  P_A1B1 <- h[1]; P_A1B2 <- h[2]; P_A2B1 <- h[3]; P_A2B2 <- h[4]

  # Allele frequencies
  p_A1 <- P_A1B1 + P_A1B2
  p_A2 <- P_A2B1 + P_A2B2
  p_B1 <- P_A1B1 + P_A2B1
  p_B2 <- P_A1B2 + P_A2B2

  # Check haplotype frequencies sum to 1
  if (abs(sum(h) - 1) > 1e-6) {
    cli::cli_warn("Haplotype frequencies do not sum to 1.")
  }

  # D
  D <- P_A1B1 - (p_A1 * p_B1)

  # D_max and D'
  D_max <- if (D >= 0) min(p_A1 * p_B2, p_A2 * p_B1) else min(p_A1 * p_B1, p_A2 * p_B2)

  D_prime <- if (D_max == 0) {
    cli::cli_warn("D_max is zero, D_prime cannot be calculated (possible fixed allele).")
    NA_real_
  } else {
    D / D_max
  }

  # R^2
  denom_r2 <- p_A1 * p_A2 * p_B1 * p_B2
  R_squared <- if (denom_r2 == 0) {
    cli::cli_warn("Denominator for R_squared is zero (possible fixed allele).")
    NA_real_
  } else {
    D^2 / denom_r2
  }
  
  # return results as row of data frame
  data.frame(
    P_A1B1, P_A1B2, P_A2B1, P_A2B2,
    p_A1, p_A2, p_B1, p_B2,
	N=nrow(df),
    D, D_max,
    D_prime,
    unsigned_D_prime = abs(D_prime),
    R_squared,
    em_iterations = iter
  )
}




