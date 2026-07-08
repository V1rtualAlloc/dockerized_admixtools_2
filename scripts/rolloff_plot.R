library(tidyverse)

OUTDIR    <- "/data/rolloff/output"
OUTFILE   <- file.path(OUTDIR, "me.out")
LOVALFIT  <- 0.45   # cM — starting distance for fit (DATES default)
HIVALFIT  <- 150.0  # cM — max distance for fit
YEARS_PER_GEN <- 29

cat("=== DATES LD dating — exponential fit ===\n\n")

read_dates_out <- function(f) {
  raw <- readLines(f)
  data_lines <- raw[!grepl("^\\s*#|^\\s*$", raw)]
  read.table(text = paste(data_lines, collapse = "\n"),
             col.names = c("dist_cM", "cov", "z", "col4", "npairs"))
}

dat <- read_dates_out(OUTFILE)
cat(sprintf("Main output: %d distance bins (%.2f – %.2f cM)\n\n",
            nrow(dat), min(dat$dist_cM), max(dat$dist_cM)))

fit_exp <- function(df, loval = LOVALFIT, hival = HIVALFIT) {
  df <- df %>% filter(dist_cM >= loval, dist_cM <= hival, is.finite(cov))
  if (nrow(df) < 10) return(NULL)
  tryCatch(
    nls(cov ~ A * exp(-n * dist_cM / 100) + B,
        data    = df,
        start   = list(A = max(df$cov), n = 100, B = min(df$cov)),
        control = nls.control(maxiter = 500, warnOnly = TRUE)),
    error = function(e) NULL
  )
}

# Main fit
fit <- fit_exp(dat)
if (is.null(fit)) stop("Exponential fit failed — check me.out format")

cf     <- coef(fit)
n_main <- cf["n"]

# Jackknife SE: fit each chromosome-leave-one-out file, then jackknife formula
jack_files <- list.files(OUTDIR, pattern = "^me\\.out:\\d+$", full.names = TRUE)
jack_n <- numeric(length(jack_files))
for (i in seq_along(jack_files)) {
  d <- tryCatch(read_dates_out(jack_files[i]), error = function(e) NULL)
  f <- if (!is.null(d)) fit_exp(d) else NULL
  jack_n[i] <- if (!is.null(f)) coef(f)["n"] else NA
}
jack_n <- jack_n[!is.na(jack_n)]
J      <- length(jack_n)
jack_mean <- mean(jack_n)
jack_se   <- sqrt((J - 1) / J * sum((jack_n - jack_mean)^2))

n_est    <- n_main
se_est   <- jack_se
n_years  <- n_est * YEARS_PER_GEN
se_years <- se_est * YEARS_PER_GEN
date_bce <- round(n_years - 2024)

cat(sprintf("Estimated generations since admixture:  %.1f ± %.1f  (jackknife SE, %d chrom)\n",
            n_est, se_est, J))
cat(sprintf("In years:                               %.0f ± %.0f years ago\n",
            n_years, se_years))
cat(sprintf("Approximate calendar date:              ~%.0f BCE  (95%% CI: ~%.0f – ~%.0f BCE)\n\n",
            date_bce,
            round((n_years + 2 * se_years) - 2024),
            round((n_years - 2 * se_years) - 2024)))
cat("Note: single-individual LD estimates have wide uncertainty.\n")
cat("Interpret as order-of-magnitude: Neolithic vs. Bronze Age vs. Iron Age.\n\n")

# Plot
fit_range <- dat %>% filter(dist_cM >= LOVALFIT, dist_cM <= HIVALFIT)
pred <- tibble(dist_cM = seq(LOVALFIT, HIVALFIT, length.out = 500),
               cov     = cf["A"] * exp(-cf["n"] * dist_cM / 100) + cf["B"])

p <- ggplot(dat %>% filter(dist_cM <= HIVALFIT), aes(dist_cM, cov)) +
  geom_point(size = 1.2, alpha = 0.5, color = "steelblue") +
  geom_line(data = pred, color = "firebrick", linewidth = 1) +
  geom_vline(xintercept = LOVALFIT, linetype = "dashed", color = "gray60") +
  annotate("text", x = HIVALFIT * 0.55, y = max(dat$cov[dat$dist_cM >= LOVALFIT]) * 0.85,
           label = sprintf("n = %.0f ± %.0f gen\n~%.0f BCE",
                           n_est, se_est, date_bce),
           hjust = 0, size = 4.5, color = "firebrick") +
  labs(
    title   = "DATES: LD decay dating (Steppe + EEF admixture)",
    x       = "Genetic distance (cM)",
    y       = "Weighted LD covariance",
    caption = "Sources: Russia_Samara_EBA_Yamnaya + Turkey_N  |  Target: me\nDashed line = start of fit range (0.45 cM)"
  ) +
  theme_bw(base_size = 12)

outfile <- file.path(OUTDIR, "dates_fit.pdf")
ggsave(outfile, p, width = 10, height = 7)
cat("Plot saved:", outfile, "\n\nDone.\n")
