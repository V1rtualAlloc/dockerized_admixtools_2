library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# f4 / D-statistics: f4(W, X; TARGET, Y)
# Tests whether TARGET has excess affinity with X vs W, relative to Y.
# A significant positive f4 means TARGET shares more drift with X than with W.
# Used to detect ancestry signals not captured by the qpAdm model.
#
# Runs one balance test per pair of the background's own qpAdm sources
# (config.R MODELS[[BACKGROUND]]$sources), e.g. for the european 3-source
# model (Steppe/Balkan HG/Anatolian EEF) this reproduces the original
# Steppe-vs-EEF and HG-vs-EEF balance tests, generalized to any BACKGROUND.

OUTGROUP <- "Mbuti"
model    <- MODELS[[BACKGROUND]]
sources  <- model$sources

all_pops <- unique(c(
  TARGET, OUTGROUP, sources,
  REFERENCES,
  WORLD_REFS
))

cat("=== f4 / D-statistics ===\n")
cat("Target:", TARGET, "\n")
cat("Sources:", paste(sources, collapse = ", "), "\n\n")

if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks (first run — may take several minutes)...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR)


for (pair in combn(sources, 2, simplify = FALSE)) {
  s3 <- pair[1]
  s4 <- pair[2]
  cat(sprintf("--- %s/%s balance: f4(%s, X; %s, %s) ---\n", s3, s4, OUTGROUP, s3, s4))
  cat(sprintf("Positive z = more %s-shifted; negative z = more %s-shifted\n\n", s4, s3))

  tests <- f4(f2,
    pop1 = OUTGROUP, pop2 = c(TARGET, REFERENCES),
    pop3 = s3,       pop4 = s4
  ) %>%
    rename(population = pop2, f4 = est) %>%
    select(population, f4, se, z) %>%
    arrange(desc(z))

  cat(sprintf("  %-38s  %10s  %8s  %6s\n", "Population", "f4", "se", "z"))
  cat(strrep("-", 70), "\n")
  for (i in seq_len(nrow(tests))) {
    r <- tests[i, ]
    marker <- if (r$population == TARGET) " <-- you" else ""
    cat(sprintf("  %-38s  %10.6f  %8.6f  %6.2f%s\n",
                r$population, r$f4, r$se, r$z, marker))
  }
  cat("\n")
}

cat("Done.\n")
