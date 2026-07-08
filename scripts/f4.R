library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# f4 / D-statistics: f4(W, X; TARGET, Y)
# Tests whether TARGET has excess affinity with X vs W, relative to Y.
# A significant positive f4 means TARGET shares more drift with X than with W.
# Used to detect ancestry signals not captured by the qpAdm model.

# Standard test: does TARGET have excess Steppe affinity compared to a reference?
# f4(Mbuti, Steppe; TARGET, EEF) — positive = more steppe-like than expected

OUTGROUP <- "Mbuti"

# Test set: compare TARGET against each reference population
# using key ancient populations as the contrast
STEPPE   <- "Russia_Samara_EBA_Yamnaya"
WHG      <- "Luxembourg_Loschbour_Mesolithic"
EEF      <- "Turkey_N"

all_pops <- c(
  TARGET, OUTGROUP, STEPPE, WHG, EEF,
  REFERENCES,
  "Han", "Yoruba", "Papuan",
  "Iran_GanjDareh_N", "Israel_Natufian"
)

cat("=== f4 / D-statistics ===\n")
cat("Target:", TARGET, "\n\n")

if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks (first run — may take several minutes)...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR)


# Test 1: Steppe vs EEF balance across populations.
# f4(Mbuti, X; Steppe, EEF): positive z = X is more EEF-shifted (Turkey_N),
# negative z = X is more Steppe-shifted (Yamnaya).
cat("--- Steppe/EEF balance: f4(Mbuti, X; Steppe, EEF) ---\n")
cat("Positive z = more EEF-shifted (Turkey_N); negative z = more Steppe-shifted\n\n")

steppe_tests <- f4(f2,
  pop1 = OUTGROUP, pop2 = c(TARGET, REFERENCES),
  pop3 = STEPPE,   pop4 = EEF
) %>%
  rename(population = pop2, f4 = est) %>%
  select(population, f4, se, z) %>%
  arrange(desc(z))

cat(sprintf("  %-38s  %10s  %8s  %6s\n", "Population", "f4", "se", "z"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(steppe_tests))) {
  r <- steppe_tests[i, ]
  marker <- if (r$population == TARGET) " <-- you" else ""
  cat(sprintf("  %-38s  %10.6f  %8.6f  %6.2f%s\n",
              r$population, r$f4, r$se, r$z, marker))
}


# Test 2: WHG affinity — f4(Mbuti, X; WHG, EEF)
cat("\n--- WHG/EEF balance: f4(Mbuti, X; WHG, EEF) ---\n")
cat("Positive z = more EEF-shifted; negative z = more WHG-shifted\n\n")

whg_tests <- f4(f2,
  pop1 = OUTGROUP, pop2 = c(TARGET, REFERENCES),
  pop3 = WHG,      pop4 = EEF
) %>%
  rename(population = pop2, f4 = est) %>%
  select(population, f4, se, z) %>%
  arrange(desc(z))

cat(sprintf("  %-38s  %10s  %8s  %6s\n", "Population", "f4", "se", "z"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(whg_tests))) {
  r <- whg_tests[i, ]
  marker <- if (r$population == TARGET) " <-- you" else ""
  cat(sprintf("  %-38s  %10.6f  %8.6f  %6.2f%s\n",
              r$population, r$f4, r$se, r$z, marker))
}

cat("\nDone.\n")
