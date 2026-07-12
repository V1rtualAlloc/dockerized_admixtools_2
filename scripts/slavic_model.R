library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Slavic / pre-Slavic ancestry model — era-matched Roman decomposition
#
# Core insight: Roman Balkans = Iron Age Balkans + Roman-era immigrants.
# We use era-matched sources for those immigrants so they don't confound
# with the Iron Age substrate or with the Neolithic-era outgroups:
#
#   Slavic migrants       — Poland_EarlyMedieval_Slav    (6th–9th c. CE)
#   Iron Age Balkans      — Croatia_EIA                  (Illyrian/Hallstatt, ~800–400 BCE)
#   Roman-era Anatolian   — Turkey_Medieval_Byzantine    (Byzantine Anatolia, ~6th–12th c. CE)
#   Roman-era Levantine   — Israel_Phoenician            (Iron Age Levant, ~900–500 BCE)
#   Roman Italian         — Italy_Lazio_ImperialRoman_Roman (1st–3rd c. CE)
#
# Outgroups: Turkey_N and Israel_Natufian are Neolithic (8000+ BCE) — too old
# to overlap with any source, so they anchor the Near Eastern rotation throughout
# without needing to swap per model.

OUTGROUPS <- SLAVIC_OUTGROUPS

# Confirmed models pool
ALL_SOURCES <- c(
  "Poland_EarlyMedieval_Slav",
  "Serbia_LateAntiquity_ImperialRoman",
  "Croatia_EIA",
  "Turkey_Medieval_Byzantine",
  "Israel_Phoenician",
  "Italy_Lazio_ImperialRoman_Roman"
)

# Reference populations for cross-population comparison
SLAVIC_REFS <- c(
  TARGET
)

# Three-role rotating pool:
#   Slavic migrant  x  Balkan Iron Age substrate  x  non-Balkan Roman/Byzantine import
# Each combination has exactly one population from each role.

SLAVIC_POOL <- SLAVIC_MIGRANT_POOL

BALKAN_IA_POOL <- SLAVIC_BALKAN_IA_CORE

ROMAN_POOL <- SLAVIC_ROMAN_FULL

ROTATE_POOL <- unique(c(SLAVIC_POOL, BALKAN_IA_POOL, ROMAN_POOL))

all_pops <- unique(c(TARGET, ALL_SOURCES, ROTATE_POOL, SLAVIC_REFS, OUTGROUPS))

SLAVIC_F2_DIR <- sub("f2/?$", "f2_slavic/", F2_DIR)

cat("=== Slavic / pre-Slavic ancestry model ===\n")
cat("Target:", TARGET, "\n\n")

cached <- list.dirs(SLAVIC_F2_DIR, full.names = FALSE, recursive = FALSE)
if (!all(all_pops %in% cached)) {
  cat("Building f2 cache (may take a few minutes)...\n\n")
  dir.create(SLAVIC_F2_DIR, showWarnings = FALSE, recursive = TRUE)
  extract_f2(MERGED_PREFIX, outdir = SLAVIC_F2_DIR, pops = all_pops, overwrite = TRUE)
} else {
  cat("f2 blocks cached, loading...\n\n")
}
f2 <- f2_from_precomp(SLAVIC_F2_DIR, pops = all_pops)


run_model <- function(sources, label) {
  r <- tryCatch(
    qpadm(f2, target = TARGET, left = sources, right = OUTGROUPS),
    error = function(e) { cat("  [error:", conditionMessage(e), "]\n\n"); NULL }
  )
  if (is.null(r)) return(invisible(NULL))

  full <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  p    <- full$p[1]

  w_col    <- intersect(c("weight", "w", "est"), names(r$weights))
  se_col   <- intersect(c("se", "se_weight"), names(r$weights))
  z_col    <- intersect(c("z", "z_weight"), names(r$weights))
  wts      <- setNames(r$weights[[w_col[1]]], r$weights$left)
  ses      <- if (length(se_col) > 0) setNames(r$weights[[se_col[1]]], r$weights$left) else rep(NA, length(wts))
  zs       <- if (length(z_col)  > 0) setNames(r$weights[[z_col[1]]],  r$weights$left) else wts / ses
  feasible <- all(wts >= -0.01 & wts <= 1.01)

  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  p = %.4f  |  %s  |  weights %s\n\n",
              p,
              if (p > 0.05) "PASS" else "fail",
              if (feasible) "feasible" else "INFEASIBLE"))
  cat(sprintf("  %-46s  %8s  %8s  %6s\n", "source", "weight", "se", "z"))
  for (nm in names(wts)) {
    cat(sprintf("  %-46s  %+7.1f%%  %7.1f%%  %+6.2f\n",
                nm, wts[nm] * 100, ses[nm] * 100, zs[nm]))
  }
  cat("\n")
  invisible(list(p = p, wts = wts, feasible = feasible))
}


cat("=== 2-source baselines ===\n\n")

run_model(
  c("Poland_EarlyMedieval_Slav", "Serbia_LateAntiquity_ImperialRoman"),
  "Model A — Slavic + Roman Balkans [bundled, for comparison]"
)
run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA"),
  "Model B — Slavic + Iron Age Balkans [cleanest baseline]"
)


cat("=== 3-source: add one era-matched Roman import ===\n\n")

run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA", "Turkey_Medieval_Byzantine"),
  "Model C — Slavic + Iron Age Balkans + Byzantine Anatolia"
)
run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA", "Israel_Phoenician"),
  "Model D — Slavic + Iron Age Balkans + Levantine [Israel_Phoenician]"
)
run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA", "Italy_Lazio_ImperialRoman_Roman"),
  "Model E — Slavic + Iron Age Balkans + Roman Italian"
)


cat("=== 4-source: two Roman import components ===\n\n")

run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA",
    "Turkey_Medieval_Byzantine", "Italy_Lazio_ImperialRoman_Roman"),
  "Model F — Slavic + Iron Age Balkans + Byzantine Anatolia + Roman Italian"
)
run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA",
    "Turkey_Medieval_Byzantine", "Israel_Phoenician"),
  "Model G — Slavic + Iron Age Balkans + Byzantine Anatolia + Levantine"
)
run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA",
    "Israel_Phoenician", "Italy_Lazio_ImperialRoman_Roman"),
  "Model H — Slavic + Iron Age Balkans + Levantine + Roman Italian"
)


cat("=== 5-source: full Roman decomposition ===\n\n")

run_model(
  c("Poland_EarlyMedieval_Slav", "Croatia_EIA",
    "Turkey_Medieval_Byzantine", "Israel_Phoenician", "Italy_Lazio_ImperialRoman_Roman"),
  "Model I — Slavic + Iron Age Balkans + Byzantine Anatolia + Levantine + Roman Italian"
)


cat("=== qpAdm rotating: Slavic x Balkan-IA x Roman/Byzantine ===\n\n")
cat("Slavic pool:    ", paste(SLAVIC_POOL,    collapse=", "), "\n")
cat("Balkan IA pool: ", paste(BALKAN_IA_POOL, collapse=", "), "\n")
cat("Roman pool:     ", paste(ROMAN_POOL,     collapse=", "), "\n\n")

# 3-source: one from each role
combos_3 <- do.call(c, lapply(SLAVIC_POOL, function(sl)
  lapply(BALKAN_IA_POOL, function(ia)
    lapply(ROMAN_POOL, function(ro) c(sl, ia, ro))))) |> unlist(recursive=FALSE)

# 2-source baselines: Slavic + IA only, and Slavic + Roman only
combos_2 <- c(
  do.call(c, lapply(SLAVIC_POOL, function(sl) lapply(BALKAN_IA_POOL, function(ia) c(sl, ia)))),
  do.call(c, lapply(SLAVIC_POOL, function(sl) lapply(ROMAN_POOL,     function(ro) c(sl, ro))))
)

combos <- c(combos_2, combos_3)
cat(sprintf("Testing %d models (%d 2-source baselines, %d 3-source)...\n\n",
            length(combos), length(combos_2), length(combos_3)))

run_one <- function(sources) {
  r <- tryCatch(
    qpadm(f2, target = TARGET, left = sources, right = OUTGROUPS),
    error = function(e) NULL
  )
  if (is.null(r) || is.null(r$popdrop) || is.null(r$weights)) return(NULL)
  full <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  if (nrow(full) == 0) return(NULL)
  w_col  <- intersect(c("weight", "w", "est"), names(r$weights))
  se_col <- intersect(c("se", "se_weight"),    names(r$weights))
  wts    <- setNames(r$weights[[w_col[1]]], r$weights$left)
  ses    <- if (length(se_col)) setNames(r$weights[[se_col[1]]], r$weights$left) else rep(NA_real_, length(wts))
  feasible <- all(wts >= -0.01 & wts <= 1.01)
  w_str  <- paste(sprintf("%s=%.0f%%±%.0f%%", names(wts), wts*100, ses*100), collapse="  ")
  tibble(n=length(sources), model=paste(sources, collapse=" + "),
         p=full$p, feasible=feasible, weights=w_str)
}

rot_res  <- map_dfr(combos, run_one) %>% arrange(desc(p))
passing  <- filter(rot_res, p > 0.05, feasible)
failing  <- filter(rot_res, p <= 0.05 | !feasible)

cat(sprintf("=== PASSING (p > 0.05, feasible): %d / %d ===\n\n",
            nrow(passing), nrow(rot_res)))
if (nrow(passing) > 0) {
  cat(sprintf("  %3s  %7s  %s\n", "src", "p", "weights"))
  cat(strrep("-", 120), "\n")
  for (i in seq_len(nrow(passing))) {
    r <- passing[i, ]
    cat(sprintf("  %3d  %7.4f  %s\n", r$n, r$p, r$weights))
  }
}

cat(sprintf("\n=== FAILING (p <= 0.05 or infeasible): %d / %d ===\n\n",
            nrow(failing), nrow(rot_res)))
cat(sprintf("  %3s  %7s  %-10s  %s\n", "src", "p", "feasible", "model"))
cat(strrep("-", 100), "\n")
for (i in seq_len(nrow(failing))) {
  r <- failing[i, ]
  cat(sprintf("  %3d  %7.4f  %-10s  %s\n", r$n, r$p,
              if (r$feasible) "" else "[infeasible]", r$model))
}
cat("\n")

cat("Done.\n  (Reference population comparison: run slavic_modern.R)\n")
