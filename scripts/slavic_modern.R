library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Tests the best Slavic/pre-Slavic models against modern Balkan and European
# reference populations. Uses a focused f2 cache (no small ancient populations
# with high missingness) to keep SEs tight.

OUTGROUPS <- SLAVIC_OUTGROUPS

ANCIENT_SOURCES <- c(
  "Poland_EarlyMedieval_Slav",
  "Croatia_EIA",
  "NorthMacedonia_IA",
  "Bulgaria_KapitanAndreevo_EIA",
  "Italy_Lazio_ImperialRoman_Roman",
  "Turkey_Medieval_Byzantine"
)

MODERN_REFS <- c(
  "Bulgarian", "Albanian", "Greek", "Greek_1", "Greek_Crete",
  "Hungarian", "Czech", "Polish", "Russian"
)

F2_MODERN_DIR <- sub("f2/?$", "f2_modern/", F2_DIR)

all_pops <- unique(c(MODERN_REFS, ANCIENT_SOURCES, OUTGROUPS))

cached <- list.dirs(F2_MODERN_DIR, full.names=FALSE, recursive=FALSE)
if (!all(all_pops %in% cached)) {
  cat("Building f2 cache...\n\n")
  dir.create(F2_MODERN_DIR, showWarnings=FALSE, recursive=TRUE)
  extract_f2(MERGED_PREFIX, outdir=F2_MODERN_DIR, pops=all_pops, overwrite=TRUE)
} else {
  cat("f2 cache loaded.\n\n")
}
f2 <- f2_from_precomp(F2_MODERN_DIR, pops=all_pops)

# Best models to test across modern populations
MODELS <- list(
  list(label="Slav + Croatia_EIA [conservative]",
       sources=c("Poland_EarlyMedieval_Slav", "Croatia_EIA")),
  list(label="Slav + NorthMacedonia_IA [central Balkans]",
       sources=c("Poland_EarlyMedieval_Slav", "NorthMacedonia_IA")),
  list(label="Slav + Bulgaria_KapitanAndreevo_EIA [Thracian]",
       sources=c("Poland_EarlyMedieval_Slav", "Bulgaria_KapitanAndreevo_EIA")),
  list(label="Slav + Croatia_EIA + Italy Roman [3-source]",
       sources=c("Poland_EarlyMedieval_Slav", "Croatia_EIA", "Italy_Lazio_ImperialRoman_Roman"))
)

run_model_on_pops <- function(sources, label) {
  src1 <- sources[1]
  src2 <- sources[2]
  src3 <- if (length(sources) == 3) sources[3] else NULL

  rows <- map_dfr(MODERN_REFS, function(pop) {
    r <- tryCatch(
      qpadm(f2, target=pop, left=sources, right=OUTGROUPS),
      error=function(e) NULL
    )
    if (is.null(r) || is.null(r$weights)) return(NULL)
    full   <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
    w_col  <- intersect(c("weight","w","est"), names(r$weights))
    se_col <- intersect(c("se","se_weight"),   names(r$weights))
    wts <- setNames(r$weights[[w_col[1]]], r$weights$left)
    ses <- setNames(r$weights[[se_col[1]]], r$weights$left)
    tibble(
      pop      = pop,
      slavic   = wts[src1], slavic_se = ses[src1],
      sub      = wts[src2], sub_se    = ses[src2],
      third    = if (!is.null(src3)) wts[src3] else NA_real_,
      third_se = if (!is.null(src3)) ses[src3] else NA_real_,
      p        = full$p[1],
      feasible = all(wts >= -0.01 & wts <= 1.01)
    )
  }) %>% arrange(desc(slavic))

  cat(sprintf("=== %s ===\n\n", label))
  if (is.null(src3)) {
    cat(sprintf("  %-20s  %8s  %6s  %8s  %6s  %7s  %s\n",
                "population", "slavic%", "+-SE", "substrate%", "+-SE", "p", "fit"))
    cat(strrep("-", 75), "\n")
    for (i in seq_len(nrow(rows))) {
      r   <- rows[i,]
      mrk <- if (r$pop == TARGET) " <-- you" else ""
      fit <- if (!r$feasible) "INFEASIBLE" else if (r$p > 0.05) "PASS" else "fail"
      cat(sprintf("  %-20s  %7.1f%%  %5.1f%%  %8.1f%%  %5.1f%%  %7.4f  %s%s\n",
                  r$pop, r$slavic*100, r$slavic_se*100,
                  r$sub*100, r$sub_se*100, r$p, fit, mrk))
    }
  } else {
    cat(sprintf("  %-20s  %8s  %6s  %8s  %6s  %8s  %6s  %7s  %s\n",
                "population", "slavic%", "+-SE", "substrate%", "+-SE",
                "roman%", "+-SE", "p", "fit"))
    cat(strrep("-", 98), "\n")
    for (i in seq_len(nrow(rows))) {
      r   <- rows[i,]
      mrk <- if (r$pop == TARGET) " <-- you" else ""
      fit <- if (!r$feasible) "INFEASIBLE" else if (r$p > 0.05) "PASS" else "fail"
      cat(sprintf("  %-20s  %7.1f%%  %5.1f%%  %8.1f%%  %5.1f%%  %7.1f%%  %5.1f%%  %7.4f  %s%s\n",
                  r$pop, r$slavic*100, r$slavic_se*100,
                  r$sub*100, r$sub_se*100,
                  r$third*100, r$third_se*100, r$p, fit, mrk))
    }
  }
  cat("\n")
}

cat("=== Slavic model: modern population comparison ===\n\n")
cat("Ancient sources: Poland_EarlyMedieval_Slav + Balkan IA proxy\n")
cat("Note: modern pops have N=1-2 individuals -- interpret SEs with caution\n\n")

for (m in MODELS) run_model_on_pops(m$sources, m$label)

cat("Done.\n")
