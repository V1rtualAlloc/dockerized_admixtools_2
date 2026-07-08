library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# MDS from pairwise f2 distances — richer population set.
# Uses a separate f2 cache (F2_MDS_DIR) so the analysis cache is not affected.

F2_MDS_DIR <- file.path(dirname(F2_DIR), "f2_mds")
dir.create(F2_MDS_DIR, showWarnings = FALSE, recursive = TRUE)

model <- MODELS[[BACKGROUND]]

# Extra ancient/medieval populations for finer regional resolution. Only
# defined for "european" (this project's original Balkan/Slavic use case) —
# add an entry here for other backgrounds if you want the same fine-grained
# view; otherwise MDS_POPS below falls back to the background's own
# sources/outgroup/pool, which is still a valid (just coarser) MDS.
REGIONAL_EXTRA <- list(
  european = c(
    "Serbia_IronGates_Mesolithic",      # Balkan Mesolithic HG
    "Serbia_EN_Starcevo",               # first Balkan farmers
    "Serbia_EBA_Maros",                 # Balkan Bronze Age
    "Serbia_ImperialRoman",             # Roman-era Serbia
    "Bulgaria_Varna_C",                 # Chalcolithic Balkans
    "NorthMacedonia_IA",                # Iron Age Macedonia
    "Croatia_MLBA",                     # Bronze Age Croatia
    "Croatia_EarlyMedieval_EarlySlav",  # Slavic migration era
    "Ukraine_N",                        # Neolithic farmers
    "Czechia_EBA_CordedWare",           # Corded Ware
    "Czechia_BellBeaker",               # Bell Beaker
    "Czechia_EBA_Unetice",              # Unetice Bronze Age
    "Poland_EarlyMedieval_Slav",        # early medieval Slavs
    "Sweden_Viking"                     # Vikings
  )
)[[BACKGROUND]]

MDS_POPS <- unique(c(
  TARGET,
  if (!is.null(model$pool)) model$pool else model$sources,
  model$outgroup,
  REGIONAL_EXTRA,
  REFERENCES,
  WORLD_REFS,
  "China_TianyuanCave_UP", "Ethiopia_MotaCave_4500BP"
))

cat("=== MDS from pairwise f2 distances (richer population set) ===\n")
cat("Target:", TARGET, "\n")
cat("Populations:", length(MDS_POPS), "\n\n")

if (length(list.files(F2_MDS_DIR, pattern = "\\.rds$")) == 0) {
  cat("Building f2 cache for MDS populations (this takes a few minutes)...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_MDS_DIR, pops = MDS_POPS, overwrite = FALSE)
  cat("Done.\n\n")
} else {
  cat("f2 cache found, loading...\n")
}
f2 <- f2_from_precomp(F2_MDS_DIR, pops = MDS_POPS)

# Average f2 over blocks → pairwise distance matrix
dmat <- apply(f2, 1:2, mean, na.rm = TRUE)
dmat <- (dmat + t(dmat)) / 2
diag(dmat) <- 0
dmat[dmat < 0] <- 0

pops_all <- rownames(dmat)

ancient <- unique(c(
  if (!is.null(model$pool)) model$pool else model$sources,
  REGIONAL_EXTRA,
  "China_TianyuanCave_UP", "Ethiopia_MotaCave_4500BP"
))
outgroups <- unique(c(model$outgroup, "Mbuti", "Yoruba", "Papuan", "Han",
                      "China_TianyuanCave_UP", "Ethiopia_MotaCave_4500BP"))

make_mds_df <- function(dm) {
  mds  <- cmdscale(dm, k = 2, eig = TRUE)
  eig_pos <- mds$eig[mds$eig > 0]
  pct  <- round(mds$eig[1:2] / sum(eig_pos) * 100, 1)
  df   <- as.data.frame(mds$points) %>%
    setNames(c("MDS1", "MDS2")) %>%
    rownames_to_column("pop") %>%
    mutate(group = case_when(
      pop == TARGET                        ~ "You",
      pop %in% REFERENCES                  ~ "Reference",
      pop %in% outgroups                   ~ "Outgroup",
      pop %in% ancient                     ~ "Ancient",
      TRUE                                 ~ "Modern/Other"
    ))
  # Remove outliers > 3 SD from center (small-n populations)
  df <- df %>%
    filter(
      abs(MDS1 - mean(MDS1)) < 3 * sd(MDS1),
      abs(MDS2 - mean(MDS2)) < 3 * sd(MDS2)
    )
  list(df = df, pct = pct)
}

make_plot <- function(df, pct, title) {
  ggplot(df, aes(MDS1, MDS2, color = group)) +
    geom_point(aes(size = (pop == TARGET), shape = group)) +
    scale_size_manual(values = c("TRUE" = 5, "FALSE" = 2.5), guide = "none") +
    geom_text(aes(label = pop), size = 2.2, hjust = -0.1,
              vjust = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c(
      "You"          = "red",
      "Reference"    = "steelblue",
      "Ancient"      = "darkorange",
      "Modern/Other" = "forestgreen",
      "Outgroup"     = "gray50"
    )) +
    scale_shape_manual(values = c(
      "You" = 18, "Reference" = 16, "Ancient" = 17,
      "Modern/Other" = 15, "Outgroup" = 4
    )) +
    labs(
      title = title,
      x     = sprintf("MDS1 (%.1f%% variance)", pct[1]),
      y     = sprintf("MDS2 (%.1f%% variance)", pct[2]),
      color = "", shape = ""
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

# Plot 1: all populations
r1 <- make_mds_df(dmat)
p1 <- make_plot(r1$df, r1$pct,
                sprintf("MDS — all populations (%d pops, f2 distances) — target: %s",
                        length(MDS_POPS), TARGET))

# Plot 2: regional focus (drop distant global outgroups to zoom in)
close_pops <- pops_all[!pops_all %in% c("Mbuti", "Yoruba", "Papuan",
                                         "Han", "China_TianyuanCave_UP",
                                         "Ethiopia_MotaCave_4500BP")]
r2 <- make_mds_df(dmat[close_pops, close_pops])
p2 <- make_plot(r2$df, r2$pct,
                sprintf("MDS — regional focus (%d pops) — target: %s",
                        length(close_pops), TARGET))

outdir <- dirname(F2_DIR)
ggsave(file.path(outdir, "mds_all.pdf"),  p1, width = 16, height = 11)
ggsave(file.path(outdir, "mds_euro.pdf"), p2, width = 16, height = 11)

cat("Plots saved:\n")
cat(" ", file.path(outdir, "mds_all.pdf"),  "— all populations\n")
cat(" ", file.path(outdir, "mds_euro.pdf"), "— regional focus\n")
cat("\nDone.\n")
