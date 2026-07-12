library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

OUTGROUPS <- SLAVIC_OUTGROUPS

SOURCES <- SLAVIC_BEST_SOURCES

FOCUSED_F2_DIR <- sub("f2/?$", "f2_focused/", F2_DIR)

pops_needed <- unique(c(TARGET, SOURCES, OUTGROUPS))
f2 <- f2_from_precomp(FOCUSED_F2_DIR, pops = pops_needed)

res <- qpadm(f2, target = TARGET, left = SOURCES, right = OUTGROUPS)

w_col  <- intersect(c("weight", "w", "est"),        names(res$weights))[1]
se_col <- intersect(c("se", "se_weight"),            names(res$weights))[1]

full_p <- res$popdrop %>% filter(!grepl("1", pat)) %>% slice(1) %>% pull(p)

label_map <- c(
  "Poland_EarlyMedieval_Slav" = "Slavic migrants\n(Poland EM, 6–9th c. CE)",
  "Croatia_EIA"               = "Iron Age Balkans\n(Croatia EIA, ~800–400 BCE)"
)
color_map <- c(
  "Poland_EarlyMedieval_Slav" = "#4472C4",
  "Croatia_EIA"               = "#C0504D"
)

df <- res$weights %>%
  mutate(
    pct       = .data[[w_col]]  * 100,
    se_pct    = .data[[se_col]] * 100,
    label     = label_map[left],
    slice_lbl = sprintf("%s\n%.0f%% ± %.0f%%\n(z = %.2f)", label, pct, se_pct, pct / se_pct)
  )

p <- ggplot(df, aes(x = "", y = pct, fill = left)) +
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  coord_polar(theta = "y", start = 0) +
  geom_text(
    aes(label = slice_lbl),
    position  = position_stack(vjust = 0.5),
    size      = 4.2,
    color     = "white",
    fontface  = "bold",
    lineheight = 1.3
  ) +
  scale_fill_manual(values = color_map) +
  labs(
    title    = sprintf("Slavic ancestry model — %s", TARGET),
    subtitle = sprintf("qpAdm  |  p = %.3f  |  %s", full_p, paste(SOURCES, collapse = " + "))
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40", margin = margin(b = 8)),
    legend.position = "none",
    plot.margin   = margin(12, 12, 12, 12)
  )

outfile <- file.path(dirname(MERGED_PREFIX), "slavic_pie.pdf")
ggsave(outfile, p, width = 6, height = 6)
cat("Saved:", outfile, "\n\n")

cat(sprintf("Model: %s + %s\n", SOURCES[1], SOURCES[2]))
cat(sprintf("p = %.4f  (%s)\n", full_p, if (full_p > 0.05) "PASS" else "FAIL"))
for (i in seq_len(nrow(df))) {
  cat(sprintf("  %-46s  %+.1f%% ± %.1f%%  (z=%.2f)\n",
              df$left[i], df$pct[i], df$se_pct[i], df$pct[i] / df$se_pct[i]))
}
