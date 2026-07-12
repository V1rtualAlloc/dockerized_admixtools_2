library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

OUTFILE <- file.path(dirname(F2_DIR), "ancestry_report.pdf")
CHR_CSV <- file.path(dirname(F2_DIR), "chr_painting.csv")

cat("=== Generating ancestry report ===\n")
cat("Output:", OUTFILE, "\n\n")

model   <- MODELS[[BACKGROUND]]
sources <- model$sources
right   <- model$outgroup
labels  <- c("Steppe", "Balkan HG", "Anatolian")

# ── Load f2 ──────────────────────────────────────────────────────────────────
all_pops <- unique(c(TARGET, sources, right, REFERENCES,
                     unlist(lapply(MODELS, `[[`, "sources"))))
all_pops <- intersect(all_pops, list.dirs(F2_DIR, full.names=FALSE, recursive=FALSE))
f2 <- f2_from_precomp(F2_DIR)

# ── Dataset stats (computed from the actual merged files, not hardcoded) ─────
cat("Reading merged dataset stats...\n")
n_snps   <- nrow(data.table::fread(paste0(MERGED_PREFIX, ".snp"), header=FALSE))
ind_tab  <- read.table(paste0(MERGED_PREFIX, ".ind"), header=FALSE,
                       col.names=c("id","sex","pop"), stringsAsFactors=FALSE)
n_ind    <- nrow(ind_tab)
n_target_inds  <- sum(ind_tab$pop == TARGET)
n_ancient_inds <- n_ind - n_target_inds

# ── Helpers ──────────────────────────────────────────────────────────────────
FILL <- c(Steppe="#E69F00", WHG="#56B4E9", "Balkan HG"="#56B4E9", Anatolian="#009E73")
BASE <- theme_bw(base_size=11) + theme(plot.title=element_text(face="bold", size=13))

voidtext <- function(lines, title="", mono=TRUE) {
  txt <- paste(lines, collapse="\n")
  fam <- if (mono) "mono" else "sans"
  ggplot() +
    annotate("text", x=0, y=0.95, label=txt, hjust=0, vjust=1,
             family=fam, size=3.4, lineheight=1.4) +
    labs(title=title) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    theme_void() +
    theme(plot.title=element_text(face="bold", size=14, margin=margin(b=8)),
          plot.margin=margin(15,15,15,15))
}

# ── Page 2 data: main qpAdm model ────────────────────────────────────────────
cat("Running qpAdm (main model)...\n")

run_qpadm <- function(ref) {
  pops <- c(TARGET, sources, right)
  f2s  <- tryCatch(f2_from_precomp(F2_DIR, pops=pops), error=function(e) NULL)
  if (is.null(f2s)) return(NULL)
  r <- tryCatch(qpadm(f2s, target=TARGET, left=sources, right=right),
                error=function(e) NULL)
  if (is.null(r) || is.null(r$weights)) return(NULL)
  full  <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  w_col <- intersect(c("weight","w","est"), names(r$weights))
  se_col <- intersect(c("se","std.error"), names(r$weights))
  wts <- setNames(r$weights[[w_col[1]]], r$weights$left)
  ses <- if (length(se_col)) setNames(r$weights[[se_col[1]]], r$weights$left) else wts*0
  tibble(reference=ref, p=full$p[1],
         source=sources, label=labels,
         weight=as.numeric(wts), se=as.numeric(ses))
}

main_res <- run_qpadm("Primary")

# ── Page 3 data: qpAdm rotating (computed early so the title page can cite it) ──
cat("Running qpAdm rotating...\n")

SOURCE_POOL <- c("Russia_Samara_EBA_Yamnaya", "Serbia_IronGates_Mesolithic",
                 "Luxembourg_Loschbour_Mesolithic", "Turkey_N",
                 "Russia_Karelia_Mesolithic_HG", "Georgia_KotiasKlde_Mesolithic")
SOURCE_POOL <- intersect(SOURCE_POOL,
                         list.dirs(F2_DIR, full.names=FALSE, recursive=FALSE))

combos <- c(
  combn(SOURCE_POOL, 2, simplify=FALSE),
  combn(SOURCE_POOL, 3, simplify=FALSE)
)

run_one_rot <- function(srcs) {
  pops <- c(TARGET, srcs, right)
  f2s  <- tryCatch(f2_from_precomp(F2_DIR, pops=pops), error=function(e) NULL)
  if (is.null(f2s)) return(NULL)
  r <- tryCatch(qpadm(f2s, target=TARGET, left=srcs, right=right),
                error=function(e) NULL)
  if (is.null(r) || is.null(r$weights)) return(NULL)
  full  <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  if (nrow(full)==0) return(NULL)
  w_col <- intersect(c("weight","w","est"), names(r$weights))
  wts   <- setNames(r$weights[[w_col[1]]], r$weights$left)
  feasible <- all(wts >= -0.01 & wts <= 1.01)
  w_str    <- paste(sprintf("%s %.0f%%", names(wts), wts*100), collapse="  ")
  tibble(n=length(srcs), model=paste(srcs, collapse=" + "),
         p=full$p, feasible=feasible, weights=w_str)
}

rot_res <- map_dfr(combos, run_one_rot) %>% arrange(desc(p))

pass <- filter(rot_res, p > 0.05, feasible)
fail <- filter(rot_res, p <= 0.05 | !feasible)

# ── Page 1: Title (built after qpAdm + rotate so it can cite real numbers) ───
build_title_page <- function(main_res) {
  if (!is.null(main_res) && nrow(main_res) == 3) {
    w  <- setNames(main_res$weight, main_res$label)
    se <- setNames(main_res$se,     main_res$label)
    z  <- w / se
    fmt <- function(lbl)
      sprintf("%3.0f%% +/- %2.0f%%  (z = %.1f)", w[lbl]*100, se[lbl]*100, z[lbl])
    anc_lines <- c(
      sprintf("    Steppe (Yamnaya, Russia ~3000 BCE)         %s", fmt("Steppe")),
      sprintf("    Balkan HG (Iron Gates, Serbia ~8000 BCE)   %s", fmt("Balkan HG")),
      sprintf("    Anatolian EEF (Turkey Neolithic ~7000 BCE) %s", fmt("Anatolian"))
    )
    hg_z <- z["Balkan HG"]
    z_note <- if (!is.na(hg_z) && hg_z < 2) c(
      "",
      "  PRECISION NOTE (Balkan HG)",
      sprintf("  z = %.1f means this estimate is close to zero -- compatible", hg_z),
      "  with a wide range at 2-sigma. This does not mean the component is",
      "  absent: the smallest of the three sources is hardest to pin down from",
      "  a single genome, since WHG/EHG-like populations are genetically closer",
      "  to each other than to the other two sources, leaving less statistical",
      "  leverage. Multiple related individuals would tighten the estimate;",
      "  N=1 cannot."
    ) else NULL
  } else {
    anc_lines <- c("    qpAdm model unavailable -- see page 2 for details.")
    z_note <- NULL
  }

  rot_note <- if (nrow(rot_res) > 0) c(
    sprintf("  Of %d alternative 2-/3-source combinations tested from the pool", nrow(rot_res)),
    sprintf("  (%s), %d pass (p > 0.05, feasible weights). Full list on page 3.",
            paste(SOURCE_POOL, collapse=", "), nrow(pass))
  ) else NULL

  c(
    "",
    sprintf("Target individual:  %s", TARGET),
    sprintf("Reference panel:    AADR v66.p1 (1240k SNP panel)"),
    sprintf("Analysis date:      %s", Sys.Date()),
    sprintf("SNPs (merged):      %s  |  Ancient individuals: %s",
            format(n_snps, big.mark=","), format(n_ancient_inds, big.mark=",")),
    "",
    strrep("-", 52),
    "",
    "  GENOME-WIDE ANCESTRY (qpAdm, 3-source European model)",
    "",
    anc_lines,
    "",
    strrep("-", 52),
    "",
    "  ALTERNATIVE MODELS",
    "",
    rot_note,
    z_note,
    "",
    strrep("-", 52),
    "",
    "  HUNTER-GATHERER PROXY",
    "",
    sprintf("  The Balkan HG component uses %s as its proxy population.", sources[2]),
    "  Alternative HG proxies can be compared with iron_gates_test.R.",
    ""
  )
}

cat("Building title page...\n")
title_lines <- build_title_page(main_res)
p_title <- ggplot() +
  annotate("text", x=0.5, y=0.82, label="Ancient DNA Ancestry Report",
           hjust=0.5, size=7, fontface="bold") +
  annotate("text", x=0.5, y=0.72, label=TARGET,
           hjust=0.5, size=5, color="grey40") +
  annotate("text", x=0.05, y=0.58,
           label=paste(Filter(Negate(is.null), tail(title_lines, -1)), collapse="\n"),
           hjust=0, vjust=1, family="mono", size=3.4, lineheight=1.45) +
  theme_void() +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))

# Also run across reference populations
ref_res <- map_dfr(names(REFERENCES), function(ref) {
  pops <- c(TARGET, sources, right)
  tryCatch({
    f2s <- f2_from_precomp(F2_DIR, pops=pops)
    r   <- qpadm(f2s, target=TARGET, left=sources, right=right)
    if (is.null(r$weights)) return(NULL)
    full  <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
    w_col <- intersect(c("weight","w","est"), names(r$weights))
    wts   <- setNames(r$weights[[w_col[1]]], r$weights$left)
    tibble(reference=ref, p=full$p[1],
           source=sources, label=labels, weight=as.numeric(wts), se=0)
  }, error=function(e) NULL)
})

qpadm_dat <- bind_rows(main_res, ref_res) %>%
  mutate(label=factor(label, levels=c("Steppe","Balkan HG","Anatolian")),
         reference=factor(reference, levels=unique(reference)))

# Error bars at top of each stacked segment for primary target only
primary_err <- qpadm_dat %>%
  filter(reference == "Primary") %>%
  arrange(label) %>%
  mutate(ytop  = cumsum(weight) * 100,
         ymin  = (ytop/100 - se) * 100,
         ymax  = (ytop/100 + se) * 100)

p_qpadm <- ggplot(qpadm_dat, aes(x=reference, y=weight*100, fill=label)) +
  geom_col(position="stack", width=0.6) +
  geom_errorbar(data=primary_err,
                aes(x=reference, y=ytop, ymin=ymin, ymax=ymax),
                inherit.aes=FALSE,
                width=0.15, linewidth=0.6, color="grey20") +
  scale_fill_manual(values=FILL) +
  scale_y_continuous(labels=function(x) paste0(x,"%"), limits=c(0,115)) +
  labs(title="qpAdm: 3-source ancestry model",
       subtitle=paste("Sources:", paste(sources, collapse=" + "),
                      "  |  Error bars (+/-1 SE) at segment tops, target only"),
       x="Reference population", y="Ancestry proportion", fill="Component") +
  BASE + theme(legend.position="bottom", axis.text.x=element_text(angle=30, hjust=1))

# ── Page 3: qpAdm rotating (text table, using rot_res computed above) ───────
rot_lines <- c(
  sprintf("%-3s  %-6s  %-8s  %s", "n", "p", "feasible", "weights / model"),
  strrep("-", 78),
  sprintf("  PASSING models (p > 0.05, feasible weights)  [%d of %d]",
          nrow(pass), nrow(rot_res)),
  strrep("-", 78)
)
if (nrow(pass) > 0) {
  rot_lines <- c(rot_lines,
    sprintf("%-3d  %-6.4f  %-8s  %s", pass$n, pass$p,
            ifelse(pass$feasible,"yes","no"), pass$weights))
}
rot_lines <- c(rot_lines, strrep("-", 78),
               sprintf("  FAILING models  [%d]", nrow(fail)), strrep("-", 78))
if (nrow(fail) > 0) {
  rot_lines <- c(rot_lines,
    sprintf("%-3d  %-6.4f  %-8s  %s", fail$n, fail$p,
            ifelse(fail$feasible,"yes","no"), fail$model))
}

p_rotate <- voidtext(rot_lines, title="qpAdm rotating: all 2- and 3-source combinations")

# ── Page 4: f3 outgroup ───────────────────────────────────────────────────────
cat("Computing outgroup f3...\n")

all_pops_in_cache <- list.dirs(F2_DIR, full.names=FALSE, recursive=FALSE)
pop_candidates <- setdiff(all_pops_in_cache, c(TARGET, right))

f2_full <- f2_from_precomp(F2_DIR)
outgroup_ref <- right[1]

f3_res <- tryCatch({
  f3(f2_full, pop1=TARGET, pop2=pop_candidates, pop3=outgroup_ref) %>%
    arrange(desc(est)) %>%
    slice_head(n=20)
}, error=function(e) NULL)

if (!is.null(f3_res)) {
  f3_dat <- f3_res %>%
    mutate(pop2 = factor(pop2, levels=rev(pop2)))
  p_f3 <- ggplot(f3_dat, aes(x=est, y=pop2)) +
    geom_col(fill="steelblue", width=0.7) +
    geom_errorbarh(aes(xmin=est-2*se, xmax=est+2*se), height=0.3, linewidth=0.5) +
    labs(title=sprintf("Outgroup f3: populations most genetically similar to %s", TARGET),
         subtitle=sprintf("f3(%s, X; %s)  --  higher = more shared drift", TARGET, outgroup_ref),
         x="f3 statistic", y=NULL) +
    BASE
} else {
  p_f3 <- voidtext("f3 outgroup computation failed.", title="Outgroup f3")
}

# ── Page 5: Chromosome painting ───────────────────────────────────────────────
cat("Loading chromosome painting results...\n")

if (file.exists(CHR_CSV)) {
  # chr_painting.R (run against the same config.R) names its weight columns
  # after the actual source populations, e.g. sources[1] == "Russia_Samara_EBA_Yamnaya".
  chr_res <- read_csv(CHR_CSV, show_col_types=FALSE)
  chr_dat <- chr_res %>%
    transmute(chr,
              Steppe    = .data[[sources[1]]],
              WHG       = .data[[sources[2]]],
              Anatolian = .data[[sources[3]]]) %>%
    pivot_longer(-chr, names_to="ancestry", values_to="prop") %>%
    mutate(ancestry=factor(ancestry, levels=c("Steppe","WHG","Anatolian")),
           chr=factor(chr, levels=1:22))

  chr_means <- chr_res %>%
    summarise(Steppe=weighted.mean(.data[[sources[1]]],n_snps),
              WHG=weighted.mean(.data[[sources[2]]],n_snps),
              Anatolian=weighted.mean(.data[[sources[3]]],n_snps)) %>%
    pivot_longer(everything(), names_to="ancestry", values_to="mean") %>%
    mutate(ancestry=factor(ancestry, levels=c("Steppe","WHG","Anatolian")))

  p_chr <- ggplot(chr_dat, aes(chr, prop*100, fill=ancestry)) +
    geom_col(width=0.8) +
    geom_hline(data=chr_means, aes(yintercept=mean*100, color=ancestry),
               linetype="dashed", linewidth=0.6, show.legend=FALSE) +
    scale_fill_manual(values=FILL) +
    scale_color_manual(values=c(Steppe="#B87D00", WHG="#2E7CA8", Anatolian="#007050")) +
    scale_y_continuous(limits=c(0,105), labels=function(x) paste0(x,"%")) +
    labs(title="Chromosome painting: ancestry proportion per autosome",
         subtitle="Method: diploid binomial MLE on allele frequencies  |  Dashed = genome-wide mean",
         x="Chromosome", y="Ancestry proportion", fill="Source") +
    BASE + theme(legend.position="bottom", panel.grid.major.x=element_blank())
} else {
  p_chr <- voidtext(
    c("chr_painting.csv not found.",
      "Run chr_painting.R first, then re-run report.R."),
    title="Chromosome painting")
}

# ── Pages 6–7: Slavic model ───────────────────────────────────────────────────
cat("Running Slavic model...\n")

SLAVIC_F2_DIR <- sub("f2/?$", "f2_slavic/", F2_DIR)
SLAVIC_OG <- c("Mbuti","Yoruba","Han","Papuan",
                "Ethiopia_MotaCave_4500BP","China_TianyuanCave_UP",
                "Iran_GanjDareh_N","Israel_Natufian","Turkey_N")

SLAVIC_MODELS <- SLAVIC_NAMED_MODELS

run_slavic_model <- function(f2s, sources) {
  r <- tryCatch(
    qpadm(f2s, target=TARGET, left=sources, right=SLAVIC_OG),
    error=function(e) NULL
  )
  if (is.null(r) || is.null(r$weights)) return(NULL)
  full  <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  w_col  <- intersect(c("weight","w","est"),      names(r$weights))
  se_col <- intersect(c("se","se_weight"),         names(r$weights))
  z_col  <- intersect(c("z","z_weight"),           names(r$weights))
  wts <- setNames(r$weights[[w_col[1]]],  r$weights$left)
  ses <- if (length(se_col)) setNames(r$weights[[se_col[1]]], r$weights$left) else rep(NA_real_, length(wts))
  zs  <- if (length(z_col))  setNames(r$weights[[z_col[1]]],  r$weights$left) else wts / ses
  list(p=full$p[1], wts=wts, ses=ses, zs=zs,
       feasible=all(wts >= -0.01 & wts <= 1.01))
}

slavic_cache_ok <- length(list.files(SLAVIC_F2_DIR, pattern="\\.rds$")) > 0

if (slavic_cache_ok) {
  slavic_pops <- unique(c(TARGET, unlist(lapply(SLAVIC_MODELS, `[[`, "sources")), SLAVIC_OG))
  f2_slav <- tryCatch(f2_from_precomp(SLAVIC_F2_DIR, pops=slavic_pops), error=function(e) NULL)
} else {
  f2_slav <- NULL
}

if (!is.null(f2_slav)) {
  slavic_results <- lapply(SLAVIC_MODELS, function(m) {
    res <- run_slavic_model(f2_slav, m$sources)
    if (!is.null(res)) res$id <- m$id
    res
  })
  names(slavic_results) <- sapply(SLAVIC_MODELS, `[[`, "id")

  # Determine the best-supported model: highest p among feasible fits.
  model_p <- sapply(slavic_results, function(r) if (!is.null(r) && r$feasible) r$p else NA_real_)
  best_id <- if (all(is.na(model_p))) NA_character_ else names(model_p)[which.max(model_p)]

  # ── Page 6: text table ──
  hdr  <- sprintf("  %-2s  %-6s  %-10s  %-50s  %s", "M", "p", "weights", "source", "weight ± se  (z)")
  sep  <- strrep("-", 100)
  tbl_lines <- c(
    "Models tested (outgroups: Turkey_N, Israel_Natufian, Iran_GanjDareh_N, + Africa/EastAsia)",
    "",
    hdr, sep
  )
  for (i in seq_along(SLAVIC_MODELS)) {
    m   <- SLAVIC_MODELS[[i]]
    res <- slavic_results[[m$id]]
    best_marker <- if (!is.na(best_id) && m$id == best_id) "  *** best ***" else ""
    if (is.null(res)) {
      tbl_lines <- c(tbl_lines, sprintf("  %s  [error]  %s", m$id, m$label))
      next
    }
    feas <- if (res$feasible) "feasible  " else "INFEASIBLE"
    tbl_lines <- c(tbl_lines,
      sprintf("  %s   %.4f  %s  %s%s", m$id, res$p, feas, m$label, best_marker)
    )
    for (nm in names(res$wts)) {
      tbl_lines <- c(tbl_lines,
        sprintf("  %-4s%-6s%-10s  %-50s  %+5.1f%% ± %4.1f%%  (z=%+.2f)",
                "", "", "", nm,
                res$wts[nm]*100, res$ses[nm]*100, res$zs[nm]))
    }
    tbl_lines <- c(tbl_lines, "")
  }

  interp_lines <- if (!is.na(best_id)) {
    bm  <- slavic_results[[best_id]]
    lbl <- SLAVIC_MODELS[[which(sapply(SLAVIC_MODELS, `[[`, "id") == best_id)]]$label
    wt_str <- paste(sprintf("%s=%.0f%%", names(bm$wts), bm$wts*100), collapse=" + ")
    n_infeasible <- sum(sapply(slavic_results, function(r) !is.null(r) && !r$feasible))
    c(
      "Interpretation:",
      sprintf("  Best-supported model: %s (p=%.3f)", lbl, bm$p),
      sprintf("  %s", wt_str),
      if (n_infeasible > 0) c("",
        sprintf("  %d of %d tested models have out-of-range (infeasible) weights --",
                n_infeasible, length(slavic_results)),
        "  the extra source could not be statistically separated from the others",
        "  with this data. This does not mean those models are wrong, just",
        "  unresolved by the available populations.") else NULL
    )
  } else {
    c("Interpretation:", "  No feasible passing model found among those tested.")
  }
  tbl_lines <- c(tbl_lines, sep, "", interp_lines)
  p_slavic_text <- voidtext(tbl_lines, title="Slavic / pre-Slavic ancestry model")

  # ── Page 7: bar chart (Models A, B, E) ──
  res_A <- slavic_results[["A"]]
  res_B <- slavic_results[["B"]]
  res_E <- slavic_results[["E"]]
  SHORT <- c(
    "Poland_EarlyMedieval_Slav"          = "Slavic (Poland EM)",
    "Serbia_LateAntiquity_ImperialRoman" = "Roman Balkans",
    "Croatia_EIA"                        = "Iron Age Balkans",
    "Turkey_Medieval_Byzantine"          = "Byzantine Anatolia",
    "Israel_Phoenician"                  = "Levantine",
    "Italy_Lazio_ImperialRoman_Roman"    = "Roman Italian"
  )
  SCOL <- c(
    "Slavic (Poland EM)"    = "#4C72B0",
    "Iron Age Balkans"      = "#DD8452",
    "Roman Balkans"         = "#937860",
    "Byzantine Anatolia"    = "#C44E52",
    "Levantine"             = "#8172B2",
    "Roman Italian"         = "#CCB974"
  )

  make_bar_dat <- function(res, model_id) {
    if (is.null(res)) return(NULL)
    tibble(
      model  = model_id,
      source = SHORT[names(res$wts)],
      weight = as.numeric(res$wts),
      se     = as.numeric(res$ses)
    )
  }

  bar_dat <- bind_rows(
    make_bar_dat(res_A, "A: Slavic +\nRoman Balkans"),
    make_bar_dat(res_B, "B: Slavic +\nIron Age Balkans"),
    make_bar_dat(res_E, "E: Slavic + Iron Age\nBalkans + Roman Italian")
  ) %>%
    mutate(source = factor(source, levels=names(SCOL)),
           model  = factor(model, levels=unique(model)))

  # error bars: position at top of each segment
  err_dat <- bar_dat %>%
    group_by(model) %>%
    arrange(model, source) %>%
    mutate(ytop = cumsum(weight) * 100,
           ymin = (ytop/100 - se) * 100,
           ymax = (ytop/100 + se) * 100)

  best_label <- if (!is.na(best_id)) sprintf("Best model: %s (p=%.2f)", best_id, model_p[best_id]) else "No passing model"

  p_slavic_chart <- ggplot(bar_dat, aes(x=model, y=weight*100, fill=source)) +
    geom_col(position="stack", width=0.55) +
    geom_errorbar(data=err_dat,
                  aes(x=model, y=ytop, ymin=ymin, ymax=ymax),
                  inherit.aes=FALSE, width=0.12, linewidth=0.6, color="grey20") +
    scale_fill_manual(values=SCOL, drop=FALSE) +
    scale_y_continuous(labels=function(x) paste0(x,"%"), limits=c(0,115)) +
    labs(title="Slavic model: ancestry proportions",
         subtitle=paste0("Target: ", TARGET, "  |  ", best_label,
                         "  |  Error bars = ±1 SE"),
         x=NULL, y="Ancestry proportion", fill="Source") +
    BASE + theme(legend.position="bottom",
                 axis.text.x=element_text(size=11))

} else {
  msg <- c("Slavic model cache not found.",
           sprintf("Run slavic_model.R first to build %s, then re-run report.R.", SLAVIC_F2_DIR))
  p_slavic_text  <- voidtext(msg, title="Slavic / pre-Slavic ancestry model")
  p_slavic_chart <- p_slavic_text
}

# ── Render ────────────────────────────────────────────────────────────────────
cat("Writing PDF...\n")

pdf(OUTFILE, width=11, height=8.5, paper="special")
print(p_title)
print(p_qpadm)
print(p_rotate)
print(p_f3)
print(p_chr)
print(p_slavic_text)
print(p_slavic_chart)
invisible(dev.off())

cat("Report saved:", OUTFILE, "\nDone.\n")
