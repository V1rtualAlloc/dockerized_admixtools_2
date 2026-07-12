# =============================================================================
# Shared configuration — sourced by all analysis scripts
# Edit TARGET and BACKGROUND for each new sample
# =============================================================================

TARGET      <- "me"          # must match the label used during PLINK conversion
BACKGROUND  <- "european"    # one of: european, south_asian, east_asian,
#                              middle_eastern, central_asian, native_american

MERGED_PREFIX <- "/data/me/merged"   # path to merged.geno/.snp/.ind (no extension)
F2_DIR        <- "/data/me/f2/"      # where to cache f2 blocks (created if needed)


# Broad modern/high-diversity populations used as extra context in the
# exploratory scripts (f3_outgroup.R, f3_admixture.R, pca.R), regardless of
# BACKGROUND — these aren't claimed to be closely related to any particular
# target, just useful fixed points for orientation.
WORLD_REFS <- c("Mbuti", "Yoruba", "Han", "Papuan",
                "French", "Sardinian", "Spanish", "Russian")

# Shared outgroups for all Slavic model scripts
SLAVIC_OUTGROUPS <- c(
  "Mbuti", "Yoruba", "Han", "Papuan",
  "Ethiopia_MotaCave_4500BP", "China_TianyuanCave_UP",
  "Iran_GanjDareh_N", "Israel_Natufian", "Turkey_N"
)

# Modern/historical populations shown alongside TARGET for comparison in
# qpadm.R, f3_outgroup.R, f3_admixture.R, f4.R and pca.R — populations from
# the same region as the target, so its results can be read in context.
# Keyed by BACKGROUND because "relevant comparison population" is inherently
# regional. Only "european" is populated (this project's original Balkan/
# Slavic use case). If you're running a different BACKGROUND, add your own
# region's comparison populations here — verify names against the AADR
# .ind file first (see CLAUDE.md's f2 cache pollution warning for N<10 pops).
REFERENCES_BY_BACKGROUND <- list(
  european        = c(
    "Serbia_Medieval",
    "Serbia_EarlyMedieval_Byzantine_Slav",
    "BosniaHerzegovina_Medieval",
    "Croatia_Medieval_Modern",
    "NorthMacedonia_Medieval",
    "Bulgarian",
    "Albanian",
    "Hungarian",
    "Romania_Medieval",
    "Greek",
    "Russian"
  ),
  south_asian     = character(0),
  east_asian      = character(0),
  middle_eastern  = character(0),
  central_asian   = character(0),
  native_american = character(0)
)
REFERENCES <- REFERENCES_BY_BACKGROUND[[BACKGROUND]]


# Ancestry models: sources (qpAdm left) and outgroups (right) per background.
# `pool` is an optional wider set of ancient-population candidates used by the
# exploratory scripts (qpadm_rotate.R, f3_admixture.R); scripts fall back to
# `sources` when a background has no dedicated pool.
# Outgroup rules for european:
#   - Excluded: EHG and CHG (direct components of Yamnaya)
#   - Included: Iran_N and Natufian (Near Eastern resolution without being sources)
MODELS <- list(

  european = list(
    sources  = c("Russia_Samara_EBA_Yamnaya",
                 "Serbia_IronGates_Mesolithic",
                 "Turkey_N"),
    outgroup = c("Mbuti", "Yoruba", "Han", "Papuan",
                 "Ethiopia_MotaCave_4500BP",
                 "China_TianyuanCave_UP",
                 "Iran_GanjDareh_N",
                 "Israel_Natufian"),
    pool     = c("Russia_Samara_EBA_Yamnaya",
                 "Serbia_IronGates_Mesolithic",
                 "Luxembourg_Loschbour_Mesolithic",
                 "Turkey_N",
                 "Russia_Karelia_Mesolithic_HG",
                 "Georgia_KotiasKlde_Mesolithic")
  ),

  south_asian = list(
    sources  = c("Russia_Samara_EBA_Yamnaya",
                 "Iran_GanjDareh_N",
                 "Papuan"),
    outgroup = c("Mbuti", "Yoruba", "Han",
                 "Turkey_N",
                 "Luxembourg_Loschbour_Mesolithic",
                 "Russia_Karelia_Mesolithic_HG",
                 "Israel_Natufian",
                 "China_TianyuanCave_UP")
  ),

  east_asian = list(
    sources  = c("China_N",
                 "Russia_PrimorskyKrai_MN_Boisman",
                 "China_TianyuanCave_UP"),
    outgroup = c("Mbuti", "Yoruba",
                 "Turkey_N",
                 "Luxembourg_Loschbour_Mesolithic",
                 "Russia_Samara_EBA_Yamnaya",
                 "Iran_GanjDareh_N",
                 "Israel_Natufian")
  ),

  middle_eastern = list(
    sources  = c("Israel_Natufian",
                 "Iran_GanjDareh_N",
                 "Georgia_KotiasKlde_Mesolithic"),
    outgroup = c("Mbuti", "Yoruba", "Han", "Papuan",
                 "Russia_Samara_EBA_Yamnaya",
                 "Luxembourg_Loschbour_Mesolithic",
                 "Turkey_N",
                 "China_TianyuanCave_UP")
  ),

  central_asian = list(
    sources  = c("Russia_Samara_EBA_Yamnaya",
                 "Iran_GanjDareh_N",
                 "China_N"),
    outgroup = c("Mbuti", "Yoruba", "Papuan",
                 "Turkey_N",
                 "Luxembourg_Loschbour_Mesolithic",
                 "Russia_Karelia_Mesolithic_HG",
                 "Israel_Natufian",
                 "China_TianyuanCave_UP")
  ),

  native_american = list(
    sources  = c("China_TianyuanCave_UP",
                 "Russia_PrimorskyKrai_MN_Boisman",
                 "Russia_Karelia_Mesolithic_HG"),
    outgroup = c("Mbuti", "Yoruba", "Han",
                 "Turkey_N",
                 "Luxembourg_Loschbour_Mesolithic",
                 "Russia_Samara_EBA_Yamnaya",
                 "Iran_GanjDareh_N")
  )
)


# =============================================================================
# Slavic ancestry model — shared population pools (case-study, european
# background only; see CLAUDE.md's case-study note). Consumed by
# scripts/slavic_*.R and report.R. Previously these were duplicated ad hoc
# across each script (with drift between them) — centralized here so a pool
# only needs editing in one place. Where scripts intentionally used different
# scope (e.g. a narrower pool to keep a pooled cache's SEs tight), that's kept
# as a separate named variant rather than silently unified.
# =============================================================================

# Slavic migrant proxies (identical across slavic_model.R and slavic_focused.R)
SLAVIC_MIGRANT_POOL <- c("Poland_EarlyMedieval_Slav", "Poland_EarlySlav")

# Best-supported 2-source model (slavic_pie.R, slavic_outgroup_test.R)
SLAVIC_BEST_SOURCES <- c("Poland_EarlyMedieval_Slav", "Croatia_EIA")

# Balkan Iron Age substrate candidates
#   CORE — slavic_model.R (rotating search), slavic_pooled.R (pooling)
#   FULL — slavic_focused.R; adds Serbia_ImperialRoman as a 4th candidate
SLAVIC_BALKAN_IA_CORE <- c("Croatia_EIA", "Bulgaria_KapitanAndreevo_EIA", "NorthMacedonia_IA")
SLAVIC_BALKAN_IA_FULL <- c(SLAVIC_BALKAN_IA_CORE, "Serbia_ImperialRoman")

# Roman/Byzantine import candidates
#   FULL    — slavic_model.R (rotating search): all era-matched variants
#   FOCUSED — slavic_focused.R: narrower set
#   MINIMAL — slavic_pooled.R: narrowest, keeps the pooled cache's SEs tight
SLAVIC_ROMAN_FULL <- c(
  "Italy_Lazio_ImperialRoman_Roman",
  "Italy_Lazio_LateAntiquity_ImperialRoman_Roman",
  "Italy_Lazio_LateAntiquity_Roman",
  "Turkey_Medieval_Byzantine",
  "Turkey_EarlyMedieval_Byzantine",
  "Turkey_LateAntiquity_Byzantine",
  "Turkey_LateAntiquity_ImperialRoman"
)
SLAVIC_ROMAN_FOCUSED <- c(
  "Turkey_Medieval_Byzantine", "Turkey_LateAntiquity_Byzantine",
  "Turkey_EarlyMedieval_Byzantine", "Italy_Lazio_ImperialRoman_Roman"
)
SLAVIC_ROMAN_MINIMAL <- c("Italy_Lazio_ImperialRoman_Roman", "Turkey_Medieval_Byzantine")

# Modern reference populations tested in slavic_modern.R
SLAVIC_MODERN_REFS <- c("Bulgarian", "Albanian", "Greek", "Greek_1", "Greek_Crete",
                        "Hungarian", "Czech", "Polish", "Russian")

# Fixed named models for the report.R Slavic comparison table
SLAVIC_NAMED_MODELS <- list(
  list(id="A", sources=c("Poland_EarlyMedieval_Slav","Serbia_LateAntiquity_ImperialRoman"),
       label="A -- Slavic + Roman Balkans [bundled, for comparison]"),
  list(id="B", sources=c("Poland_EarlyMedieval_Slav","Croatia_EIA"),
       label="B -- Slavic + Iron Age Balkans"),
  list(id="C", sources=c("Poland_EarlyMedieval_Slav","Croatia_EIA","Turkey_Medieval_Byzantine"),
       label="C -- Slavic + IA Balkans + Byzantine Anatolia"),
  list(id="D", sources=c("Poland_EarlyMedieval_Slav","Croatia_EIA","Israel_Phoenician"),
       label="D -- Slavic + IA Balkans + Levantine"),
  list(id="E", sources=c("Poland_EarlyMedieval_Slav","Croatia_EIA","Italy_Lazio_ImperialRoman_Roman"),
       label="E -- Slavic + IA Balkans + Roman Italian")
)
