# =============================================================================
# Shared configuration — sourced by all analysis scripts
# Edit TARGET and BACKGROUND for each new sample
# =============================================================================

TARGET      <- "me"          # must match the label used during PLINK conversion
BACKGROUND  <- "european"    # one of: european, south_asian, east_asian,
#                              middle_eastern, central_asian, native_american

MERGED_PREFIX <- "/data/me/merged"   # path to merged.geno/.snp/.ind (no extension)
F2_DIR        <- "/data/me/f2/"      # where to cache f2 blocks (created if needed)


# Reference populations included in every analysis for comparison
# Shared outgroups for all Slavic model scripts
SLAVIC_OUTGROUPS <- c(
  "Mbuti", "Yoruba", "Han", "Papuan",
  "Ethiopia_MotaCave_4500BP", "China_TianyuanCave_UP",
  "Iran_GanjDareh_N", "Israel_Natufian", "Turkey_N"
)

REFERENCES <- c(
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
)


# Ancestry models: sources (qpAdm left) and outgroups (right) per background.
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
                 "Israel_Natufian")
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
