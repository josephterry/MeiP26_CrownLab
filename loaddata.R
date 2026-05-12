suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(glmmTMB); library(blme);  library(lme4);  library(lmerTest);
  library(emmeans); library(performance); library(car); library(MuMIn); library(DHARMa);
  library(readr); library(pROC); library(MASS); library(scales); library(patchwork)
})
select <- dplyr::select  # namespace conflicts...

# ── Settings ─────────────────────────────────────────────────
# sum-to-zero contrasts
options(contrasts = c("contr.sum", "contr.poly"))

#Sets up dirtecotry tree and globals
infile      <- "./inputs/cyst_csv.csv"
outdir      <- "./outputs/H2Av_output"
stage_order <- c("GSC","CB","2cc","4cc","8cc","16cc")
endpoints   <- c("Continuous","Discontinuous","Foci")
geno_col    <- c("WT" = "#00B7C3", "meiP26" = "#C3007A")
c3g_colors <- c("Foci" = "#E69F00", "Discontinuous" = "#56B4E9", "Continuous" = "#009E73")


if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
int_dir <- file.path(outdir, "int") #intermediate files
if (!dir.exists(int_dir)) dir.create(int_dir, recursive = TRUE)

#ggplot theme
base_theme <- theme_bw(base_size = 12) +
  theme(text = element_text(size=12),
        legend.title     = element_blank(),
        panel.grid.minor = element_blank())

# Star function
p_star <- function(p) dplyr::case_when(
  is.na(p) ~ "", p<0.001 ~ "***", p<0.01 ~ "**", p<0.05 ~ "*", p<0.1 ~ ".", TRUE ~ "")

# Data Wrangling
#base data
dat <- read.csv(infile, stringsAsFactors=FALSE, check.names=FALSE)
names(dat) <- trimws(names(dat))

# ── H2AV Wrangling ─────────────────────────────────────────────
# updated_cyst_csv.csv has Cyst_ID, Genotype, MStage, and New_ID
# casts cat fields as factors
# Removes na data, anything out of possible intensity ranges, and duplicates
# Cleans imageId and calcs log of normalized data
intensityDat <- dat %>%
  filter(!is.na(Cyst_ID), !is.na(Genotype), !is.na(MStage),
         !is.na(Normalized_H2Av), Normalized_H2Av > 0,
         Genotype %in% c("WT","meiP26"), MStage %in% stage_order) %>%
  distinct(Cyst_ID, .keep_all=TRUE) %>% #removes duplicates, preserving first row only
  mutate(
    imageId  = factor(trimws(New_ID)),
    Genotype = factor(Genotype, levels=c("WT","meiP26")),
    MStage   = factor(MStage,   levels=stage_order),
    log_H2Av = log(Normalized_H2Av),
    log_H2Av_masked = log(Normalized_H2Av_masked),
    log_C3G = log(Normalized_C3G)
  )

# ── C3G Wrangling ─────────────────────────────────────────────
phenoDat <- dat %>%
  rename_with(trimws) %>%
  filter(!is.na(Cyst_ID), !is.na(Genotype), !is.na(MStage),
         Genotype %in% c("WT","meiP26"), MStage %in% stage_order) %>%
  distinct(Cyst_ID, .keep_all = TRUE) %>%
  mutate(
    imageId       = factor(trimws(New_ID)),
    Genotype      = factor(Genotype, levels = c("WT","meiP26")),
    MStage        = factor(MStage,   levels = stage_order),
    Continuous    = as.integer(c3G_continuous),
    Discontinuous = as.integer(c3G_discontinuous),
    Foci          = as.integer(c3G_foci)
  )
phenoDat <- phenoDat %>%mutate(c3g_sum = Foci + Discontinuous + Continuous)

# ── Zero-fill per-germarium counts ───────────────────────────
# A germarium with no cysts of a given stage gets 0 — not a
# missing row — so variability stats include all observations.
dat <- read.csv(infile, stringsAsFactors=FALSE, check.names=FALSE)
names(dat) <- trimws(names(dat))

dat3 <- dat %>%
  filter(!is.na(Cyst_ID), !is.na(Genotype), !is.na(MStage),
         Genotype %in% c("WT","meiP26"), MStage %in% stage_order) %>%
  mutate(germariumId = trimws(New_ID)) %>%
  distinct(germariumId, Cyst_ID, .keep_all=TRUE) %>%
  mutate(
    Genotype = factor(Genotype, levels=c("WT","meiP26")),
    MStage   = factor(MStage,   levels=stage_order)
  )

cyst_raw <- dat3 %>%
  group_by(germariumId, Genotype, MStage) %>%
  summarise(Cyst_Count=n(), .groups="drop")

cystPG <- cyst_raw %>%
  complete(nesting(germariumId, Genotype),
           MStage=factor(stage_order, levels=stage_order),
           fill=list(Cyst_Count=0)) %>% #map zeros
  mutate(MStage=factor(MStage, levels=stage_order))




