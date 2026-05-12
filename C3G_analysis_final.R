# ============================================================
# C3G intensity analysis
# ============================================================

# Summary
###########################
# C3G protein intensity was analyzed using a zero-inflated Gamma generalized linear mixed model
# (GLMM), incorporating Image ID as a random intercept to account for hierarchical nesting.
# To address the bimodal distribution of the data, signal initiation was independently modeled
# through a zero-inflation component (treating cysts ≤ 80 units as non-expressing; Wildtype floor),
# while the intensity of the expressing population was analyzed via a Gamma distribution with a
# log-link to accommodate right-skewness. Significance was determined using Type III Wald
# Chi-square tests for the intensity component. Post-hoc inference was performed on estimated
# marginal means (EMMs) with Tukey-adjusted pairwise comparisons, and model integrity was
# confirmed using DHARMa residual simulations.

source("./loadData.R", chdir = TRUE)


# ── Image-level summary & Wilcoxon sanity check ──────────────
# Data are pseudo-replicated at image level.
# This aggregates data to genotype to assess differences at that level
# Data are non-parametric (see boxplot) so use Wilcoxon sign-rank
# This is a high-level assessment of genotype differences, but is only marginally
# informative.
img_dat <- intensityDat %>%
  group_by(imageId, Genotype) %>%
  summarise(image_mean_C3G = mean(Normalized_C3G), n_cysts=n(), .groups="drop")

boxplot(img_dat$image_mean_C3G~img_dat$Genotype)
cat("\nWilcoxon on image means (conservative check only):\n")
print(wilcox.test(image_mean_C3G ~ Genotype, data=img_dat))

# ── Model Data ────────────────────────────────────────
# Model selection rationale:
# Initial Gaussian, log-Gaussian, Gamma, and Tweedie mixed models all showed
# persistent residual violations under DHARMa diagnostics (not shown), even when modeling
# heterogeneous dispersion. Inspection of genotype-stratified distributions
# revealed a distinct low-intensity subpopulation unique to meiP26 cysts,
# suggesting a biologically distinct low-recruitment state rather than simple
# variance inflation. Because WT cysts showed a stable empirical lower bound
# near 80 intensity units and scatterplots revealed a clear discontinuity below
# this range, values <=80 were classified as a low-recruitment state and mapped
# to zero. A zero-inflated Gamma GLMM substantially improved model diagnostics,
# supporting a two-process interpretation of the phenotype: recruitment failure
# probability plus conditional C3G intensity among recruited cysts.
intensityDat <- intensityDat %>%
  mutate(z = ifelse(Normalized_C3G > 80, Normalized_C3G, 0))

model <- glmmTMB(
  z ~ Genotype * MStage + (1 | imageId),
  # Do not use interaction here. Previous test (not shown) had strong collinearity in the interaction
  ziformula = ~ Genotype + MStage, # Interaction terms were excluded from the zero-inflation component because low-intensity states were confined to the meiP26 genotype
  family = ziGamma(link = "log"),
  data = intensityDat
)

#Model Diags
# Test: Optimizer convergence and model stability
# Look for: No "failed to converge" errors; Singularity = FALSE
check_convergence(model)
check_singularity(model)

# Test: Simulated quantile residuals (DHARMa)
# Look for: p > 0.05; no significant outliers; points follow diagonal line
res = simulateResiduals(model, n=1000)
plotResiduals(res)
testUniformity(res)
testDispersion(res)
testOutliers(res, type = "bootstrap")
testZeroInflation(res)

#Homogeneity of variance within groups
plotResiduals(res, form = intensityDat$Genotype)
plotResiduals(res, form = intensityDat$MStage)

# Test: Multicollinearity (VIF)
# Look for: Adjusted VIF < 10 (confirms Genotype/Stage are independent)
check_collinearity(model)

# Explanatory Power
r2_nakagawa(model)

# Test: Partitioning and normality of image-level noise
# Look for: Standard Deviation (VarCorr) and Shapiro-Wilk p > 0.05
VarCorr(model) #image noise

re <- ranef(model)$cond$imageId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)
icc(model)

# Model Results
print(summary(model))
# Fixed Effects
print(car::Anova(model, type = 3, test.statistic = "Chisq"))
# Zero-inflation
# The wald test fails here at genotype p ~ 1
# However, AIC favors genotype in model (see below)
# Treating zero-inflation as a correction for low-expressing cysts
print(car::Anova(model, type = 3, component = 'zi', test.statistic = "Chisq"))
model_no_mstage <- update(model, ziformula = ~ Genotype)
model_no_geno <- update(model, ziformula = ~ MStage)
model_int <- update(model, ziformula = ~ Genotype*MStage)
anova(model, model_no_mstage, model_no_geno, model_int)

# ── model estimated means ─────────────────────────
# Fixed Effects
emm <- emmeans(model, ~ Genotype | MStage, type='response')
emm_pairs <- pairs(emm, reverse = TRUE)
emm_resp   <- summary(emm_pairs)
emm_df <- as.data.frame(emm)
# Adds p-value stars to plots
star_df <- as.data.frame(emm_pairs) %>%
  mutate(star = p_star(p.value)) %>%
  filter(star != "") %>%
  left_join(emm_df %>% group_by(MStage) %>%
              summarise(y_pos=max(asymp.UCL)*1.10, .groups="drop"),
            by="MStage")
#Emmeans
p_emm = ggplot(emm_df, aes(x = MStage, y = response, color = Genotype, group = Genotype)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df)>0)
    geom_text(data=star_df, aes(x=MStage, y=y_pos, label=star),
              inherit.aes=FALSE, color="black", size=5) } +
  scale_color_manual(values = geno_col) +
  labs(y = "C3G Intensity (Units)", x = "Cyst Stage") +
  base_theme
p_emm

emm_resp
