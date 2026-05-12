# ============================================================
# H2Av intensity analysis
# ============================================================

# Summary
# -------------
# H2Av intensity was analyzed using a log-linked Generalized Linear Mixed Model (GLMM)
# via the glmmTMB package in R. The model included Genotype, Developmental Stage, and
# their interaction as fixed effects, with Image ID as a random intercept to
# account for technical variation between imaging sessions. To address heteroscedasticity
# and quantify changes in developmental precision, residual variance was explicitly
# modeled using a dispersion formula (~Genotype*Stage). Statistical significance was
# determined via Type III Wald Chi-square tests, followed by post-hoc comparisons
# of Estimated Marginal Means (EMMeans) with p-values adjusted for multiple testing.

source("./loadData.R", chdir = TRUE)

# ── Image-level summary & Wilcoxon sanity check ──────────────
# Data are pseudo-replicated at image level.
# This aggregates data to genotype to assess differences at that level
# Data are non-parametric (see boxplot) so use Wilcoxon sign-rank
# This is a high-level assessment of genotype differences, but is only marginally
# informative.
img_dat <- intensityDat %>%
  group_by(imageId, Genotype) %>%
  summarise(image_mean_H2Av = mean(Normalized_H2Av), n_cysts=n(), .groups="drop")

boxplot(img_dat$image_mean_H2Av~img_dat$Genotype)
cat("\nWilcoxon on image means (conservative check only):\n")
print(wilcox.test(image_mean_H2Av ~ Genotype, data=img_dat))

# ── Model Data ────────────────────────────────────────
# Assess via linear model. Most parsimonious
model <- lmer(Normalized_H2Av ~ Genotype * MStage + (1 | imageId),
              data=intensityDat, REML=TRUE)

# Residual plots Normalized Data FAIL
par(mfrow=c(1,2))
plot(fitted(model), residuals(model), xlab="Fitted", ylab="Residuals")
abline(h=0, lty=2)
qqnorm(residuals(model)); qqline(residuals(model))
par(mfrow=c(1,1))

# Log data instead to control variance
model <- lmer(log_H2Av ~ Genotype * MStage + (1 | imageId),
              data=intensityDat)

# Residual plots FAIL (Heteroskedastic)
par(mfrow=c(1,2))
plot(fitted(model), residuals(model), xlab="Fitted", ylab="Residuals")
abline(h=0, lty=2)
qqnorm(residuals(model)); qqline(residuals(model))
par(mfrow=c(1,1))

# Variance in data is not handled by model, resulting heterscedasticity in residuals
# Switch to generalized linear mixed models and model heterscedasticity as a process

# Heterogeneous Variance Generalized Linear Mixed Model (GLMM)
model <- glmmTMB(
  formula = log_H2Av ~ Genotype * MStage + (1 | imageId),
  data = intensityDat,
  # Dispersion Formula: Explicitly models unequal variance across Genotypes and Stages;
  # necessary to resolve significant heteroscedasticity detected.
  dispformula = ~ MStage*Genotype
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

# Test: Residual spread across groups (checking if dispersion formula worked)
# Look for: Similar "heights" and spreads across all genotypes and stages
plotResiduals(res, form = intensityDat$Genotype)
plotResiduals(res, form = intensityDat$MStage)

# Test: Multicollinearity (VIF)
# Look for: Adjusted VIF < 10 (confirms Genotype/Stage are independent)
check_collinearity(model)

# Test: Likelihood-ratio R2 (Cox & Snell, 1989)
# Dispersion model made Nakagawa method unavailable
# Calculated the Likelihood-ratio r2. However, this assesses entire fit
# It does not distinguish between fixed, dispersed, and random effects
model_null <- glmmTMB(log_H2Av ~ 1 + (1 | imageId), data = intensityDat, family = gaussian())
r.squaredLR(model, null = model_null)

# Test: Partitioning and normality of image-level noise
# Look for: Standard Deviation (VarCorr) and Shapiro-Wilk p > 0.05
VarCorr(model) #image noise

re <- ranef(model)$cond$imageId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)

# Model Results
print(summary(model))
# Fixed Effects
print(car::Anova(model, type = 3, test.statistic = "Chisq"))
# Dispersion
print(car::Anova(model, type = 3, component = 'disp', test.statistic = "Chisq"))

# ── model estimated means ─────────────────────────
# Fixed Effects
emm <- emmeans(model, ~ Genotype | MStage, tran = "log", type = "response")
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
  labs(y = "H2Av Intensity (Units)", x = "Cyst Stage") +
  base_theme
p_emm

emm_resp

# ── model estimated residual variance ─────────────────────────
# Residual Variance
emm_disp <- emmeans(model, ~ Genotype | MStage, component = "disp", type = "response")
disp_pairs <- pairs(emm_disp, reverse = TRUE)
plot_disp <- as.data.frame(emm_disp)

# Adds p-value stars to plots
star_df <- as.data.frame(disp_pairs) %>%
  mutate(star = p_star(p.value)) %>%
  filter(star != "") %>%
  left_join(as.data.frame(emm_disp) %>% group_by(MStage) %>%
              summarise(y_pos=max(asymp.UCL)*1.10, .groups="drop"),
            by="MStage")
#Emmeans
ggplot(plot_disp, aes(x = MStage, y = response, color = Genotype, group = Genotype)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df)>0)
    geom_text(data=star_df, aes(x=MStage, y=y_pos, label=star),
              inherit.aes=FALSE, color="black", size=5) } +
  labs(title = "System Noise (Dispersion)",
       y = "Estimated Variance Component",
       x = "Cyst Stage") +
  scale_color_manual(values = geno_col) +
  base_theme

summary(disp_pairs)


# ── Model 16cc Regions ────────────────────────────────────────
##################################################################
intensityDat16cc = filter(intensityDat, MStage == '16cc')

# Generalized Linear Mixed Model (GLMM)
# We assessed the normalized and log normalized h2av
# Both were heteroskedastic (not shown here) so using the dispersion model again
model <- glmmTMB(
  formula = log_H2Av_masked ~ Genotype * Region + (1 | imageId),
  data = intensityDat16cc,
  dispformula = ~ Genotype*Region
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

# Test: Multicollinearity (VIF)
# Look for: Adjusted VIF < 10 (confirms Genotype/Stage are independent)
check_collinearity(model)

# Test: Likelihood-ratio R2 (Cox & Snell, 1989)
# Dispersion model made Nakagawa method unavailable
# Calculated the Likelihood-ratio r2. However, this assesses entire fit
# It does not distinguish between fixed, dispersed, and random effects
model_null <- glmmTMB(log_H2Av_masked ~ 1 + (1 | imageId), data = intensityDat16cc)
r.squaredLR(model, null = model_null)

# Test: Partitioning and normality of image-level noise
# Look for: Standard Deviation (VarCorr) and Shapiro-Wilk p > 0.05
VarCorr(model) #image noise

re <- ranef(model)$cond$imageId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)


# Model Results
print(summary(model))
# Fixed Effects
print(car::Anova(model, type = 3, test.statistic = "Chisq"))
# Dispersion
print(car::Anova(model, type = 3, component = 'disp', test.statistic = "Chisq"))

# ── model estimated means ─────────────────────────
# Fixed Effects
emm <- emmeans(model, ~ Genotype | Region, tran = "log", type = "response")
emm_pairs <- pairs(emm, reverse = TRUE)
emm_resp   <- summary(emm_pairs)
emm_df <- as.data.frame(emm)
# Adds p-value stars to plots
star_df <- as.data.frame(emm_pairs) %>%
  mutate(star = p_star(p.value)) %>%
  filter(star != "") %>%
  left_join(emm_df %>% group_by(Region) %>%
              summarise(y_pos=max(asymp.UCL)*1.10, .groups="drop"),
            by="Region")
#Emmeans
p_emm = ggplot(emm_df, aes(x = Region, y = response, color = Genotype, group = Genotype)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df)>0)
    geom_text(data=star_df, aes(x=Region, y=y_pos, label=star),
              inherit.aes=FALSE, color="black", size=5) } +
  scale_color_manual(values = geno_col) +
  labs(y = "H2Av Intensity (Units)", x = "Cyst Stage") +
  base_theme
p_emm

emm_resp

# ── model estimated residual variance ─────────────────────────
# Residual Variance
emm_disp <- emmeans(model, ~ Genotype | Region, component = "disp", type = "response")
disp_pairs <- pairs(emm_disp, reverse = TRUE)
plot_disp <- as.data.frame(emm_disp)

# Adds p-value stars to plots
star_df <- as.data.frame(disp_pairs) %>%
  mutate(star = p_star(p.value)) %>%
  filter(star != "") %>%
  left_join(as.data.frame(emm_disp) %>% group_by(Region) %>%
              summarise(y_pos=max(asymp.UCL)*1.10, .groups="drop"),
            by="Region")
#Emmeans
ggplot(plot_disp, aes(x = Region, y = response, color = Genotype, group = Genotype)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df)>0)
    geom_text(data=star_df, aes(x=Region, y=y_pos, label=star),
              inherit.aes=FALSE, color="black", size=5) } +
  labs(title = "System Noise (Dispersion)",
       y = "Estimated Variance Component",
       x = "Cyst Stage") +
  scale_color_manual(values = geno_col) +
  base_theme

summary(disp_pairs)


