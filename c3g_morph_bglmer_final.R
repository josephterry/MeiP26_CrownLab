# ============================================================
# C3G morphology — Bayesian-penalized binomial mixed models
# Endpoints: Continuous, Discontinuous, Foci (binary 0/1)
# ============================================================

# Summary
# -----------------
# C3G morphological phenotypes were analyzed using Bayesian-penalized binomial
# mixed models (blme) to provide stable estimates in the presence of complete
# separation (0% or 100% frequencies). The models included Genotype, Stage,
# and their interaction as fixed effects, with Image ID as a random intercept
# and a weakly informative Normal prior (Gelman et al., 2008) applied to fixed-effect
# coefficients. Significance was determined via Type III Wald Chi-square tests
# followed by Tukey-adjusted pairwise comparisons of Estimated Marginal Means.

source("./loadData.R", chdir = TRUE)

# ── FOCI ───────────────────────────────────────────────
p = ncol(model.matrix(~ MStage * Genotype, data = phenoDat)) #number of coeffients in the model
model = bglmer(Foci ~ MStage * Genotype + (1 | imageId),
               data = phenoDat,
               family = binomial, # Used for binary/dichotomous data (0/1)

               # fixef.prior: Handles "Complete Separation" (0% or 100% cells) by applying a weak Bayesian penalty to keep estimates finite.
               # We choose a Normal prior because it is the most mathematically conservative way to apply a symmetric "penalty" that pulls extreme, infinite values back toward a stable, finite center without biasing the rest of the phenoData.
               # Gelman et al. (2008) recommend 2.5 as a starting point as it is weakly informative.
               fixef.prior = normal(cov = diag(2.5, p)),

               # This is a prior on the random effects. It keeps the model from assuming variance between photos is zero, which would ignore pseudoreplication.
               # Bglmer must apply a prior, so we choose gamma since it is better for single fixed effects
               cov.prior = gamma,

               # control: Specifies the "bobyqa" optimizer, which is more robust than the default for complex mixed models and unbalanced designs.
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6))
)

#Model Diags
# Test: Optimizer convergence and model stability
# Look for: No "failed to converge" errors; Singularity = FALSE
check_convergence(model)
check_singularity(model)

# Test: Simulated quantile residuals (DHARMa)
# Look for: Non-significant p-values (p > 0.05) and points on the diagonal
# NOTE: DHARMa does not formally support bglmer objects, but bglmer inherits from the
# lme4 framework used by supported GLMMs. Residual simulations are therefore
# used here as an approximate check for dispersion and gross model misspecification.
res <- simulateResiduals(model, n = 1000)
plotResiduals(res)
testUniformity(res)
testDispersion(res)

# Test: Multicollinearity (VIF) and Bayesian Prior impact
# Look for: Adj. VIF < 10
check_collinearity(model)

# Test: Discrimination (AUC) and Explanatory power (R2)
# Look for: AUC > 0.8; High Marginal R2 (Signal vs. Noise)
prob <- predict(model, type = "response")
roc_obj <- roc(phenoDat$Foci, prob)
plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 2)))

r2_nakagawa(model)

# Test: Inter-image noise distribution and impact (ICC)
# Look for: ICC < 0.2 (indicates image noise doesn't drown out biology)
VarCorr(model)
# Test: Normality and impact of noise; Significant p-value is okay if ICC is low
re <- ranef(model)$imageId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)
icc(model)


#Model Results
print(summary(model))
print(car::Anova(model, type = 3, test.statistic = "Chisq"))

#Emmeans
stageOrder <- c('GSC', 'CB', '2cc', '4cc', '8cc', '16cc')
emm_out <- emmeans(model, ~ Genotype | MStage, type = "response")
emmDf   <- as.data.frame(emm_out) # Contains 'prob', 'asymp.LCL', 'asymp.UCL'
emmDf$MStage <- factor(emmDf$MStage, levels = stageOrder, ordered = TRUE)
geno_pairs <- pairs(emm_out, reverse = TRUE)

star_df <- as.data.frame(geno_pairs) %>%
  mutate(star = p_star(p.value)) %>%
  filter(star != "")
star_df$MStage <- factor(star_df$MStage, levels = stageOrder, ordered = TRUE)
star_df <- star_df %>%
  left_join(emmDf %>% group_by(MStage) %>%
              summarise(y_pos = max(asymp.UCL) + 0.05, .groups="drop"),
            by = "MStage")

ggplot(emmDf, aes(x = MStage, y = prob, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df) > 0)
    geom_text(data = star_df, aes(x = MStage, y = y_pos, label = star),
              inherit.aes = FALSE, color = "black", size = 6, vjust = 0) } +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  scale_color_manual(values = geno_col) +
  labs(y = "Probability of Foci C3G (Odds Ratio)", x = "Stage") +
  base_theme

summary(geno_pairs)

# ── Discontinuous ───────────────────────────────────────────────
# Analysis was restricted to the maturation window (8cc–16cc),
# where synaptonemal complex elongation is biologically relevant.
# This focused approach also prevents statistical instability (complete separation)
# associated with the structural zeros inherent to early germline development.
phenoDat_discontinuous <- phenoDat %>%
  filter(MStage %in% c("8cc", "16cc")) %>%
  mutate(MStage = factor(MStage),     # This drops GSC, CB, etc.
         Genotype = factor(Genotype))
p = ncol(model.matrix(~ MStage * Genotype, data = phenoDat_discontinuous))
model = bglmer(Discontinuous ~ MStage * Genotype + (1 | imageId),
               data = phenoDat_discontinuous,
               family = binomial,
               fixef.prior = normal(cov = diag(2.5, p)),
               cov.prior = gamma
)

#Model Diags
# Test: Optimizer convergence and model stability
# Look for: No "failed to converge" errors; Singularity = FALSE
check_convergence(model)
check_singularity(model)

# Test: Simulated quantile residuals (DHARMa)
# Look for: Non-significant p-values (p > 0.05) and points on the diagonal
res <- simulateResiduals(model, n = 1000)
plotResiduals(res)
testUniformity(res)
testDispersion(res)

# Test: Multicollinearity (VIF) and Bayesian Prior impact
# Look for: Adj. VIF < 10
check_collinearity(model)

# Test: Discrimination (AUC) and Explanatory power (R2)
# Look for: AUC > 0.8; High Marginal R2 (Signal vs. Noise)
prob <- predict(model, type = "response")
roc_obj <- roc(phenoDat_discontinuous$Discontinuous, prob)
plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 2)))

r2_nakagawa(model)

# Test: Inter-image noise distribution and impact (ICC)
# Look for: ICC < 0.2 (indicates image noise doesn't drown out biology)
VarCorr(model)
# Test: Normality and impact of noise; Significant p-value is okay if ICC is low
re <- ranef(model)$imageId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)
icc(model)

#Model Results
print(summary(model))
print(car::Anova(model, type = 3))

#Emmeans
stageOrder <- c('GSC', 'CB', '2cc', '4cc', '8cc', '16cc')
emm_out <- emmeans(model, ~ Genotype | MStage, type = "response")
emmDf   <- as.data.frame(emm_out) # Contains 'prob', 'asymp.LCL', 'asymp.UCL'
emmDf$MStage <- factor(emmDf$MStage, levels = stageOrder, ordered = TRUE)
geno_pairs <- pairs(emm_out, reverse = TRUE)

star_df <- as.data.frame(geno_pairs) %>%
  mutate(star = p_star(p.value)) %>%
  filter(star != "")
star_df$MStage <- factor(star_df$MStage, levels = stageOrder, ordered = TRUE)
star_df <- star_df %>%
  left_join(emmDf %>% group_by(MStage) %>%
              summarise(y_pos = max(asymp.UCL) + 0.05, .groups="drop"),
            by = "MStage")

ggplot(emmDf, aes(x = MStage, y = prob, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df) > 0)
    geom_text(data = star_df, aes(x = MStage, y = y_pos, label = star),
              inherit.aes = FALSE, color = "black", size = 6, vjust = 0) } +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  scale_color_manual(values = geno_col) +
  labs(y = "Probability of Discontinuous C3G", x = "Stage") +
  base_theme

summary(geno_pairs)

# ── Continuous ───────────────────────────────────────────────
# Analysis was restricted to the maturation window (8cc–16cc),
# where synaptonemal complex elongation is biologically relevant.
# This focused approach also prevents statistical instability (complete separation)
# associated with the structural zeros inherent to early germline development.
phenoDat_continuous <- phenoDat %>%
  filter(MStage %in% c("8cc", "16cc")) %>%
  mutate(MStage = factor(MStage),     # This drops GSC, CB, etc.
         Genotype = factor(Genotype))
p = ncol(model.matrix(~ MStage * Genotype, data = phenoDat_continuous)) #number of coeffients in the model; we need to know this for the fixef priors
model = bglmer(Continuous ~ MStage * Genotype + (1 | imageId),
               data = phenoDat_continuous,
               family = binomial,
               fixef.prior = normal(cov = diag(2.5, p)), #use 2.5, 1, or 0.1... each decrease, increases chance of Type II error
               cov.prior = gamma,
)

#Model Diags
# Test: Optimizer convergence and model stability
# Look for: No "failed to converge" errors; Singularity = FALSE
check_convergence(model)
check_singularity(model)

# Test: Simulated quantile residuals (DHARMa)
# Look for: Non-significant p-values (p > 0.05) and points on the diagonal
res <- simulateResiduals(model, n = 1000)
plotResiduals(res)
testUniformity(res)
testDispersion(res)

# Test: Multicollinearity (VIF) and Bayesian Prior impact
# Look for: Adj. VIF < 10
check_collinearity(model)

# Test: Discrimination (AUC) and Explanatory power (R2)
# Look for: AUC > 0.8; High Marginal R2 (Signal vs. Noise)
prob <- predict(model, type = "response")
roc_obj <- roc(phenoDat_continuous$Continuous, prob)
plot(roc_obj, main = paste("AUC =", round(auc(roc_obj), 2)))

r2_nakagawa(model)

# Test: Inter-image noise distribution and impact (ICC)
# Look for: ICC < 0.2 (indicates image noise doesn't drown out biology)
VarCorr(model)
# Test: Normality and impact of noise; Significant p-value is okay if ICC is low
re <- ranef(model)$imageId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)
icc(model)


#Model Results
print(summary(model))
print(car::Anova(model, type = 3))


#Emmeans. Interaction is non-significant. Plotting emmeans but not showing significance
# There is no genotype contrast within each stage.
# There are stage contrasts and there are genotype contrasts, but not genotype by stage constrasts.
stageOrder <- c('GSC', 'CB', '2cc', '4cc', '8cc', '16cc')
emm_out <- emmeans(model, ~ Genotype | MStage, type = "response")
emmDf   <- as.data.frame(emm_out) # Contains 'prob', 'asymp.LCL', 'asymp.UCL'
emmDf$MStage <- factor(emmDf$MStage, levels = stageOrder, ordered = TRUE)


ggplot(emmDf, aes(x = MStage, y = prob, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  scale_color_manual(values = geno_col) +
  labs(y = "Probability of Continuous C3G", x = "Stage") +
  base_theme

# Note that while emmeans are significant, the interaction is not.
# These should be interpreted with caution. They are likely from the genotype main effect.
summary(pairs(emm_out, reverse = TRUE))

