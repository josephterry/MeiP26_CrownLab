# ============================================================
# Cyst composition analysis
# Count models
# ============================================================

source("./loadData.R", chdir = TRUE)

# ── Model the count data ───────────────────────────
# We used a single-process conditional GLMM because cyst counts represent a continuous developmental trajectory
# rather than distinct discrete states requiring a hurdle or zero-inflation process.
#
# Initial Poisson and negative binomial GLMMs showed systematic lack of fit in DHARMa diagnostics,
# indicating mismatch between assumed and observed mean–variance relationships.
#
# We therefore compared flexible count distributions, including generalized Poisson and Conway-Maxwell Poisson,
# both of which allow relaxation of the Poisson mean–variance constraint.
#
# The COM-Poisson model provided the best overall fit based on AIC and DHARMa residual diagnostics,
# and was selected as the final model.

#Poisson
model = glmmTMB(
  Cyst_Count ~ Genotype * MStage + (1 | germariumId),
  data = cystPG,
  family = poisson(),
)

# Residuals test FAIL
res = simulateResiduals(model)
plot(res)

#Approxiimates quasi-poisson by adding obs as a random effect
cystPG$obs <- 1:nrow(cystPG)
model = glmmTMB(
  Cyst_Count ~ Genotype * MStage + (1 | germariumId) + (1 | obs),
  data = cystPG,
  family = poisson(),
)

# Residuals test FAIL
res = simulateResiduals(model)
plot(res)


# Negative binomial; note nb1 fails to converge
model = glmmTMB(
  Cyst_Count ~ Genotype * MStage + (1 | germariumId),
  data = cystPG,
  family = nbinom2(),
)

# Residuals test FAIL
res = simulateResiduals(model)
plot(res)

#Switch to a model with more flexible mean-variance relationships
#Assessing generalized vs conway-maxwell poisson
genPoisModel = glmmTMB(
  Cyst_Count ~ Genotype * MStage + (1 | germariumId),
  data = cystPG,
  family = genpois(),
)

# Residuals test PASS with some outlier effects
res = simulateResiduals(genPoisModel)
plot(res)
plotResiduals(res)
testUniformity(res)
testDispersion(res)
testOutliers(res, type = "bootstrap")

cwPoisModel = glmmTMB(
  Cyst_Count ~ Genotype * MStage + (1 | germariumId),
  data = cystPG,
  family = compois(),
)

# Residuals test PASS
res = simulateResiduals(cwPoisModel)
plot(res)
plotResiduals(res)
testUniformity(res)
testDispersion(res)
testOutliers(res, type = "bootstrap")

# Residual Checks pass for CW and GenPois (marginally)
# Compare via AIC and LLH
anova(genPoisModel, cwPoisModel)

# Model comparison favors the Conway Maxwell Poisson distribution
model = cwPoisModel

#Model Diags
# Test: Optimizer convergence and model stability
# Look for: No "failed to converge" errors; Singularity = FALSE
check_convergence(model)
check_singularity(model) #It is singular... meaning the RE effect is nearly 0

# Test: Simulated quantile residuals (DHARMa)
# Look for: p > 0.05; no significant outliers; points follow diagonal line
res = simulateResiduals(model, n=1000)
plotResiduals(res)
testUniformity(res)
testDispersion(res)
testOutliers(res, type = "bootstrap")
testZeroInflation(res)

# Levene Tests
plotResiduals(res, cystPG$Genotype)
plotResiduals(res, cystPG$MStage)

# Test: Multicollinearity (VIF)
# Look for: Adjusted VIF < 10 (confirms Genotype/Stage are independent)
check_collinearity(model)

# Test: Partitioning and normality of image-level noise
# Note these are not needed with RE variance collapsed to 0. Keeping for completeness
VarCorr(model) #image noise
re <- ranef(model)$cond$germariumId
qqnorm(re$`(Intercept)`); qqline(re$`(Intercept)`)
shapiro.test(re$`(Intercept)`)

# Explanatory Power
r2_nakagawa(model) #Fails.... there random effect variance is 0
# Calculated the Likelihood-ratio r2. However, this assesses entire fit
# It does not distinguish between fixed, dispersed, and random effects
model_null <- glmmTMB(Cyst_Count ~ 1, data = cystPG, family = compois())
r.squaredLR(model, null = model_null)

#Compare compois without the random effect
model_no_re <- glmmTMB(
  Cyst_Count ~ Genotype * MStage,
  data = cystPG,
  family = compois()
)
anova(model, model_no_re)
# A random intercept for germarium identity was included to reflect the experimental design.
# However, its variance was estimated to be near zero, suggesting minimal additional clustering beyond the fixed effects;
# results were robust to models fitted with and without this term (AIC diff is 2; P value is 1).
# Maintaining the random effect to reflect the experimental design

#Model Results
print(summary(model))
print(car::Anova(model, type = 3))

#Emmeans.
emm <- emmeans(model, ~ Genotype | MStage, type = "response")
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

p_emm = ggplot(emm_df, aes(x = MStage, y = response, color = Genotype, group = Genotype)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  { if (nrow(star_df)>0)
    geom_text(data=star_df, aes(x=MStage, y=y_pos, label=star),
              inherit.aes=FALSE, color="black", size=5) } +
  scale_color_manual(values = geno_col) +
  labs(y = "Cyst Count", x = "Cyst Stage") +
  base_theme
p_emm

emm_resp
