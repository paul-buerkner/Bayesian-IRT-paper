# load required pacakges
library(tidyverse)
library(brms)
# for comparison with brms
library(lme4)
library(TAM)

# set ggplot theme
theme_set(bayesplot::theme_default())

# set rstan options
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))


# ----------- Code for Section 5.1 ------------
# Analysis of the VerbAgg data set using dichotomous IRT models
data("VerbAgg", package = "lme4")

# get an overview of the data
head(VerbAgg, 10)

# ---------- 1PL models ----------------------
# specify a 1PL model in brms
formula_va_1pl <- bf(r2 ~ 1 + (1 | item) + (1 | id))

# specify some weakly informative priors
prior_va_1pl <- 
  prior("normal(0, 3)", class = "sd", group = "id") + 
  prior("normal(0, 3)", class = "sd", group = "item")

# fit the 1PL model
fit_va_1pl <- brm(
  formula = formula_va_1pl,
  data = VerbAgg,
  family = brmsfamily("bernoulli", "logit"),
  prior = prior_va_1pl
)

# obtain basic summaries
summary(fit_va_1pl)
plot(fit_va_1pl, ask = FALSE)

# extract item parameters
(item_pars_va_1pl <- coef(fit_va_1pl)$item)

# plot person parameters
item_pars_va_1pl[, , "Intercept"] %>%
	as_tibble() %>%
	rownames_to_column() %>%
	rename(item = "rowname") %>%
	mutate(item = as.numeric(item)) %>%
	ggplot(aes(item, Estimate, ymin = Q2.5, ymax = Q97.5)) +
	geom_pointrange() +
	coord_flip() +
	labs(x = "Item Number")

# extract person parameters
(person_pars_va_1pl <- ranef(fit_va_1pl)$id)

# plot person parameters
person_pars_va_1pl[, , "Intercept"] %>%
	as_tibble() %>%
	rownames_to_column() %>%
	arrange(Estimate) %>%
	mutate(id = seq_len(n())) %>%
	ggplot(aes(id, Estimate, ymin = Q2.5, ymax = Q97.5)) +
	geom_pointrange(alpha = 0.7) +
	coord_flip() +
	labs(x = "Person Number (Sorted)")

# specify a 1PL model with lme4 for comparison
lme4_va_1pl <- glmer(
	r2 ~ 1 + (1 | item) + (1 | id),
	data = VerbAgg,
	family = binomial()
)
summary(lme4_va_1pl)

# person and item parameters are similar to those obtained by brms
coef(lme4_va_1pl)$item
ranef(lme4_va_1pl)$id

# specify a 1PL model with TAM for comparison
# bring the data in wide structure first
VerbAgg_wide <- VerbAgg %>%
	select(item, id, r2) %>%
	mutate(r2 = ifelse(r2 == "Y", 1, 0)) %>%
	spread(key = "item", value = "r2") %>%
	select(-id)

# fit the model with TAM 
tam_va_1pl <- tam(VerbAgg_wide, irtmodel = "1PL", verbose = FALSE)

# person and item parameters are similar to those obtained by brms
summary(tam_va_1pl)
IRT.factor.scores(tam_va_1pl)


# ---------- 2PL models ----------------------
## specify a 2PL model
formula_va_2pl <- bf(
  r2 ~ exp(logalpha) * eta,
  eta ~ 1 + (1 |i| item) + (1 | id),
  logalpha ~ 1 + (1 |i| item),
  nl = TRUE
)

# specify some weakly informative priors
prior_va_2pl <- 
  prior("normal(0, 5)", class = "b", nlpar = "eta") +
  prior("normal(0, 1)", class = "b", nlpar = "logalpha") +
  prior("normal(0, 3)", class = "sd", group = "id", nlpar = "eta") + 
  prior("normal(0, 3)", class = "sd", group = "item", nlpar = "eta") +
  prior("normal(0, 1)", class = "sd", group = "item", nlpar = "logalpha")

# fit the 2PL model
fit_va_2pl <- brm(
  formula = formula_va_2pl,
  data = VerbAgg,
  family = brmsfamily("bernoulli", "logit"),
  prior = prior_va_2pl,
)

# obtain some basic summaries
summary(fit_va_2pl)
plot(fit_va_2pl, ask = FALSE)

# extract item parameters
(item_pars_va_1pl <- coef(fit_va_2pl)$item)

# plot item parameters
# difficulties
eta <- item_pars_va_1pl[, , "eta_Intercept"] %>%
	as_tibble() %>%
	rownames_to_column()

# discriminations
alpha <- item_pars_va_1pl[, , "logalpha_Intercept"] %>%
	exp() %>%
	as_tibble() %>%
	rownames_to_column()

# plot difficulties and discrimination next to each other
bind_rows(eta, alpha, .id = "nlpar") %>%
	rename(item = "rowname") %>%
	mutate(item = as.numeric(item)) %>%
	mutate(nlpar = factor(nlpar, labels = c("Easiness", "Discrimination"))) %>%
	ggplot(aes(item, Estimate, ymin = Q2.5, ymax = Q97.5)) +
	facet_wrap("nlpar", scales = "free_x") +
	geom_pointrange() +
	coord_flip() +
	labs(x = "Item Number")


# extract person parameters
(ranef_va_2pl <- ranef(fit_va_2pl)$id)

# plot person parameters
ranef_va_2pl[, , "eta_Intercept"] %>%
	as_tibble() %>%
	rownames_to_column() %>%
	select(-Est.Error) %>%
	arrange(Estimate) %>%
	mutate(id = seq_len(n())) %>%
	ggplot(aes(id, Estimate, ymin = Q2.5, ymax = Q97.5)) +
	geom_pointrange(alpha = 0.7) +
	coord_flip() +
	labs(x = "Person Number (Sorted)")


# specify a model with constant but estimated discrimination
formula_va_1pl_alpha <- bf(
	r2 ~ exp(logalpha) * eta,
	eta ~ 1 + (1 | item) + (1 | id),
	logalpha ~ 1,
	nl = TRUE
)

# specify some weakly informative priors
prior_va_1pl_alpha <- 
	prior("normal(0, 5)", class = "b", nlpar = "eta") +
	prior("normal(0, 1)", class = "b", nlpar = "logalpha") +
	prior("normal(0, 3)", class = "sd", group = "id", nlpar = "eta") + 
	prior("normal(0, 3)", class = "sd", group = "item", nlpar = "eta")

# fit the model
fit_va_1pl_alpha <- brm(
	formula = formula_va_1pl_alpha,
  data = VerbAgg, 
	family = brmsfamily("bernoulli", "logit"),
	prior = prior_va_1pl_alpha,
	inits = 0
)

# perform model comparison via approximate LOO-CV
loo_va_1pl <- loo(fit_va_1pl)
loo_va_1pl_alpha <- loo(fit_va_1pl_alpha)
loo_va_2pl <- loo(fit_va_2pl)
loo_va_compare <- loo_compare(loo_va_1pl, loo_va_1pl_alpha, loo_va_2pl)
print(loo_va_compare, simplify = FALSE)


# ---------- 1PL models with covariates ----------------------
# specify a model including item covariates
formula_va_1pl_cov1 <- bf(
  r2 ~ btype + situ + mode + (1 | item) + (0 + mode | id)
)

fit_va_1pl_cov1 <- brm(
  formula = formula_va_1pl_cov1,
  data = VerbAgg,
  family = brmsfamily("bernoulli", "logit"),
  prior = prior_va_1pl
)

summary(fit_va_1pl_cov1)
marginal_effects(fit_va_1pl_cov1, "mode")

# compare standard deviations
hyp <- "modedo - modewant > 0"
hypothesis(fit_va_1pl_cov1, hyp, class = "sd", group = "id")

# fit a more complex covariate model
formula_va_1pl_cov2 <- bf(
	r2 ~ Anger + Gender + btype + situ + mode + mode:Gender +
	(0 + Gender | item) + (0 + mode | id)
)
fit_va_1pl_cov2 <- brm(
  formula = formula_va_1pl_cov2,
  data = VerbAgg, 
  family = brmsfamily("bernoulli", "logit"),
  prior = prior_va_1pl
)

summary(fit_va_1pl_cov2)
plot(marginal_effects(fit_va_1pl_cov2, c("Anger", "mode:Gender")), ask = FALSE)


# compare convergence of lme4 and brms for a complex covariate model
# does not converge well
glmer_va_1pl_cov_full <- lme4::glmer(
  r2 ~ 1 + Anger + Gender + btype + situ + mode +
	  (1 + Anger + Gender | item) + (1 + btype + situ + mode  | id),
	data = VerbAgg, family = binomial("logit")
)
summary(glmer_va_1pl_cov_full)

# converges nicely and shows sensible results
fit_va_1pl_cov_full <- brm(
  r2 ~ 1 + Anger + Gender + btype + situ + mode +
	  (1 + Anger + Gender | item) + (1 + btype + situ + mode  | id),
  data = VerbAgg,
  family = brmsfamily("bernoulli", "logit"),
  prior = prior_va_1pl
)
summary(fit_va_1pl_cov_full)


# ---------- 3PL models ----------------------

# 3PL model with known guessing parameter
formula_va_3pl <- bf(
  r2 ~ 0.25 + 0.75 * inv_logit(exp(logalpha) * eta),
  eta ~ 1 + (1 |i| item) + (1 | id),
  logalpha ~ 1 + (1 |i| item),
  nl = TRUE
)
family_va_3pl <- brmsfamily("bernoulli", link = "identity")

# 3PL model with unknown guessing parameters
formula_va_3pl <- bf(
  r2 ~ gamma + (1 - gamma) * inv_logit(exp(logalpha) * eta),
  eta ~ 1 + (1 |i| item) + (1 | id),
  logalpha ~ 1 + (1 |i| item),
  logitgamma ~ 1 + (1 |i| item),
  nlf(gamma ~ inv_logit(logitgamma)),
  nl = TRUE
)



# ----------- Code for Section 5.2 ------------
# fit ordinal models to the VerbAgg data

# specify a basic GRM
formula_va_ord_1pl <- bf(resp ~ 1 + (1 | item) + (1 | id))
fit_va_ord_1pl <- brm(
  formula = formula_va_ord_1pl,
  data = VerbAgg,
  family = brmsfamily("cumulative", "logit"),
  prior = prior_va_1pl
)

summary(fit_va_ord_1pl)
plot(fit_va_ord_1pl, ask = FALSE)

# extract item and person parameters
(ranef_va_ord_1pl <- ranef(fit_va_ord_1pl))

# plot person parameters
ranef_va_ord_1pl$id[, , "Intercept"] %>%
	as_tibble() %>%
	rownames_to_column() %>%
	arrange(Estimate) %>%
	mutate(id = seq_len(n())) %>%
	ggplot(aes(id, Estimate, ymin = Q2.5, ymax = Q97.5)) +
	geom_pointrange(alpha = 0.7) +
	coord_flip() +
	labs(x = "Person Number (Sorted)")

# -------------- ordinal 2PL models ---------------
# specify a GRM with varying discriminations
formula_va_ord_2pl <- bf(
  resp ~ 1 + (1 |i| item) + (1 | id),
  disc ~ 1 + (1 |i| item)	
)

# some weakly informative priors
prior_va_ord_2pl <- 
  prior("normal(0, 3)", class = "sd", group = "id") + 
  prior("normal(0, 3)", class = "sd", group = "item") +
  prior("normal(0, 1)", class = "sd", group = "item", dpar = "disc")

# fit the model
fit_va_ord_2pl <- brm(
  formula = formula_va_ord_2pl,
  data = VerbAgg,
  family = brmsfamily("cumulative", "logit"),
  prior = prior_va_ord_2pl
)
summary(fit_va_ord_2pl)

# extract item and person parameters
(ranef_va_ord_2pl <- ranef(fit_va_ord_2pl))

# plot person parameters
# item easinesses (deviations from thresholds)
eta <- ranef_va_ord_2pl$item[, , "Intercept"] %>%
	as_tibble() %>%
	rownames_to_column()

# discriminations
alpha <- ranef_va_ord_2pl$item[, , "disc_Intercept"] %>%
	exp() %>%
	as_tibble() %>%
	rownames_to_column()

# put easinesses and discriminations together
bind_rows(eta, alpha, .id = "nlpar") %>%
	rename(item = "rowname") %>%
	mutate(item = as.numeric(item)) %>%
	mutate(nlpar = factor(nlpar, labels = c("Easiness", "Discrimination"))) %>%
	ggplot(aes(item, Estimate, ymin = Q2.5, ymax = Q97.5)) +
	facet_wrap("nlpar", scales = "free_x") +
	geom_pointrange() +
	coord_flip() +
	labs(x = "Item Number")

# compute correlations between person parameters across models
cbind(
	va_1pl = ranef_va_1pl$id[, "Estimate", "Intercept"],
	va_2pl = ranef_va_2pl$id[, "Estimate", "eta_Intercept"],
	va_ord_1pl = ranef_va_ord_1pl$id[, "Estimate", "Intercept"],
	va_ord_2pl = ranef_va_ord_2pl$id[, "Estimate", "Intercept"]
) %>%
	cor() %>%
	round(3)


# ------- ordinal models with covariates -------------------
# fit a GRM with person and item covariates
formula_va_ord_cov1 <- bf(
	resp ~ Anger + Gender + btype + situ + mode + mode:Gender +
	(0 + Gender | item) + (0 + mode | id)
)
fit_va_ord_cov1 <- brm(
  formula = formula_va_ord_cov1,
  data = VerbAgg, 
  family = brmsfamily("cumulative", "logit"),
  prior = prior_va_1pl
)
summary(fit_va_ord_cov1)

# plot effects of Anger
marginal_effects(fit_va_ord_cov1, effects = "Anger", categorical = TRUE)


# fit a PCM with covariates and a category specific effect of 'Anger'
formula_va_ord_cov2 <- bf(
  resp ~ cs(Anger) + Gender + btype + situ + mode + mode:Gender +
	(0 + Gender | item) + (0 + mode | id)
)

# fit the model
fit_va_ord_cov2 <- brm(
  formula = formula_va_ord_cov2,
  data = VerbAgg, 
  family = brmsfamily("acat", "logit"),
  prior = prior_va_1pl
)

# summarize the results
summary(fit_va_ord_cov2)
marginal_effects(fit_va_ord_cov2, effects = "Anger", categorical = TRUE)




# ----------- Code for Section 5.2 ------------
# Analysis of the rotation data set using response times IRT models
data("rotation", package = "diffIRT")
rotation <- rotation %>%
	as_tibble() %>%
	mutate(person = seq_len(n())) %>%
	gather("key", "value", -person) %>%
	extract("key", into = c("type", "item"), regex = "(.)\\[(.+)\\]") %>%
	spread("type", "value") %>%
	rename(time = T, resp = X) %>%
	mutate(
		rotate = factor(case_when(
		  item %in% c(2, 5, 8) ~ 50,
		  item %in% c(3, 6, 10) ~ 100,
		  item %in% c(1, 4, 7, 9) ~ 150
	  )),
		item = as.numeric(item)
	)

# get an overview of the data
head(rotation, 10)


# specify a distributional exgaussian model
bform_exg1 <- bf(
  time ~ rotate + (1 |p| person) + (1 |i| item),
  sigma ~ rotate + (1 |p| person) + (1 |i| item),
  beta ~ rotate + (1 |p| person) + (1 |i| item)
)

# fit the model
fit_exg1 <- brm(
  bform_exg1, data = rotation,
  family = brmsfamily("exgaussian", link_sigma = "log", link_beta = "log"),
  chains = 4, cores = 4, inits = 0,
  control = list(adapt_delta = 0.99)
)

# summarize the results
summary(fit_exg1)
pp_check(fit_exg1)

# visualize effects of 'rotate'
marginal_effects(fit_exg1, "rotate", dpar = "mu")
marginal_effects(fit_exg1, "rotate", dpar = "sigma")
marginal_effects(fit_exg1, "rotate", dpar = "beta")



# specify a 3-parameter drift diffusion model
bform_drift1 <- bf(
  time | dec(resp) ~ rotate + (1 |p| person) + (1 |i| item),
  bs ~ rotate + (1 |p| person) + (1 |i| item),
  ndt ~ rotate + (1 |p| person) + (1 |i| item),
  bias = 0.5
)

# specify initial values to help the model start sampling
chains <- 4
inits_drift <- list(temp_ndt_Intercept = -3)
inits_drift <- replicate(chains, inits_drift, simplify = FALSE)

# fit the model
fit_drift1 <- brm(
  bform_drift, data = rotation,
  family = brmsfamily("wiener", "log", link_bs = "log", link_ndt = "log"),
  chains = chains, cores = chains,
  inits = inits_drift, init_r = 0.05,
  control = list(adapt_delta = 0.99)
)

# summarize the model
summary(fit_drift1)

# extract item specific parameters
coef(fit_drift1)$item

# plot the effect of 'rotate'
marginal_effects(fit_drift1, "rotate", dpar = "mu")
marginal_effects(fit_drift1, "rotate", dpar = "bs")
marginal_effects(fit_drift1, "rotate", dpar = "ndt")


# specify a drift diffusion model without 
# 'rotate' affecting the boundary separation
bform_drift2 <- bf(
	time | dec(resp) ~ rotate + (1 |p| person) + (1 |i| item),
	bs ~ 1 + (1 |p| person),
	ndt ~ rotate + (1 |p| person) + (1 |i| item),
	bias = 0.5
)

# fit the model
fit_drift2 <- brm(
	bform_drift2, rotation,
	family = wiener("log", link_bs = "log", link_ndt = "log"),
	chains = chains, cores = chains,
	inits = inits_diff, init_r = 0.05,
	control = list(adapt_delta = 0.99)
)

# perform model comparison via approximate LOO-CV
loo_drift1 <- loo(fit_drift1)
loo_drift2 <- loo(fit_drift2)
loo_drift_compare <- loo_compare(loo_drift1, loo_drift2)
print(loo_drift_compare, simplify = FALSE)
