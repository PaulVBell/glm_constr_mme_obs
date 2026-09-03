
library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(bayesplot)

# ---- Configuration ---------------------------------------------------------

CSV      <- "data/exeter_annmax.csv"


END_HIST <- 2014


FAMILY   <- Gamma(link = "log")

CHAINS   <- 4
WARMUP   <- 1000
SAMPLING <- 1000
SEED     <- 1234

MAKE_PLOTS <- TRUE
FIG_DIR    <- "figures"

# ---- 1. Data ---------------------------------------------------------------

dat <- readr::read_csv(CSV, col_types = readr::cols(
  year = "i", member = "c", role = "c", precip = "d"
))

start_year <- min(dat$year)
end_year   <- max(dat$year)
nf         <- end_year - start_year + 1L


dat <- mutate(dat, t_sc = (year - start_year) / (nf - 1))

obs_df  <- filter(dat, role == "pseudo_obs", year <= END_HIST)
sim_df  <- filter(dat, role == "ensemble")
members <- sort(unique(sim_df$member))
m       <- length(members)


holdout <- filter(dat, role == "pseudo_obs", year > END_HIST)

# ---- 2. Stage one: a GLM per series ----------------------------------------

fits <- vector("list", m)
names(fits) <- members
for (j in seq_len(m)) {
  fits[[j]] <- glm(precip ~ t_sc, family = FAMILY,
                   data = filter(sim_df, member == members[j]))
}
a_fit <- glm(precip ~ t_sc, family = FAMILY, data = obs_df)

# ---- 3. Reduce each fit to an MLE and a Cholesky factor --------------------


len_a <- 2L
b_hat <- matrix(NA_real_, nrow = m, ncol = len_a)
J     <- array(0, dim = c(m, len_a, len_a))

for (j in seq_len(m)) {
  V <- vcov(fits[[j]])
  stopifnot(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
  b_hat[j, ] <- coef(fits[[j]])
  J[j, , ]   <- t(chol(V))
}

V_a   <- vcov(a_fit)

a_hat <- coef(a_fit)
J_a   <- t(chol(V_a))


# ---- 4. Model A: fixed delta = (1, 1) --------------------------------------

prior_hyper_param <- numeric()     

data_list <- list(
  len_a = len_a, m = m,
  b_hat = b_hat, a_hat = a_hat, J = J, J_a = J_a,
  len_hyper_param = length(prior_hyper_param),
  hyper_param = prior_hyper_param,
  delta_hyper_param = c(1, 1, 1, 1)
)

data_list_fixed_delta <- c(data_list, list(delta = c(1, 1)))

mod_fixed <- cmdstan_model("stan/glm_constr_mme_obs_fixed_delta.stan")
fit_fixed <- mod_fixed$sample(
  data = data_list_fixed_delta, seed = SEED,
  iter_warmup = WARMUP, iter_sampling = SAMPLING,
  chains = CHAINS, parallel_chains = CHAINS,
  refresh = 0, show_exceptions = FALSE
)

fit_fixed$summary(c("a", "a_ind", "mu", "tau"))

if (MAKE_PLOTS) {
  dir.create(FIG_DIR, showWarnings = FALSE)
  p_fixed <- fit_fixed$draws(c("a", "b_mme", "a_ind"), format = "df") |>
    pivot_longer(contains("["), names_to = c("param", "comp"),
                 names_pattern = "(.+)\\[(.)\\]") |>
    mutate(comp = c("Intercept", "Slope")[as.integer(comp)]) |>
    ggplot(aes(x = value, colour = param, fill = param)) +
    geom_density(alpha = .2) +
    scale_colour_manual(name = NULL, values = c(a = "blue", a_ind = "orange",
                                                b_mme = "black")) +
    scale_fill_manual(name = NULL, values = c(a = "blue", a_ind = "orange",
                                              b_mme = "black")) +
    facet_wrap(~comp, scales = "free_x") +
    labs(x = NULL, y = NULL, title = "Exeter, delta = (1, 1)")
  print(p_fixed)
  ggsave(file.path(FIG_DIR, "exeter-fixed-delta.png"), p_fixed,
         width = 8, height = 3.2, dpi = 150)
}

# ---- 5. Model B: a grid of fixed deltas ------------------------------------

deltas    <- c(0.001, seq(.25, 1, .25))
delta_mat <- cbind(rep(deltas, length(deltas)), rep(deltas, each = length(deltas)))


prior_hyper_param_B <- c(0, 0, 10, 5, 2, 2.5)

data_list_mult_delta <- list(
  len_a = len_a, m = m,
  b_hat = b_hat, a_hat = a_hat, J = J, J_a = J_a,
  len_hyper_param = length(prior_hyper_param_B),
  hyper_param = prior_hyper_param_B,
  n_delta = nrow(delta_mat),
  delta = delta_mat
)

mod_vec <- cmdstan_model("stan/glm_constr_mme_obs_fixed_delta_vector.stan")
fit_vec <- mod_vec$sample(
  data = data_list_mult_delta, seed = SEED,
  iter_warmup = WARMUP, iter_sampling = SAMPLING,
  chains = CHAINS, parallel_chains = CHAINS,
  refresh = 0, show_exceptions = FALSE
)


pct_change <- function(slope) 100 * (exp(slope) - 1)

if (MAKE_PLOTS) {
  p_sweep <- fit_vec$draws("a", format = "df") |>
    pivot_longer(contains("["), values_to = "a", names_to = c("index", "comp"),
                 names_pattern = "\\[(.+)\\,(.)\\]",
                 names_transform = list(index = as.integer, comp = as.integer)) |>
    mutate(delta_1 = delta_mat[index, 1], delta_2 = delta_mat[index, 2]) |>
    filter(comp == 2) |>
    mutate(pcnt = pct_change(a)) |>
    ggplot(aes(x = pcnt)) +
    geom_density() +
    xlim(-100, 200) +
    facet_grid(factor(delta_1, levels = rev(deltas)) ~ delta_2,
               labeller = label_both) +
    labs(x = "% change in annual maximum, 1981 to 2080", y = NULL,
         title = "Sensitivity to the similarity parameter")
  print(p_sweep)
  ggsave(file.path(FIG_DIR, "exeter-delta-sweep.png"), p_sweep,
         width = 8, height = 7, dpi = 150)
}

# ---- 6. Model C: a prior on delta ------------------------------------------


mod_prior <- cmdstan_model("stan/glm_constr_mme_obs.stan")
fit_prior <- mod_prior$sample(
  data = data_list, seed = SEED,
  iter_warmup = WARMUP, iter_sampling = SAMPLING,
  chains = CHAINS, parallel_chains = CHAINS,
  refresh = 0, show_exceptions = FALSE
)

fit_prior$summary(c("a", "a_ind", "delta", "mu"))

if (MAKE_PLOTS) {
  p_delta <- mcmc_areas(fit_prior$draws("delta")) +
    labs(title = "Posterior for the similarity parameter")
  print(p_delta)
  ggsave(file.path(FIG_DIR, "exeter-delta-posterior.png"), p_delta,
         width = 6, height = 3, dpi = 150)
}

# ---- 7. Scoring against the withheld member ---------------------------------


truth_fit   <- glm(precip ~ t_sc, family = FAMILY, data = filter(dat, role == "pseudo_obs"))
truth_slope <- unname(coef(truth_fit)[2])

slope_draws <- function(fit, var = "a") {
  d <- posterior::as_draws_matrix(fit$draws(var))
  as.numeric(d[, grep("\\[2\\]$", colnames(d))[1]])
}

comparison <- tibble::tibble(
  quantity = c("pseudo-obs alone (stage one)",
               "constrained, delta = (1,1)",
               "unconstrained control (a_ind)",
               "constrained, prior on delta",
               "WITHHELD TRUTH (member 01, full record)"),
  slope = c(unname(a_hat[2]),
            mean(slope_draws(fit_fixed, "a")),
            mean(slope_draws(fit_fixed, "a_ind")),
            mean(slope_draws(fit_prior, "a")),
            truth_slope),
  sd = c(sqrt(diag(V_a))[2],
         sd(slope_draws(fit_fixed, "a")),
         sd(slope_draws(fit_fixed, "a_ind")),
         sd(slope_draws(fit_prior, "a")),
         sqrt(diag(vcov(truth_fit)))[2])
) |>
  mutate(pct_change = pct_change(slope),
         covers_truth = abs(slope - truth_slope) < 1.645 * sd)

print(as.data.frame(comparison), digits = 3)

cat(sprintf("\nWithheld truth: slope %.3f  =>  %.1f%% change over 1981-2080\n",
            truth_slope, pct_change(truth_slope)))
