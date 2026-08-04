# Function to build different IPMs to play with

make_test_ipm <- function(type = "simple_di_det") {

  inv_logit <- function(int, slope, x) {
    1 / (1 + exp(-(int + slope * x)))
  }

  if (type == "simple_di_det") {

    data_list <- list(
      s_int = 2.2, s_slope = 0.25,
      g_int = 0.2, g_slope = 1.02, sd_g = 0.7,
      f_r_int = 0.003, f_r_slope = 0.015,
      f_s_int = 1.3, f_s_slope = 0.075,
      mu_fd = 2, sd_fd = 0.3
    )

    ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
      define_kernel(
        "P",
        formula = s * g,
        family = "CC",
        s = inv_logit(s_int, s_slope, dbh_1),
        g = dnorm(dbh_2, mu_g, sd_g),
        mu_g = g_int + g_slope * dbh_1,
        data_list = data_list,
        states = list(c("dbh")),
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "g")
      ) %>%
      define_kernel(
        "F",
        formula = f_r * f_s * f_d,
        family = "CC",
        f_r = inv_logit(f_r_int, f_r_slope, dbh_1),
        f_s = exp(f_s_int + f_s_slope * dbh_1),
        f_d = dnorm(dbh_2, mu_fd, sd_fd),
        data_list = data_list,
        states = list(c("dbh")),
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "f_d")
      ) %>%
      define_impl(
        make_impl_args_list(
          c("P", "F"),
          int_rule = rep("midpoint", 2),
          state_start = rep("dbh", 2),
          state_end = rep("dbh", 2)
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      define_pop_state(n_dbh = runif(100)) %>%
      make_ipm(
        usr_funs = list(inv_logit = inv_logit),
        iterate = TRUE,
        iterations = 50,
        normalize_pop_size = FALSE
      )

  } else if (type == "simple_dd_det") {

    data_list <- list(
      s_int = 2.2, s_slope = 0.25, dd_slope = -0.022,
      g_int = 0.2, g_slope = 1.02, sd_g = 0.7,
      f_r_int = 0.03, f_r_slope = 0.015,
      f_s_int = 1.3, f_s_slope = 0.075,
      mu_fd = 0.5, sd_fd = 0.2
    )

    ipm <- init_ipm(sim_gen = "simple", di_dd = "dd", det_stoch = "det") %>%
      define_kernel(
        "P",
        formula = s * g,
        family = "CC",
        s = plogis(s_int + s_slope * dbh_1 + dd_slope * sum(n_dbh_t)),
        g = dnorm(dbh_2, mu_g, sd_g),
        mu_g = g_int + g_slope * dbh_1,
        data_list = data_list,
        states = list(c("dbh")),
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "g")
      ) %>%
      define_kernel(
        "F",
        formula = f_r * f_s * f_d,
        family = "CC",
        f_r = plogis(f_r_int + f_r_slope * dbh_1),
        f_s = exp(f_s_int + f_s_slope * dbh_1 + dd_slope * sum(n_dbh_t)),
        f_d = dnorm(dbh_2, mu_fd, sd_fd),
        data_list = data_list,
        states = list(c("dbh")),
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "f_d")
      ) %>%
      define_impl(
        make_impl_args_list(
          c("P", "F"),
          int_rule = rep("midpoint", 2),
          state_start = rep("dbh", 2),
          state_end = rep("dbh", 2)
        )
      ) %>%
      define_domains(dbh = c(0, 50, 100)) %>%
      define_pop_state(n_dbh = runif(100)) %>%
      make_ipm(
        iterate = TRUE,
        iterations = 50,
        normalize_pop_size = FALSE
      )

  } else if (type == "simple_di_stoch") {

    inv_logit <- function(int, slope, x) {
      1 / (1 + exp(-(int + slope * x)))
    }

    mvt_wrapper <- function(r_means, r_sigma, nms) {
      out <- mvtnorm::rmvnorm(1, r_means, r_sigma)
      out <- as.list(out)
      names(out) <- nms
      out
    }

    data_list <- list(
      s_slope = 0.2,
      g_slope = 0.99,
      g_sd = 0.2,
      f_r_slope = 0.003,
      f_s_slope = 0.01,
      f_d_mu = 2,
      f_d_sd = 0.75
    )

    r_means <- c(
      s_int_yr = 0.8,
      g_int_yr = 0.1,
      f_r_int_yr = 0.3,
      f_s_int_yr = 0.01
    )

    r_sigma <- diag(4) * 0.1

    ipm <- init_ipm(
      sim_gen = "simple",
      di_dd = "di",
      det_stoch = "stoch",
      kern_param = "param"
    ) %>%
      define_kernel(
        "P",
        formula = s * g,
        family = "CC",
        g_mu = g_int_yr + g_slope * surf_area_1,
        s = inv_logit(s_int_yr, s_slope, surf_area_1),
        g = dnorm(surf_area_2, g_mu, g_sd),
        data_list = data_list,
        states = list(c("surf_area")),
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "g")
      ) %>%
      define_kernel(
        "F",
        formula = f_r * f_s * f_d,
        family = "CC",
        f_r = inv_logit(f_r_int_yr, f_r_slope, surf_area_1),
        f_s = exp(f_s_int_yr + f_s_slope * surf_area_1),
        f_d = dnorm(surf_area_2, f_d_mu, f_d_sd),
        data_list = data_list,
        states = list(c("surf_area")),
        evict_cor = TRUE,
        evict_fun = truncated_distributions("norm", "f_d")
      ) %>%
      define_impl(
        make_impl_args_list(
          c("P", "F"),
          int_rule = rep("midpoint", 2),
          state_start = rep("surf_area", 2),
          state_end = rep("surf_area", 2)
        )
      ) %>%
      define_domains(surf_area = c(0, 10, 100)) %>%
      define_env_state(
        env_params = mvt_wrapper(
          r_means, r_sigma,
          names(r_means)
        ),
        data_list = list(r_means = r_means, r_sigma = r_sigma)
      ) %>%
      define_pop_state(n_surf_area = runif(100)) %>%
      make_ipm(
        usr_funs = list(
          inv_logit = inv_logit,
          mvt_wrapper = mvt_wrapper
        ),
        iterate = TRUE,
        iterations = 50,
        normalize_pop_size = TRUE
      )

  } else {
    stop("Unknown type")
  }

  return(ipm)
}

ipm <- make_test_ipm("simple_di_det")
ipm_dd <- make_test_ipm("simple_dd_det")
ipm_stoch <- make_test_ipm("simple_di_stoch")

test_that("sensitivity errors on missing ipm", {
  expect_error(
    sensitivity(pars = NULL, kernels = NULL, delta = 1e-4),
    "`ipm` must be provided and cannot be NULL."
  )
})

test_that("sensitivity rejects non-list ipm", {
  expect_error(
    sensitivity(ipm = 5, pars = NULL, kernels = NULL),
    "must be a valid ipmr object"
  )
})

test_that("sensitivity runs on simple_di_det", {

  ipm <- make_test_ipm("simple_di_det")

  res <- sensitivity(ipm, kernels = c("P", "F"))

  expect_type(res, "list")
  expect_true(all(c("sensitivity", "elasticity") %in% names(res)))

  expect_true(is.numeric(res$sensitivity))
  expect_true(is.numeric(res$elasticity))
})

test_that("sensitivity respects parameter subset", {

  ipm <- make_test_ipm("simple_di_det")

  res <- sensitivity(
    ipm,
    pars = c("s_int", "g_int"),
    kernels = c("P", "F")
  )

  expect_equal(names(res$sensitivity), c("s_int", "g_int"))
  expect_equal(length(res$sensitivity), 2)
})

test_that("sensitivity rejects density-dependent IPMs", {

  ipm <- make_test_ipm("simple_dd_det")

  expect_error(
    sensitivity(ipm, kernels = c("P", "F")),
    "Only density-independent IPMs are currently supported."
  )
})

test_that("sensitivity rejects stochastic IPMs", {

  ipm <- make_test_ipm("simple_di_stoch")

  expect_error(
    sensitivity(ipm, kernels = c("P", "F")),
    "Only deterministic IPMs are currently supported."
  )
})

test_that("sensitivity errors on unknown parameter", {

  ipm <- make_test_ipm("simple_di_det")

  expect_error(
    sensitivity(ipm, pars = "not_a_param", kernels = c("P", "F")),
    "Unknown parameter"
  )
})

test_that("elasticity is consistent with sensitivity definition", {

  ipm <- make_test_ipm("simple_di_det")

  res <- sensitivity(ipm, kernels = c("P", "F"))

  pars_all <- parameters(ipm)
  lambda <- {
    ker <- make_iter_kernel(ipm, mega_mat = c("P", "F"))
    Re(eigen(ker[[1]])$values[1])
  }

  expected_elas <- res$sensitivity * abs(as.numeric(pars_all[names(res$sensitivity)]) / lambda)

  expect_equal(
    res$elasticity,
    expected_elas,
    tolerance = 1e-6
  )
})

test_that("duplicate kernels trigger warning", {

  ipm <- make_test_ipm("simple_di_det")

  expect_warning(
    sensitivity(ipm, kernels = c("P", "P")),
    "Duplicate kernel name"
  )
})

test_that("kernel length mismatch triggers warning", {

  ipm <- make_test_ipm("simple_di_det")

  expect_warning(
    sensitivity(ipm, kernels = c("P")),
    "has length"
  )
})

test_that("unknown kernel errors clearly", {

  ipm <- make_test_ipm("simple_di_det")

  expect_error(
    sensitivity(ipm, kernels = c("P", "NOT_REAL")),
    "Unknown kernel"
  )
})
