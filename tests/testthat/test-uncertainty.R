# Function to build different IPMs to play with --------------------------------

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

samples <- as.data.frame(
  lapply(parameters(ipm), function(x) rep(x, 2))
)

vr <- data.frame(
  parameter = names(parameters(ipm)),
  vital_rate = "test"
)

ker <- c("P", "F")

# Tests for structural correctness ---------------------------------------------

test_that("uncertainty errors on missing IPM", {

  expect_error(
    uncertainty(pars = NULL, samples = NULL, mega_mat = NULL, vr_table = NULL),
    "`ipm` must be provided"
  )
})

test_that("uncertainty rejects non-list IPM", {

  expect_error(
    uncertainty(ipm = 5, samples = data.frame(), mega_mat = NULL, vr_table = NULL),
    "valid ipmr object"
  )
})

test_that("uncertainty rejects DD IPMs", {

  ipm <- make_test_ipm("simple_dd_det")

  expect_error(
    uncertainty(ipm, samples = data.frame(), mega_mat = c("P","F"), vr_table = NULL),
    "density-independent"
  )
})

test_that("uncertainty rejects stochastic IPMs", {

  ipm <- make_test_ipm("simple_di_stoch")

  expect_error(
    uncertainty(ipm, samples = data.frame(), mega_mat = c("P","F"), vr_table = NULL),
    "deterministic"
  )
})

test_that("uncertainty returns correct structure", {

  ipm <- make_test_ipm("simple_di_det")

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P","F"),
    vr_table = vr
  )

  expect_s3_class(res, "ipmr_uncertainty")

  expect_true(all(c("lambdas", "mod_uncert", "params_uncert") %in%
                    names(res)))
})

test_that("uncertainty respects parameter subset", {

  ipm <- make_test_ipm("simple_di_det")

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    pars = c("s_int", "g_int"),
    samples = samples,
    mega_mat = c("P","F"),
    vr_table = vr
  )

  expect_true(all(res$params_uncert$parameter %in%
                    c("s_int", "g_int")))
})

test_that("uncertainty rejects non-list bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      vr_table = vr,
      bounds = c(0, 1)
    ),
    "`bounds` must be NULL or a named list"
  )

})

test_that("uncertainty rejects unnamed bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(c(0, 1)),
      vr_table = vr
    ),
    "must be a named list"
  )

})

test_that("uncertainty rejects unknown parameter names in bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        fake_parameter = c(0, 1)
      ),
      vr_table = vr
    ),
    "Unknown parameter"
  )

})

test_that("uncertainty rejects bounds of incorrect length 1", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c(0)
      ),
      vr_table = vr
    ),
    "length 2"
  )

})

test_that("uncertainty rejects bounds of incorrect length 3", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c(0, 1, 2)
      ),
      vr_table = vr
    ),
    "length 2"
  )

})

test_that("uncertainty rejects non-numeric bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c("a", "b")
      ),
      vr_table = vr
    ),
    "must be numeric"
  )

})

test_that("uncertainty rejects missing values in bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c(NA, 1)
      ),
      vr_table = vr
    ),
    "may not contain NA"
  )

})

test_that("uncertainty rejects reversed bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c(1, 0)
      ),
      vr_table = vr
    ),
    "lower bound"
  )

})

test_that("uncertainty rejects parameter values outside supplied bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_error(
    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c(10, 20)
      ),
      vr_table = vr
    ),
    "outside the supplied bounds"
  )

})

test_that("uncertainty accepts infinite bounds", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  expect_no_error(

    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      bounds = list(
        s_int = c(-Inf, Inf),
        f_s_int = c(0, Inf)
      ),
      vr_table = vr
    )

  )

})

test_that("uncertainty warns when delta exceeds feasible region", {

  ipm <- make_test_ipm("simple_di_det")
  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"
  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  pars <- parameters(ipm)

  expect_error(

    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      delta = 5,
      bounds = list(
        s_int = c(pars$s_int - 1,
                  pars$s_int + 1)
      ),
      vr_table = vr
    ),

    "delta"

  )

})

test_that("uncertainty returns difference method", {

  ipm <- make_test_ipm("simple_di_det")

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 10)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P", "F"),
    vr_table = vr
  )

  expect_true("difference_method" %in% names(res$params_uncert))
})

test_that("uncertainty records forward differences for bounded parameters", {

  ipm <- make_test_ipm("simple_di_det")

  pars <- parameters(ipm)

  samples <- as.data.frame(
    lapply(pars, function(x) rep(x, 10))
  )

  samples$s_int[] <- pars$s_int

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P", "F"),
    vr_table = vr,
    bounds = list(s_int = c(pars$s_int, Inf))
  )

  expect_equal(
    res$params_uncert$difference_method[
      res$params_uncert$parameter == "s_int"
    ],
    "forward"
  )

})

test_that("uncertainty records backward differences for bounded parameters", {

  ipm <- make_test_ipm("simple_di_det")

  pars <- parameters(ipm)

  samples <- as.data.frame(
    lapply(pars, function(x) rep(x, 10))
  )

  samples$s_int[] <- pars$s_int

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P", "F"),
    vr_table = vr,
    bounds = list(s_int = c(-Inf, pars$s_int))
  )

  expect_equal(
    res$params_uncert$difference_method[
      res$params_uncert$parameter == "s_int"
    ],
    "backward"
  )

})

test_that("uncertainty records mixed difference methods", {

  ipm <- make_test_ipm("simple_di_det")

  pars <- parameters(ipm)

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 10)))

  ## Half of the samples sit on the lower bound, half do not
  samples$s_int <- c(rep(pars$s_int, 5), rep(pars$s_int + 0.5, 5))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P", "F"),
    vr_table = vr,
    bounds = list(
      s_int = c(pars$s_int, Inf)
    )
  )

  expect_match(
    res$params_uncert$difference_method[
      res$params_uncert$parameter == "s_int"
    ],
    "^mixed"
  )

})

test_that("uncertainty errors when delta exceeds feasible interval", {

  ipm <- make_test_ipm("simple_di_det")

  pars <- parameters(ipm)

  samples <- as.data.frame(
    lapply(pars, function(x) rep(x, 5))
  )

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  expect_error(

    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      vr_table = vr,
      delta = 1,
      bounds = list(
        s_int = c(
          pars$s_int - 0.25,
          pars$s_int + 0.25
        )
      )
    ),

    "Choose a smaller value of delta"

  )

})

test_that("uncertainty errors when sampled values cross parameter bounds", {

  ipm <- make_test_ipm("simple_di_det")

  pars <- parameters(ipm)

  samples <- as.data.frame(
    lapply(pars, function(x) rep(x, 10))
  )

  samples$s_int <- c(
    rep(pars$s_int, 5),
    rep(pars$s_int - 0.05, 5)
  )

  vr <- data.frame(
    parameter = names(parameters(ipm)),
    vital_rate = "test" )

  expect_error(

    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P","F"),
      vr_table = vr,
      delta = 0.1,
      bounds = list(
        s_int = c(pars$s_int, Inf)
      )
    ),

    "supplied bounds"
  )

})

test_that("uncertainty rejects sampled values outside bounds", {

  pars <- parameters(ipm)

  samples <- as.data.frame(
    lapply(pars, function(x) rep(x, 10))
  )

  samples$s_int[1] <- pars$s_int - 0.5

  expect_error(

    uncertainty(
      ipm,
      samples = samples,
      mega_mat = c("P", "F"),
      vr_table = vr,
      bounds = list(
        s_int = c(pars$s_int, Inf)
      )
    ),

    "Sampled value\\(s\\) fall outside the supplied bounds"

  )

})

test_that("uncertainty template mode returns expected table", {

  ipm <- make_test_ipm("simple_di_det")

  vr <- uncertainty(ipm, vr_table = "template")

  expect_s3_class(vr, "data.frame")

  expect_equal(
    names(vr),
    c("parameter", "vital_rate")
  )

  expect_equal(
    vr$parameter,
    names(parameters(ipm))
  )

  expect_true(all(is.na(vr$vital_rate)))

})



# Tests for numerical consistency ----------------------------------------------

test_that("uncertainty computes lambda variance", {

  ipm <- make_test_ipm("simple_di_det")

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P","F"),
    vr_table = vr
  )

  expect_true(is.numeric(res$mod_uncert))
  expect_true(res$mod_uncert >= 0)
})

test_that("uncertainty sensitivities match sensitivity()", {

  ipm <- make_test_ipm("simple_di_det")

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  unc <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P","F"),
    vr_table = vr
  )

  sens <- sensitivity(ipm, mega_mat = c("P","F"))$sensitivity

  expect_equal(
    unname(unc$params_uncert$sensitivity),
    unname(sens[unc$params_uncert$parameter]),
    tolerance = 1e-6
  )
})

test_that("bounded uncertainty sensitivities match bounded sensitivity()", {

  ipm <- make_test_ipm("simple_di_det")

  pars <- parameters(ipm)

  samples <- as.data.frame(lapply(pars, function(x) rep(x, 25)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "test"

  bounds <- list(
    s_int = c(pars$s_int, Inf)
  )

  unc <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P", "F"),
    vr_table = vr,
    bounds = bounds
  )

  sens <- sensitivity(
    ipm,
    mega_mat = c("P", "F"),
    bounds = bounds
  )

  expect_equal(
    unname(unc$params_uncert$sensitivity),
    unname(sens$sensitivity[unc$params_uncert$parameter]),
    tolerance = 1e-6
  )

})

test_that("vr_table join works correctly", {

  ipm <- make_test_ipm("simple_di_det")

  samples <- as.data.frame(lapply(parameters(ipm), function(x) rep(x, 50)))

  vr <- uncertainty(ipm, vr_table = "template")
  vr$vital_rate <- "growth"

  res <- uncertainty(
    ipm,
    samples = samples,
    mega_mat = c("P","F"),
    vr_table = vr
  )

  expect_true("vital_rate" %in% names(res$params_uncert))
})
