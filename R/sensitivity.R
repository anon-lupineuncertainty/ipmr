#' @title Sensitivity and elasticity analysis of population growth rate
#'
#' @description Computes numerical sensitivities and elasticities of the
#' dominant eigenvalue (\eqn{\lambda}) with respect to model parameters.
#'
#' @param ipm An \code{ipmr} IPM object.
#' @param pars A character vector of parameter names to perturb. Defaults to
#'  all parameters in \code{ipm}.
#' @param kernels A character vector specifying the kernel structure in row
#'  major order.
#' @param delta A numeric scalar specifying the perturbation size used for
#'  numerical differentiation. Default is 1e-4.
#'
#'@details
#' This function uses a central difference approximation to estimate the
#' sensitivity of \eqn{\lambda} to each parameter in \code{pars}.
#'
#' Currently implemented for deterministic ("det"), density-independent
#' ("di") IPMs constructed with \code{ipmr}. Stochastic ("stoch") and
#' density-dependent ("dd") IPMs are not supported.
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item \code{sensitivity}: Numeric vector of sensitivities of
#'     \eqn{\lambda} to each parameter.
#'   \item \code{elasticity}: Numeric vector of elasticities of
#'     \eqn{\lambda} to each parameter.
#' }
#'
#' @export
sensitivity <- function(ipm, pars = NULL, kernels, delta = 1e-4) {

  ## ---- Input checks

  val <- validate_ipm_sensitivity(ipm, pars, kernels, delta)

  ipm      <- val$ipm
  pars     <- val$pars
  pars_all <- val$pars_all
  kernels  <- val$kernels
  delta    <- val$delta

  ## ---- Setup

  npar <- length(pars)
  spar <- numeric(npar)
  vpar <- pars_all

  kernel <- make_iter_kernel(ipm, mega_mat = kernels)

  if (!is.list(kernel) || is.null(kernel[[1]])) {
    stop("Kernel construction failed.")
  }

  eig <- tryCatch(
    eigen(kernel[[1]]),
    error = function(e) stop("Eigen decomposition failed for base kernel.")
  )

  lambda_all <- dominant_eigen(kernel[[1]])

  ## ---- Loop over each parameter

  for(i in seq_len(npar)) {

    par_now <- pars[i]

    # ---- DOWN perturbation
    vpar <- pars_all
    vpar[[par_now]] <- vpar[[par_now]] - delta

    ipm_down <- ipm$proto_ipm
    parameters(ipm_down) <- vpar

    ker_down <- ipm_down %>%
      make_ipm(iterate = FALSE) %>%
      make_iter_kernel(mega_mat = kernels)

    lambda_down <- dominant_eigen(ker_down[[1]])

    # ---- UP perturbation
    vpar <- pars_all
    vpar[[par_now]] <- vpar[[par_now]] + delta

    ipm_up <- ipm$proto_ipm
    parameters(ipm_up) <- vpar

    ker_up <- ipm_up %>%
      make_ipm(iterate = FALSE) %>%
      make_iter_kernel(mega_mat = kernels)

    lambda_up <- dominant_eigen(ker_up[[1]])

    # ---- Central difference
    spar[i] <- (lambda_up - lambda_down) / (2 * delta)
  }

  ## ---- Output

  sens <- spar
  elas <- spar * abs(as.numeric(pars_all[pars]) / lambda_all)

  names(sens) <- pars
  names(elas) <- pars

  return(list(
    sensitivity = sens,
    elasticity = elas
  ))
}
