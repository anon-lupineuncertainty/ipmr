#' @title Sensitivity and elasticity analysis of population growth rate
#'
#' @description Computes numerical sensitivities and elasticities of the
#' dominant eigenvalue (\eqn{\lambda}) with respect to model parameters using
#' finite differences.
#'
#' @param ipm An \code{ipmr} IPM object.
#' @param pars A character vector of parameter names to perturb. Defaults to
#'  all parameters in \code{ipm}.
#' @param kernels A character vector specifying the kernel structure in row
#'  major order.
#' @param delta A numeric scalar specifying the perturbation size used for
#'  numerical differentiation. Default is 1e-4.
#'  @param bounds An optional named list specifying lower and/or upper bounds
#'   for model parameters when calculating numerical sensitivities. Each element
#'   should be a numeric vector of length two giving the lower and upper bounds,
#'   respectively (e.g., \code{list(surv = c(0, 1), recruit = c(0, Inf))}).
#'   Parameters omitted from this list are assumed to be unbounded and are
#'   perturbed using the central finite difference. For parameters with
#'   specified bounds, the function automatically uses a forward or backward
#'   finite difference whenever perturbation by \code{delta} would cross a
#'   boundary; otherwise, the central finite difference is used.
#'
#'@details
#' This function uses perturbation to estimate the sensitivity of \eqn{\lambda}
#' to each parameter in \code{pars}.
#'
#' Currently implemented for deterministic ("det"), density-independent
#' ("di") IPMs constructed with \code{ipmr}. Stochastic ("stoch") and
#' density-dependent ("dd") IPMs are not supported.
#'
#' ## Finite differences and bounded parameters
#'
#' Sensitivities are estimated numerically by perturbing each parameter by
#' \code{delta}. By default, a central finite difference is used for all
#' parameters. However, some demographic parameters are biologically constrained
#' (e.g., probabilities bounded between 0 and 1, or fecundity parameters bounded
#' below by 0). If perturbing a parameter by \code{delta} would cross one of
#' these user-specified bounds, the function automatically computes the
#' appropriate one-sided finite difference (forward or backward) rather than
#' evaluating the model outside the feasible parameter space.
#'
#' One-sided finite differences approximate the same local derivative as the
#' central finite difference, but are used only when a symmetric perturbation
#' would leave the biologically feasible parameter space. Consequently,
#' sensitivities calculated using one-sided and central finite differences are
#' directly comparable and may be used together in downstream uncertainty
#' analyses.
#'
#' Supplying \code{bounds} is recommended whenever model parameters are subject
#' to biological constraints and numerical perturbation could otherwise produce
#' invalid parameter values or biologically impossible kernels.
#'
#' @return A named list with three elements:
#' \itemize{
#'   \item \code{sensitivity}: A named numeric vector giving the numerical
#'     sensitivities of the dominant eigenvalue (\eqn{\lambda}) to each
#'     parameter.
#'   \item \code{elasticity}: A named numeric vector giving the elasticities
#'     of the dominant eigenvalue (\eqn{\lambda}) to each parameter.
#'   \item \code{difference_method}: A named character vector indicating the
#'     finite-difference scheme used to calculate the sensitivity for each
#'     parameter. This component is provided to facilitate reproducibility and
#'     to identify parameters whose sensitivities were evaluated using one-sided
#'     finite differences because of user-specified parameter bounds. Values are
#'     \code{"central"}, \code{"forward"}, or \code{"backward"}.
#' }
#'
#' @examples
#' # Probability bounded between 0 and 1
#' # bounds <- list(
#' #   surv = c(0, 1)
#' # )
#'
#' # Recruitment bounded below by zero
#' # bounds <- list(
#' #   recruit = c(0, Inf)
#' # )
#'
#' # Multiple constrained parameters
#' # bounds <- list(
#' #   surv = c(0, 1),
#' #   recruit = c(0, Inf),
#' #   fec = c(0, Inf)
#' # )
#'
#' # Compute sensitivities using defined parameter boundaries
#' # sens <- sensitivity(
#' #   ipm,
#' #   kernels = c("P", "F"),
#' #   bounds = bounds )
#'
#' @export
sensitivity <- function(ipm, pars = NULL, kernels, delta = 1e-4,
                        bounds = NULL) {

  ## ---- Input checks

  val <- validate_ipm_sensitivity(ipm, pars, kernels, delta, bounds)

  ipm      <- val$ipm
  pars     <- val$pars
  pars_all <- val$pars_all
  kernels  <- val$kernels
  delta    <- val$delta
  bounds   <- val$bounds

  ## ---- Setup

  npar        <- length(pars)
  spar        <- numeric(npar)
  diff_method <- character(npar)

  kernel <- make_iter_kernel(ipm, mega_mat = kernels)

  if (!is.list(kernel) || is.null(kernel[[1]])) {
    stop("Kernel construction failed.")
  }

  lambda_all <- dominant_eigen(kernel[[1]])

  ## ---- Helper to perturb a single parameter and calculate lambda

  perturb_lambda <- function(par, value) {
    vpar <- pars_all
    vpar[[par]] <- value
    ipm_tmp <- ipm$proto_ipm
    parameters(ipm_tmp) <- vpar
    ker <- ipm_tmp %>%
      make_ipm(iterate = FALSE) %>%
      make_iter_kernel(mega_mat = kernels)
    dominant_eigen(ker[[1]])
  }

  ## ---- Loop over each parameter

  for (i in seq_len(npar)) {

    par_now <- pars[i]
    theta <- pars_all[[par_now]]

    ## ---- Define bounds (default none)

    if (!is.null(bounds) && par_now %in% names(bounds)) {
      lower <- bounds[[par_now]][1]
      upper <- bounds[[par_now]][2]
    } else {
      lower <- -Inf
      upper <- Inf
    }

    use_forward <- theta - delta < lower
    use_backward <- theta + delta > upper

    ## ---- Forward difference

    if (use_forward) {
      lambda_up <- perturb_lambda(par_now, theta + delta)

      spar[i] <- (lambda_up - lambda_all)/delta
      diff_method[i] <- "forward"

      ## ---- Backward difference
    } else if (use_backward) {
      lambda_down <- perturb_lambda(par_now, theta - delta)

      spar[i] <- (lambda_all - lambda_down)/delta
      diff_method[i] <- "backward"

      ## ---- Central difference
    } else {
      lambda_down <- perturb_lambda(par_now, theta - delta)
      lambda_up <- perturb_lambda(par_now, theta + delta)

      spar[i] <- (lambda_up - lambda_down)/(2*delta)
      diff_method[i] <- "central"
    }
  }

  ## ---- Output

  sens <- spar
  elas <- spar * abs(as.numeric(pars_all[pars]) / lambda_all)

  names(sens) <- pars
  names(elas) <- pars
  names(diff_method) <- pars

  return(list(
    sensitivity = sens,
    elasticity = elas,
    difference_method = diff_method
  ))
}
