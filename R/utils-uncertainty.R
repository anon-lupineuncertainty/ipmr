#' Extract dominant eigenvalue
#'
#' Returns the dominant eigenvalue using the same ordering assumption
#' currently used in sensitivity calculations.
#'
#' @noRd
dominant_eigen <- function(mat) {
  eig <- eigen(mat)
  Re(eig$values[which.max(Mod(eig$values))])
}


#' Validate inputs for sensitivity and uncertainty functions
#'
#' Returns usable objects to be passed to the function-specific validations.
#'
#' @noRd
validate_ipm_base <- function(ipm) {

  if (missing(ipm) || is.null(ipm)) {
    stop("`ipm` must be provided and cannot be NULL.")
  }

  if (!is.list(ipm)) {
    stop("`ipm` must be a valid ipmr object.")
  }

  cls <- strsplit(class(ipm)[1], "_")[[1]]

  if (length(cls) < 3) {
    stop("Malformed IPM class structure.")
  }

  di_dd <- cls[2]
  det_stoch <- cls[3]

  if (is.na(di_dd) || di_dd != "di") {
    stop("Only density-independent IPMs are currently supported.")
  }

  if (is.na(det_stoch) || det_stoch != "det") {
    stop("Only deterministic IPMs are currently supported.")
  }

  pars_all <- tryCatch(
    parameters(ipm),
    error = function(e) {
      stop("Could not extract parameters from IPM.")
    }
  )

  if (!is.list(pars_all) || is.null(names(pars_all))) {
    stop("`parameters(ipm)` must return a named list.")
  }

  list(
    ipm = ipm,
    pars_all = pars_all
  )
}

#' Warn if finite-difference perturbation may change parameter sign
#'
#' Internal helper used by sensitivity and uncertainty validation.
#'
#' @param par_names Character vector of parameter names.
#' @param par_values Numeric vector of parameter values corresponding to
#'   `par_names`.
#' @param delta Finite-difference perturbation size.
#'
#' @noRd
warn_delta_crossing <- function(par_names, par_values = NULL, delta) {

  if (!is.null(par_values)) {
    par_names <- unique(par_names[abs(par_values) <= delta])
  }

  if (length(par_names) > 0) {
    warning(
      "The following parameter(s) are within `delta` of zero and may change ",
      "sign during perturbation: ",
      paste(par_names, collapse = ", "),
      ". This may result in failed or invalid subkernel construction or failed ",
      "sensitivity/uncertainty calculations."
    )
  }

  invisible(NULL)
}



#' Validate inputs to sensitivity function
#'
#' Returns usable objects for the sensitivity function.
#'
#' @noRd
validate_ipm_sensitivity <- function(ipm, pars, kernels, delta) {

  base <- validate_ipm_base(ipm)

  pars_all <- base$pars_all
  ipm <- base$ipm

  if (is.null(names(pars_all))) {
    stop("Parameters must be named.")
  }

  if (is.null(pars)) {
    pars <- names(pars_all)
  } else {
    if (!is.character(pars)) {
      stop("`pars` must be a character vector of parameter names.")
    }

    if (length(pars) == 0) {
      stop("`pars` cannot be empty.")
    }

    if (!all(pars %in% names(pars_all))) {
      bad <- pars[!pars %in% names(pars_all)]
      stop("Unknown parameter(s): ", paste(bad, collapse = ", "))
    }

    pars <- unique(pars)
  }

  if (any(pars_all[pars] == 0)) {
    warning("Zero-valued parameters detected; elasticity may be undefined.")
  }

  ## ---- kernels

  if (missing(kernels) || is.null(kernels)) {
    stop("`kernels` must be provided.")
  }

  if (!is.character(kernels) || length(kernels) == 0) {
    stop("`kernels` must be a non-empty character vector.")
  }

  # Check for duplicates
  dup_kerns <- unique(kernels[duplicated(kernels)])

  if (length(dup_kerns) > 0) {
    warning(
      "Duplicate kernel name(s) detected in `kernels`: ",
      paste(dup_kerns, collapse = ", "),
      ". This may lead to incorrect kernel construction if not intentional."
    )
  }

  available_kernels <- names(ipm$sub_kernels)

  if (!all(kernels %in% available_kernels)) {
    bad <- setdiff(kernels, available_kernels)
    stop(
      "Unknown kernel name(s): ",
      paste(bad, collapse = ", "),
      ". Available kernels are: ",
      paste(available_kernels, collapse = ", ")
    )
  }

  n_expected <- length(ipm$sub_kernels)

  if (length(kernels) != n_expected) {
    warning(
      "`kernels` has length ", length(kernels),
      " but IPM contains ", n_expected, " subkernels. ",
      "Ensure kernels are provided in correct row-major order."
    )
  }

  ## ---- delta

  if (!is.numeric(delta) || length(delta) != 1 || is.na(delta)) {
    stop("`delta` must be a single numeric value.")
  }

  if (delta <= 0) {
    stop("`delta` must be > 0.")
  }

  ## ---- check if delta changes sign on any parameters

  warn_delta_crossing(
    par_names = pars,
    par_values = pars_all[pars],
    delta = delta
  )

  ## return cleaned objects for main function
  list(
    ipm = ipm,
    pars = pars,
    pars_all = pars_all,
    kernels = kernels,
    delta = delta
  )
}


#' Create a new ipmr uncertainty object
#'
#' Internal constructor for ipmr uncertainty outputs.
#'
#' @param ipm The original IPM
#' @param lambdas Data frame of lambda values
#' @param mod_uncert Numeric variance in lambda
#' @param params_uncert Data frame of parameter summaries
#' @param vr_uncert Data frame of uncertainty contributions summed by vital rate
#'
#' @noRd
new_ipmr_uncertainty <- function(ipm, lambdas, mod_uncert, params_uncert, vr_uncert) {

  structure(
    list(
      ipm           = ipm,
      lambdas       = lambdas,
      mod_uncert    = mod_uncert,
      params_uncert = params_uncert,
      vr_uncert     = vr_uncert
    ),
    class = "ipmr_uncertainty"
  )
}

#' Validate inputs to uncertainty function
#'
#'
#' @noRd
validate_ipm_uncertainty <- function(ipm, pars, samples, kernels, vr_table,
                                     delta, cores) {

  base <- validate_ipm_base(ipm)

  pars_all <- base$pars_all
  ipm <- base$ipm

  if (is.null(pars)) {
    pars <- names(pars_all)
  }

  if (!is.character(pars)) {
    stop("`pars` must be a character vector.")
  }

  if (!all(pars %in% names(pars_all))) {
    stop("Unknown parameter(s): ",
         paste(pars[!pars %in% names(pars_all)], collapse = ", "))
  }

  ## ---- samples
  if (missing(samples) || !is.data.frame(samples)) {
    stop("`samples` must be a data.frame.")
  }

  if (nrow(samples) < 1) {
    stop("`samples` must have at least one row. We suggest using no fewer than 100 samples.")
  }

  missing_in_samples <- setdiff(names(pars_all), colnames(samples))
  if (length(missing_in_samples) > 0) {
    stop("Samples missing required parameter columns: ",
         paste(missing_in_samples, collapse = ", "))
  }

  extra_in_samples <- setdiff(colnames(samples), names(pars_all))
  if (length(extra_in_samples) > 0) {
    warning("Extra columns in `samples` not used: ",
            paste(extra_in_samples, collapse = ", "))
  }

  ## ---- kernels
  if (missing(kernels) || is.null(kernels)) {
    stop("`kernels` must be provided.")
  }

  if (!is.character(kernels) || length(kernels) == 0) {
    stop("`kernels` must be a non-empty character vector.")
  }

  # Check for duplicates
  dup_kerns <- unique(kernels[duplicated(kernels)])

  if (length(dup_kerns) > 0) {
    warning(
      "Duplicate kernel name(s) detected in `kernels`: ",
      paste(dup_kerns, collapse = ", "),
      ". This may lead to incorrect kernel construction if not intentional."
    )
  }

  available_kernels <- names(ipm$sub_kernels)

  if (!all(kernels %in% available_kernels)) {
    bad <- setdiff(kernels, available_kernels)
    stop(
      "Unknown kernel name(s): ",
      paste(bad, collapse = ", "),
      ". Available kernels are: ",
      paste(available_kernels, collapse = ", ")
    )
  }

  n_expected <- length(ipm$sub_kernels)

  if (length(kernels) != n_expected) {
    warning(
      "`kernels` has length ", length(kernels),
      " but IPM contains ", n_expected, " subkernels. ",
      "Ensure kernels are provided in correct row-major order."
    )
  }

  ## ---- vr_table
  if (!identical(vr_table, "template")) {

    if (!is.data.frame(vr_table)) {
      stop("`vr_table` must be a data.frame or 'template'.")
    }

    if (!all(c("parameter", "vital_rate") %in% colnames(vr_table))) {
      stop("`vr_table` must contain columns: parameter, vital_rate")
    }

    if (any(duplicated(vr_table$parameter))) {
      stop("Each parameter must appear only once in `vr_table`.")
    }

    if (!all(pars %in% vr_table$parameter)) {
      stop("Missing parameters in `vr_table`.")
    }
  }

  ## ---- delta
  if (!is.numeric(delta) || length(delta) != 1 || is.na(delta) || delta <= 0) {
    stop("`delta` must be a single numeric value > 0.")
  }

  ## ---- check if delta changes sign on any parameters
  problem_pars <- pars[
    vapply(
      pars,
      function(p) any(abs(samples[[p]]) <= delta),
      logical(1)
    )
  ]

  warn_delta_crossing(
    par_names = problem_pars,
    delta = delta
  )

  ## ---- cores
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1) {
    stop("`cores` must be a single integer >= 1.")
  }

  cores <- as.integer(cores)

  if (cores < 1) {
    stop("`cores` must be >= 1.")
  }

  max_cores <- parallel::detectCores(logical = TRUE)
  if (cores > max_cores) {
    stop(
      "`cores` (", cores, ") exceeds available logical cores (",
      max_cores, ")."
    )
  }

  slurm <- Sys.getenv("SLURM_JOB_ID") != ""
  pbs   <- Sys.getenv("PBS_JOBID") != ""

  if ((slurm || pbs) && cores == max_cores) {
    warning(
      "Detected HPC scheduler environment (Slurm/PBS). ",
      "Using all detected cores may exceed allocated resources."
    )
  }

  list(
    ipm      = ipm,
    pars     = pars,
    samples  = samples,
    kernels  = kernels,
    vr_table = vr_table,
    delta    = delta,
    cores    = as.integer(cores)
  )
}

#' Plot three-panel uncertainty plot
#'
#' @noRd
.plot_uncert_param <- function(x, vr_colors, ...) {

  df <- x$params_uncert

  # Order parameters by contribution
  df$parameter <- factor(
    df$parameter,
    levels = rev(df$parameter[order(df$vital_rate)])
  )

  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5)
    )

  title_wrap <- function(x) stringr::str_wrap(x, width = 28)

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cv, y = .data$parameter,
                                         color = .data$vital_rate)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$cv,
                                       yend = .data$parameter), linewidth = 1) +
    ggplot2::labs(x = "Coefficient of variation",
                  y = "Parameter",
                  title = ggplot2::title_wrap("Parameter uncertainty")) +
    ggplot2::scale_color_manual(values = vr_colors, guide = "none") +
    base_theme

  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$elasticity, y = .data$parameter,
                                         color = .data$vital_rate)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$elasticity,
                                       yend = .data$parameter), linewidth = 1) +
    ggplot2::labs(x = "Elasticity",
                  y = NULL,
                  title = title_wrap("Model sensitivity to parameter")) +
    ggplot2::scale_color_manual(values = vr_colors, guide = "none") +
    base_theme

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$variance_prop, y = .data$parameter,
                                         color = .data$vital_rate)) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$variance_prop,
                                       yend = .data$parameter), linewidth = 1) +
    ggplot2::labs(x = "Uncertainty contribution",
                  y = NULL,
                  title = title_wrap("Parameter contribution to uncertainty in \u03BB"),
                  color = "Vital rate") +
    ggplot2::scale_color_manual(values = vr_colors) +
    base_theme

  p1 + p2 + p3 + patchwork::plot_layout(ncol = 3, axes = 'collect_y')
}

#' Plot uncertainty aggregated by vital rate
#'
#' @noRd
.plot_uncert_vr <- function(x, vr_colors, ...) {

  df <- x$vr_uncert

  df$vital_rate <- factor(
    df$vital_rate,
    levels = df$vital_rate[order(df$variance_sum)]
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$variance_sum, y = .data$vital_rate,
                                   fill = .data$vital_rate)) +
    ggplot2::geom_col() +
    ggplot2::labs( x = "Total uncertainty contribution",
                   y = "Vital rate" ) +
    ggplot2::scale_fill_manual(values = vr_colors) +
    ggplot2::theme_minimal() +
    ggplot2::guides(fill = "none")
}

#' Plot model uncertainty and vital rate uncertainty contributions
#'
#' @noRd
.plot_uncert_stacked <- function(x, vr_colors, ...) {

  df <- x$vr_uncert

  total_var <- x$mod_uncert


  p <- ggplot2::ggplot(df, ggplot2::aes(x = 1, y = .data$variance_sum,
                                        fill = .data$vital_rate)) +
    ggplot2::geom_col(width = 0.6) +
    # Outline box for total lambda variance
    ggplot2::geom_rect(
      ggplot2::aes( xmin = 0.7,
                    xmax = 1.3,
                    ymin = 0,
                    ymax = total_var,
                    color = "Total variance in lambda"),
      inherit.aes = FALSE,
      fill = NA,
      linewidth = 1 ) +
    ggplot2::labs( x = NULL,
                   y = "Variance in \u03BB",
                   title = "Uncertainty decomposition by vital rate",
                   fill = "Vital rate" ) +
    ggplot2::scale_fill_manual(values = vr_colors) +
    ggplot2::scale_color_manual( name = NULL,
                                 values = c("Total variance in lambda" = "black") ) +
    ggplot2::guides( fill = ggplot2::guide_legend(order = 1),
                     color = ggplot2::guide_legend(order = 2) ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( axis.text.x = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank() )

  p
}
