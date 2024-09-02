


#' Custom Time-Dependent Drift Diffusion Model
#'
#' Density (PDF), distribution function (CDF), and random sampler for a custom time-dependent (CSTM_T) drift diffusion model.
#'
#' @param rt vector of response times
#' @param resp vector of responses ("upper" and "lower")
#' @param n number of samples
#' @param phi parameter vector in your specified order
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @param dt step size of time. We recommend 0.00001 (1e-5)
#' @return For the density a list of PDF values, log-PDF values, and the sum of the
#'   log-PDFs, for the distribution function a list of of CDF values, log-CDF values,
#'   and the sum of the log-CDFs, and for the random sampler a list of response
#'   times (rt) and response thresholds (resp).
#' @examples
#' # Probability density function
#' dCSTM_T(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'       phi = c(0.3, 0.5, 1.0, 1.0, 0.75, 0.0, 0.0, 1.0))
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter
#'   inference for a wide class of binary evidence accumulation models.
#'   \emph{Behavior Research Methods}, 1-21.
#' @author Raphael Hartmann & Matthew Murrow
#' @name CSTM_T
NULL




########### PDF ###########



#' @rdname CSTM_T
#' @useDynLib "ream", .registration=TRUE
#' @export
dCSTM_T <- function(rt,
                  resp,
                  phi,
                  x_res = "default",
                  t_res = "default") {


  # constants
  modelname <- "CSTM_T"
  phi_len <- length(phi)
  phi <- c(phi, rep(0, 100-phi_len))
  Nphi <- 100


  # check
  dist_checks(rt, resp, phi, Nphi, x_res, t_res, modelname)


  # more specific checks


  # setting options
  opt <- dist_options(rt, x_res, t_res)


  # get separated RTs for lower and upper response and get order
  len_rt <- length(rt)
  ind_l <- which(resp=="lower")
  RTL <- rt[ind_l]
  order_l <- order(RTL)
  ind_u <- which(resp=="upper")
  RTU <- rt[ind_u]
  order_u <- order(RTU)


  # prepare arguments for .Call
  dt_scale <- N_deps <- NULL
  REAL <- c(dt_scale = opt[[3]], rt_max = opt[[1]], phi = phi)
  REAL_RTL <- as.double(RTL[order_l])
  REAL_RTU <- as.double(RTU[order_u])
  INTEGER <- c(N_deps = opt[[2]], N_rtl = length(REAL_RTL), N_rtu = length(REAL_RTU), Nphi = length(phi))
  CHAR <- modelname


  # call C++ function
  out <- .Call("PDF",
               as.double(REAL),
               as.integer(INTEGER),
               as.double(REAL_RTL),
               as.double(REAL_RTU),
               as.character(CHAR))


  # transform output
  out$pdf <- numeric(length = len_rt)
  out$pdf[ind_l] <- out$likl[order_l]
  out$pdf[ind_u] <- out$liku[order_u]
  out$log_pdf <- numeric(length = len_rt)
  out$log_pdf[ind_l] <- out$loglikl[order_l]
  out$log_pdf[ind_u] <- out$logliku[order_u]
  out$likl <- out$liku <- out$loglikl <- out$logliku <- NULL


  return(out)

}




########### CDF ###########



#' @rdname CSTM_T
#' @useDynLib "ream", .registration=TRUE
#' @export
pCSTM_T <- function(rt,
                  resp,
                  phi,
                  x_res = "default",
                  t_res = "default") {


  # constants
  modelname <- "CSTM_T"
  phi_len <- length(phi)
  phi <- c(phi, rep(0, 100-phi_len))
  Nphi <- 100

  # check
  dist_checks(rt, resp, phi, Nphi, x_res, t_res, modelname)


  # more specific checks


  # setting options
  opt <- dist_options(rt, x_res, t_res)


  # get separated RTs for lower and upper response and get order
  len_rt <- length(rt)
  ind_l <- which(resp=="lower")
  RTL <- rt[ind_l]
  order_l <- order(RTL)
  ind_u <- which(resp=="upper")
  RTU <- rt[ind_u]
  order_u <- order(RTU)


  # prepare arguments for .Call
  dt_scale <- N_deps <- NULL
  REAL <- c(dt_scale = opt[[3]], rt_max = opt[[1]], phi = phi)
  REAL_RTL <- as.double(RTL[order_l])
  REAL_RTU <- as.double(RTU[order_u])
  INTEGER <- c(N_deps = opt[[2]], N_rtl = length(REAL_RTL), N_rtu = length(REAL_RTU), Nphi = length(phi))
  CHAR <- modelname


  # call C++ function
  out <- .Call("CDF",
               as.double(REAL),
               as.integer(INTEGER),
               as.double(REAL_RTL),
               as.double(REAL_RTU),
               as.character(CHAR))


  # transform output
  out$cdf <- numeric(length = len_rt)
  out$cdf[ind_l] <- out$CDFlow[order_l]
  out$cdf[ind_u] <- out$CDFupp[order_u]
  out$log_cdf <- numeric(length = len_rt)
  out$log_cdf[ind_l] <- out$logCDFlow[order_l]
  out$log_cdf[ind_u] <- out$logCDFupp[order_u]
  out$CDFlow <- out$CDFupp <- out$logCDFlow <- out$logCDFupp <- NULL


  return(out)

}




########### RAND ###########



#' @rdname CSTM_T
#' @useDynLib "ream", .registration=TRUE
#' @export
rCSTM_T <- function(n,
                    phi,
                    dt = 0.00001) {

  # constants
  modelname <- "CSTM_T"
  phi_len <- length(phi)
  phi <- c(phi, rep(0, 100-phi_len))
  Nphi <- 100

  # check arguments
  sim_checks(n, phi, Nphi, dt, modelname)


  # more checks needed for limits etc.


  # prepare arguments for .Call
  REAL <- c(dt = dt, phi = phi)
  INTEGER <- c(N = n, Nphi = length(phi))
  CHAR <- modelname


  # call C++ function
  out <- .Call("SIM",
               as.double(REAL),
               as.integer(INTEGER),
               as.character(CHAR))


  # transform output
  out$resp <- ifelse(out$rt >= 0, "upper", "lower")
  out$rt <- abs(out$rt)


  return(out)

}




########### GRID PDF ###########



#' Generate Grid for PDF of Custom Time-Dependent Drift Diffusion Model
#'
#' Beschreibung.
#'
#' @param rt_max maximal response time <- max(rt)
#' @param phi parameter vector in your order
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @return such and such
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter inference for a wide class of binary evidence accumulation models.
#'   Behavior Research Methods.
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "ream", .registration=TRUE
#' @export
dCSTM_T_grid <- function(rt_max = 10.0,
                         phi,
                         x_res = "default",
                         t_res = "default") {


  # constants
  modelname <- "CSTM_T"
  phi_len <- length(phi)
  phi <- c(phi, rep(0, 100-phi_len))
  Nphi <- 100


  # checking input
  grid_checks(rt_max, phi, Nphi, x_res, t_res, modelname)


  # more specific checks


  # setting options
  opt <- grid_options(x_res, t_res)


  # prepare arguments for r
  dt_scale <- N_deps <- NULL

  CHAR <- modelname

  REAL <- c(dt_scale = dt_scale, rt_max = rt_max, phi = phi)

  INTEGER <- c(N_deps = N_deps, N_phi = length(phi))


  # call C++ function
  out <- .Call("grid_pdf",
               as.double(REAL),
               as.integer(INTEGER),
               as.character(CHAR))



  return(out)

}
