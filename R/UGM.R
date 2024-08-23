


#' Urgency Gating Model
#'
#' Density (PDF), distribution function (CDF), and random sampler for the urgency gating
#'   model (UGM) by Cisek, Puskas, and El-Murr (2009).
#'
#' @param rt vector of response times
#' @param resp vector of responses ("upper" and "lower")
#' @param n number of samples
#' @param phi parameter vector in the following order:
#'   \itemize{
#'     \item non-decision time
#'     \item relative starting point
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item diffusion rate
#'     \item upper threshold (equals negative lower threshold)
#'     \item contamination strength
#'     \item contamination probability for the lower response
#'     \item contamination probability for the upper response
#'   }
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @param dt step size of time. We recommend 0.00001 (1e-5)
#' @return For the density a list of PDF values, log-PDF values, and the sum of the
#'   log-PDFs, for the distribution function a list of of CDF values, log-CDF values,
#'   and the sum of the log-CDFs, and for the random sampler a list of response
#'   times (rt) and response thresholds (resp).
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter
#'   inference for a wide class of binary evidence accumulation models.
#'   \emph{Behavior Research Methods}, 1-21.
#'
#' Cisek, P., Puskas, G. A., & El-Murr, S. (2009). Decisions in changing conditions:
#'   the urgency-gating model. \emph{Journal of Neuroscience, 29}(37), 11560-11571.
#' @examples
#' # Probability density function
#' dUGM(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'      phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0))
#'
#' # Cumulative distribution function
#' pUGM(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'      phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0))
#'
#' # Random sampling
#' rUGM(n = 100, phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0))
#' @author Raphael Hartmann & Matthew Murrow
#' @name UGM
NULL




########### PDF ###########



#' @rdname UGM
#' @useDynLib "ream", .registration=TRUE
#' @export
dUGM <- function(rt,
                 resp,
                 phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                 x_res = "default",
                 t_res = "default") {


  # constants
  modelname <- "UGM"
  Nphi <- 10


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



#' @rdname UGM
#' @useDynLib "ream", .registration=TRUE
#' @export
pUGM <- function(rt,
                 resp,
                 phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                 x_res = "default",
                 t_res = "default") {


  # constants
  modelname <- "UGM"
  Nphi <- 10


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



#' @rdname UGM
#' @useDynLib "ream", .registration=TRUE
#' @export
rUGM <- function(n,
                 phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                 dt = 0.00001) {

  # constants
  modelname <- "UGM"
  Nphi <- 10


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



#' Generate Grid for PDF of the Urgency Gating Model
#'
#' Beschreibung.
#'
#' @param rt_max maximal response time <- max(rt)
#' @param phi parameter vector in the following order:
#'   \itemize{
#'     \item non-decision time
#'     \item relative starting point
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item diffusion rate
#'     \item upper threshold (equals negative lower threshold)
#'     \item contamination strength
#'     \item contamination probability for the lower response
#'     \item contamination probability for the upper response
#'   }
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @return such and such
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter inference for a wide class of binary evidence accumulation models.
#'   Behavior Research Methods.
#'
#' Cisek, P., Puskas, G. A., & El-Murr, S. (2009). Decisions in changing conditions:
#'   the urgency-gating model. \emph{Journal of Neuroscience, 29}(37), 11560-11571.
#'
#' @examples
#' # here come some examples
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "ream", .registration=TRUE
#' @export
dUGM_grid <- function(rt_max = 10.0,
                      phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                      x_res = "default",
                      t_res = "default") {


  # constants
  modelname <- "UGM"
  Nphi <- 10


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
