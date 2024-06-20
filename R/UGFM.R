



########### PDF ###########



#' PDF of the Urgency Gating Flip Model
#'
#' Calculation the (Log-)PDF of the urgency gating flip model (UGFM)
#'
#' @param rt vector of response times
#' @param resp vector of responses ("upper" and "lower")
#' @param phi parameter vector in the following order:
#'   \itemize{
#'     \item non-decision time
#'     \item relative starting point
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item diffusion rate
#'     \item upper threshold (equals negative lower threshold)
#'     \item contamination strength
#'     \item contamination probability for the lower response
#'     \item contamination probability for the upper response
#'   }
#' @param rt_max maximal response time <- max(rt)
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @return a list of PDF values, log-PDF values, and the sum of the log-PDFs
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter
#'   inference for a wide class of binary evidence accumulation models.
#'   \emph{Behavior Research Methods}, 1-21.
#'
#' @examples
#' # here come some examples
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "rbeam", .registration=TRUE
#' @export
dUGFM <- function(rt,
                 resp,
                 phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                 x_res = "A",
                 t_res = "A") {


  # constants
  char_res <- c("A", "B", "C", "D")

  # checking input
  if (any(rt < 0)) stop("rt must be larger than 0.")
  if (!all(resp %in% c("lower", "upper"))) stop("resp must be either \"upper\" or \"lower\".")
  if (length(phi) != 11) stop("phi must be of length 11 for the UGFM")
  if (!x_res %in% char_res) stop("x_res has not a valid entry")
  if (!t_res %in% char_res) stop("t_res has not a valid entry")
  if (length(resp) != length(rt) & length(resp) != 1) stop("resp must be the same length as rt or of length one")
  if (length(resp) == 1) resp <- rep(resp, length(rt))

  # more checks needed for limits etc.


  # setting options
  x_ind <- which(char_res == x_res)
  t_ind <- which(char_res == t_res)

  N_deps <- 151 + c(0, 100, 200, 300)[x_ind]
  dt_scale <- 0.025 * c(1, 0.75, .5, 0.25)[t_ind]

  rt_max <- max(rt)

  # get separated RTs for lower and upper response and get order
  len_rt <- length(rt)
  ind_l <- which(resp=="lower")
  RTL <- rt[ind_l]
  order_l <- order(RTL)
  ind_u <- which(resp=="upper")
  RTU <- rt[ind_u]
  order_u <- order(RTU)


  # prepare arguments for .Call
  REAL <- c(dt_scale = dt_scale, rt_max = rt_max, phi = phi)
  REAL_RTL <- as.double(RTL[order_l])
  REAL_RTU <- as.double(RTU[order_u])
  INTEGER <- c(N_deps = N_deps, N_rtl = length(REAL_RTL), N_rtu = length(REAL_RTU), Nphi = length(phi))
  CHAR <- "UGFM"


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



#' CDF of the Urgency Gating Flip Model
#'
#' Calculation the (Log-)CDF of the urgency gating flip model (UGFM)
#'
#' @param rt vector of response times
#' @param resp vector of responses ("upper" and "lower")
#' @param phi parameter vector in the following order:
#'   \itemize{
#'     \item non-decision time
#'     \item relative starting point
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item diffusion rate
#'     \item upper threshold (equals negative lower threshold)
#'     \item contamination strength
#'     \item contamination probability for the lower response
#'     \item contamination probability for the upper response
#'   }
#' @param rt_max maximal response time <- max(rt)
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @return a list of PDF values, log-PDF values, and the sum of the log-PDFs
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter
#'   inference for a wide class of binary evidence accumulation models.
#'   \emph{Behavior Research Methods}, 1-21.
#'
#' @examples
#' # here come some examples
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "rbeam", .registration=TRUE
#' @export
pUGFM <- function(rt,
                 resp,
                 phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                 x_res = "A",
                 t_res = "A") {


  # constants
  char_res <- c("A", "B", "C", "D")

  # checking input
  if (any(rt < 0)) stop("rt must be larger than 0.")
  if (!all(resp %in% c("lower", "upper"))) stop("resp must be either \"upper\" or \"lower\".")
  if (length(phi) != 11) stop("phi must be of length 11 for the UGFM")
  if (!x_res %in% char_res) stop("x_res has not a valid entry")
  if (!t_res %in% char_res) stop("t_res has not a valid entry")
  if (length(resp) != length(rt) & length(resp) != 1) stop("resp must be the same length as rt or of length one")
  if (length(resp) == 1) resp <- rep(resp, length(rt))

  # more checks needed for limits etc.


  # setting options
  x_ind <- which(char_res == x_res)
  t_ind <- which(char_res == t_res)

  N_deps <- 151 + c(0, 100, 200, 300)[x_ind]
  dt_scale <- 0.025 * c(1, 0.75, .5, 0.25)[t_ind]

  rt_max <- max(rt)

  # get separated RTs for lower and upper response and get order
  len_rt <- length(rt)
  ind_l <- which(resp=="lower")
  RTL <- rt[ind_l]
  order_l <- order(RTL)
  ind_u <- which(resp=="upper")
  RTU <- rt[ind_u]
  order_u <- order(RTU)


  # prepare arguments for .Call
  REAL <- c(dt_scale = dt_scale, rt_max = rt_max, phi = phi)
  REAL_RTL <- as.double(RTL[order_l])
  REAL_RTU <- as.double(RTU[order_u])
  INTEGER <- c(N_deps = N_deps, N_rtl = length(REAL_RTL), N_rtu = length(REAL_RTU), Nphi = length(phi))
  CHAR <- "UGFM"

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



#' Sampling From the Urgency Gating Flip Model
#'
#' Sampling from the urgency gating flip model (UGFM)
#'
#' @param n number of samples
#' @param phi parameter vector in the following order:
#'   \itemize{
#'     \item non-decision time
#'     \item relative starting point
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item ?
#'     \item diffusion rate
#'     \item upper threshold (equals negative lower threshold)
#'     \item contamination strength
#'     \item contamination probability for the lower response
#'     \item contamination probability for the upper response
#'   }
#' @param dt step size of time. We recommend 0.00001 for the UGFM
#' @return list of response times (rt) and response threholds (resp)
#' @references
#' Murrow, M., & Holmes, W. R. (2023). PyBEAM: A Bayesian approach to parameter
#'   inference for a wide class of binary evidence accumulation models.
#'   \emph{Behavior Research Methods}, 1-21.
#'
#' @examples
#' rdmc(n = 100, phi = c(0.3, 0.5, -1.0, 0.2, 0.05, 2.5, 3.0, 1.0, 0.5, 0.0, 0.0, 1.0))
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "rbeam", .registration=TRUE
#' @export
rUGFM <- function(n,
                 phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                 dt = 0.00001) {

  # check arguments
  if (!is.numeric(n) | n %% 1 != 0) stop("n must be a whole number")
  if (length(phi) != 11) stop("phi must be of length 11 for the RDMC")
  if (!is.numeric(dt)) stop("dt must be a numeric value")

  # more checks needed for limits etc.


  # prepare arguments for .Call
  REAL <- c(dt = dt, phi = phi)
  INTEGER <- c(N = n, Nphi = length(phi))
  CHAR <- "UGFM"


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



#' Generate Grid for PDF of the Urgency Gating Flip Model
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
#' @examples
#' # here come some examples
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "rbeam", .registration=TRUE
#' @export
dUGFM_grid <- function(rt_max = 10.0,
                      phi = c(0.3, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0, 1.5, 0.0, 0.0, 1.0),
                      x_res = "A",
                      t_res = "A") {


  # constants
  char_res <- c("A", "B", "C", "D")

  # checking input
  if (any(rt < 0)) stop("rt must be larger than 0.")
  if (!all(resp %in% c("lower", "upper"))) stop("resp must be either \"upper\" or \"lower\".")
  if (length(phi) != 11) stop("phi must be of length 11 for the UGFM")
  if (!x_res %in% char_res) stop("x_res has not a valid entry")
  if (!t_res %in% char_res) stop("t_res has not a valid entry")

  # more checks needed for limits etc.


  # setting options
  x_ind <- which(char_res == x_res)
  t_ind <- which(char_res == t_res)

  N_deps <- 151 + c(0, 100, 200, 300)[x_ind]
  dt_scale <- 0.025 * c(1, 0.75, .5, 0.25)[t_ind]

  # prepare arguments for r
  CHAR <- "UGFM"

  REAL <- c(dt_scale = dt_scale, rt_max = rt_max, phi = phi)

  INTEGER <- c(N_deps = N_deps, N_phi = length(phi))


  # call C++ function
  out <- .Call("grid_pdf",
               as.double(REAL),
               as.integer(INTEGER),
               as.character(CHAR))



  return(out)

}


