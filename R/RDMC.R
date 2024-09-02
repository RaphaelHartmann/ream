


#' Revised Diffusion Model of Conflict Tasks
#'
#' A DMC-like model which modifies the shape of the controlled and automatic processes
#'   to ensure consistent stimulus representation across the task. It maintains all SDDM
#'   parameters outside the drift rate which is \eqn{v(x,t) = w_a(t)*d_a + w_c(t)*d_c}, where
#'   \eqn{w_a(t) = A_0*exp(-k*t)} and \eqn{w_c(t) = 1 - w_a(t)}.
#'
#' @param rt vector of response times
#' @param resp vector of responses ("upper" and "lower")
#' @param n number of samples
#' @param phi parameter vector in the following order:
#'   \enumerate{
#'     \item Non-decision time (\eqn{t_{nd}}). Time for non-decision processes such as stimulus
#'       encoding and response execution. Total decision time t is the sum of the decision
#'       and non-decision times.
#'     \item Relative start (\eqn{w}). Sets the start point of accumulation as a ratio of
#'       the two decision thresholds. Related to the absolute start z point via equation
#'       \eqn{z = b_l + w*(b_u - b_l)}.
#'     \item Automatic process amplitude (\eqn{A_0}). Max value of automatic process.
#'     \item Attention shift parameter (\eqn{k}). Encodes congruency and thus differs between
#'       congruent and incongruent trials.
#'     \item Base drift rate of the automatic channel (\eqn{d_a}).
#'     \item Base drift rate of the controlled channel (\eqn{d_c}).
#'     \item Noise scale (\eqn{\sigma}). Model noise scale parameter.
#'     \item Decision thresholds (\eqn{b}). Sets the location of each decision threshold. The
#'       upper threshold \eqn{b_u} is above 0 and the lower threshold \eqn{b_l} is below 0 such that
#'       \eqn{b_u = -b_l = b}. The threshold separation \eqn{a = 2b}.
#'     \item Contamination (\eqn{g}). Sets the strength of the contamination process. Contamination
#'       process is a uniform distribution \eqn{f_c(t)} where \eqn{f_c(t) = 1/(g_u-g_l)}
#'       if \eqn{g_l <= t <= g_u} and \eqn{f_c(t) = 0} if \eqn{t < g_l} or \eqn{t > g_u}. It is
#'       combined with PDF \eqn{f_i(t)} to give the final combined distribution
#'       \eqn{f_{i,c}(t) = g*f_c(t) + (1-g)*f_i(t)}, which is then output by the program.
#'       If \eqn{g = 0}, it just outputs \eqn{f_i(t)}.
#'     \item Lower bound of contamination distribution (\eqn{g_l}). See parameter \eqn{g}.
#'     \item Upper bound of contamination distribution (\eqn{g_u}). See parameter \eqn{g}.
#'   }
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @param dt step size of time. We recommend 0.00001 (1e-5)
#' @return For the density a list of PDF values, log-PDF values, and the sum of the
#'   log-PDFs, for the distribution function a list of of CDF values, log-CDF values,
#'   and the sum of the log-CDFs, and for the random sampler a list of response
#'   times (rt) and response thresholds (resp).
#' @references
#' Lee, P.-S., & Sewell, D. K. (2023). A revised diffusion model for conflict tasks.
#'   \emph{Psychonomic Bulletin &amp; Review, 31}(1), 1–31.
#' @examples
#' # Probability density function
#' dRDMC(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'      phi = c(0.35, 0.5, 7.5, 40.0, 5.0, 5.0, 1.0, 0.5, 0.0, 0.0, 1.0))
#'
#' # Cumulative distribution function
#' pRDMC(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'      phi = c(0.35, 0.5, 7.5, 40.0, 5.0, 5.0, 1.0, 0.5, 0.0, 0.0, 1.0))
#'
#' # Random sampling
#' rRDMC(n = 100, phi = c(0.35, 0.5, 7.5, 40.0, 5.0, 5.0, 1.0, 0.5, 0.0, 0.0, 1.0))
#' @author Raphael Hartmann & Matthew Murrow
#' @name RDMC
NULL




########### PDF ###########



#' @rdname RDMC
#' @useDynLib "ream", .registration=TRUE
#' @export
dRDMC <- function(rt,
                  resp,
                  phi,
                  x_res = "default",
                  t_res = "default") {


  # constants
  modelname <- "RDMC"
  Nphi <- 11


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



#' @rdname RDMC
#' @useDynLib "ream", .registration=TRUE
#' @export
pRDMC <- function(rt,
                 resp,
                 phi,
                 x_res = "default",
                 t_res = "default") {


  # constants
  modelname <- "RDMC"
  Nphi <- 11


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



#' @rdname RDMC
#' @useDynLib "ream", .registration=TRUE
#' @export
rRDMC <- function(n,
                  phi,
                  dt = 0.00001) {

  # constants
  modelname <- "RDMC"
  Nphi <- 11


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



#' Generate Grid for PDF of the Revised Diffusion Model of Conflict Tasks
#'
#' Generate a grid of response-time values and the corresponding PDF values.
#'   For more details on the model see, for example, \code{\link{dRDMC}}.
#'
#' @param rt_max maximal response time <- max(rt)
#' @param phi parameter vector in the following order:
#'   \enumerate{
#'     \item Non-decision time (\eqn{t_{nd}}). Time for non-decision processes such as stimulus
#'       encoding and response execution. Total decision time t is the sum of the decision
#'       and non-decision times.
#'     \item Relative start (\eqn{w}). Sets the start point of accumulation as a ratio of
#'       the two decision thresholds. Related to the absolute start z point via equation
#'       \eqn{z = b_l + w*(b_u - b_l)}.
#'     \item Automatic process amplitude (\eqn{A_0}). Max value of automatic process.
#'     \item Attention shift parameter (\eqn{k}). Encodes congruency and thus differs between
#'       congruent and incongruent trials.
#'     \item Base drift rate of the automatic channel (\eqn{d_a}).
#'     \item Base drift rate of the controlled channel (\eqn{d_c}).
#'     \item Noise scale (\eqn{\sigma}). Model noise scale parameter.
#'     \item Decision thresholds (\eqn{b}). Sets the location of each decision threshold. The
#'       upper threshold \eqn{b_u} is above 0 and the lower threshold \eqn{b_l} is below 0 such that
#'       \eqn{b_u = -b_l = b}. The threshold separation \eqn{a = 2b}.
#'     \item Contamination (\eqn{g}). Sets the strength of the contamination process. Contamination
#'       process is a uniform distribution \eqn{f_c(t)} where \eqn{f_c(t) = 1/(g_u-g_l)}
#'       if \eqn{g_l <= t <= g_u} and \eqn{f_c(t) = 0} if \eqn{t < g_l} or \eqn{t > g_u}. It is
#'       combined with PDF \eqn{f_i(t)} to give the final combined distribution
#'       \eqn{f_{i,c}(t) = g*f_c(t) + (1-g)*f_i(t)}, which is then output by the program.
#'       If \eqn{g = 0}, it just outputs \eqn{f_i(t)}.
#'     \item Lower bound of contamination distribution (\eqn{g_l}). See parameter \eqn{g}.
#'     \item Upper bound of contamination distribution (\eqn{g_u}). See parameter \eqn{g}.
#'   }
#' @param x_res spatial/evidence resolution
#' @param t_res time resolution
#' @return list of RTs and corresponding defective PDFs at lower and upper threshold
#' @references
#' Lee, P.-S., & Sewell, D. K. (2023). A revised diffusion model for conflict tasks.
#'   \emph{Psychonomic Bulletin &amp; Review, 31}(1), 1–31.
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "ream", .registration=TRUE
#' @export
dRDMC_grid <- function(rt_max = 10.0,
                       phi,
                       x_res = "default",
                       t_res = "default") {


  # constants
  modelname <- "RDMC"
  Nphi <- 11


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
