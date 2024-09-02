


#' Sequential Dual Process Model
#'
#' The Sequential Dual Process Model (SDPM) is similar in principle to the DSTP, but instead
#'   of simultaneous accumulators, it contains sequential accumulator s. Its drift rate is given by
#'   \eqn{v(x,t) = w(t)*\mu} where \eqn{w(t)} is 0 if the second process hasn't crossed a
#'   threshold yet and 1 if it has. The noise scale has a similar structure \eqn{D(x,t) = w(t)*\sigma}.
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
#'     \item Relative start of the target selection process (\eqn{w_{ts}}). Sets the start point
#'       of accumulation for the target selection process as a ratio of the two decision
#'       thresholds. Related to the absolute start \eqn{z_{ts}} point via equation
#'       \eqn{z_{ts} = b_{lts} + w_ts*(b_{uts} – b_{lts})}.
#'     \item Stimulus strength (\eqn{\mu}).
#'     \item Stimulus strength of process 2 (\eqn{\mu_2}).
#'     \item Noise scale (\eqn{\sigma}). Model scaling parameter.
#'     \item Effective noise scale of continuous approximation (\eqn{\sigma_{eff}}). See ream
#'       publication for full description.
#'     \item Decision thresholds (\eqn{b}). Sets the location of each decision threshold. The
#'       upper threshold \eqn{b_u} is above 0 and the lower threshold \eqn{b_l} is below 0 such that
#'       \eqn{b_u = -b_l = b}. The threshold separation \eqn{a = 2b}.
#'     \item Target selection decision thresholds (\eqn{b_{ts}}). Sets the location of each decision
#'       threshold for the target selection process. The upper threshold \eqn{b_{uts}} is above 0
#'       and the lower threshold \eqn{b_{lts}} is below 0 such that \eqn{b_{uts} = -b_{lts} = b_{ts}}. The
#'       threshold separation \eqn{a_{ts} = 2b_{ts}}.
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
#' Hübner, R., Steinhauser, M., & Lehle, C. (2010). A dual-stage two-phase model of
#'   selective attention. \emph{Psychological Review, 117}(3), 759-784.
#' @examples
#' # Probability density function
#' dSDPM(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'       phi = c(0.3, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.0, 0.0, 1.0))
#'
#' # Cumulative distribution function
#' pSDPM(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'       phi = c(0.3, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.0, 0.0, 1.0))
#'
#' # Random sampling
#' rSDPM(n = 100, phi = c(0.3, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.0, 0.0, 1.0),
#'       dt = 0.001)
#' @author Raphael Hartmann & Matthew Murrow
#' @name SDPM
NULL




########### PDF ###########



#' @rdname SDPM
#' @useDynLib "ream", .registration=TRUE
#' @export
dSDPM <- function(rt,
                 resp,
                 phi,
                 x_res = "default",
                 t_res = "default") {


  # constants
  modelname <- "SDPM"
  Nphi <- 12


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



#' @rdname SDPM
#' @useDynLib "ream", .registration=TRUE
#' @export
pSDPM <- function(rt,
                 resp,
                 phi,
                 x_res = "default",
                 t_res = "default") {


  # constants
  modelname <- "SDPM"
  Nphi <- 12


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



#' @rdname SDPM
#' @useDynLib "ream", .registration=TRUE
#' @export
rSDPM <- function(n,
                  phi,
                  dt = 0.00001) {

  # constants
  modelname <- "SDPM"
  Nphi <- 12


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



#' Generate Grid for PDF of the Sequential Dual Process Model
#'
#' Generate a grid of response-time values and the corresponding PDF values.
#'   For more details on the model see, for example, \code{\link{dSDPM}}.
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
#'     \item Relative start of the target selection process (\eqn{w_{ts}}). Sets the start point
#'       of accumulation for the target selection process as a ratio of the two decision
#'       thresholds. Related to the absolute start \eqn{z_{ts}} point via equation
#'       \eqn{z_{ts} = b_{lts} + w_ts*(b_{uts} – b_{lts})}.
#'     \item Stimulus strength (\eqn{\mu}).
#'     \item Stimulus strength of process 2 (\eqn{\mu_2}).
#'     \item Noise scale (\eqn{\sigma}). Model scaling parameter.
#'     \item Effective noise scale of continuous approximation (\eqn{\sigma_{eff}}). See ream
#'       publication for full description.
#'     \item Decision thresholds (\eqn{b}). Sets the location of each decision threshold. The
#'       upper threshold \eqn{b_u} is above 0 and the lower threshold \eqn{b_l} is below 0 such that
#'       \eqn{b_u = -b_l = b}. The threshold separation \eqn{a = 2b}.
#'     \item Target selection decision thresholds (\eqn{b_{ts}}). Sets the location of each decision
#'       threshold for the target selection process. The upper threshold \eqn{b_{uts}} is above 0
#'       and the lower threshold \eqn{b_{lts}} is below 0 such that \eqn{b_{uts} = -b_{lts} = b_{ts}}. The
#'       threshold separation \eqn{a_{ts} = 2b_{ts}}.
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
#' Hübner, R., Steinhauser, M., & Lehle, C. (2010). A dual-stage two-phase model of
#'   selective attention. \emph{Psychological Review, 117}(3), 759-784.
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "ream", .registration=TRUE
#' @export
dSDPM_grid <- function(rt_max = 10.0,
                       phi,
                       x_res = "default",
                       t_res = "default") {


  # constants
  modelname <- "SDPM"
  Nphi <- 12


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
