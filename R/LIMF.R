


#' Leaky Integration Model With Flip
#'
#' LIM with time varying drift rate. Specifically, the stimulus strength changes from
#'   \eqn{\mu_1} to \eqn{\mu_2} at time \eqn{t_0}. Identified by (Evans et al., 2020; Trueblood et al., 2021)
#'   as a way to improve recovery of the leakage rate. Drift rate becomes
#'   \eqn{v(x,t) = \mu_1 - L*x} if \eqn{t < t_0} and \eqn{v(x,t) = \mu_2 - L*x} if \eqn{t >= t_0.}
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
#'     \item Stimulus strength 1 (\eqn{\mu_1}). Strength of the stimulus prior to \eqn{t_0}.
#'     \item Stimulus strength 2 (\eqn{\mu_2}). Strength of the stimulus after \eqn{t_0}.
#'     \item Log10-leakage (\eqn{log_{10}(L)}). Rate of leaky integration.
#'     \item Flip-time (\eqn{t_0}). Time when stimulus strength changes.
#'     \item Noise scale (\eqn{\sigma}). Model scaling parameter.
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
#' Evans, N. J., Trueblood, J. S., & Holmes, W. R. (2019). A parameter recovery assessment of
#'   time-variant models of decision-making. \emph{Behavior Research Methods, 52}(1), 193-206.
#'
#' Trueblood, J. S., Heathcote, A., Evans, N. J., & Holmes, W. R. (2021). Urgency, leakage,
#'   and the relative nature of information processing in decision-making.
#'   \emph{Psychological Review, 128}(1), 160-186.
#' @examples
#' # Probability density function
#' dLIMF(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'      phi = c(0.3, 0.5, 1.0, 0.9, 0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0))
#'
#' # Cumulative distribution function
#' pLIMF(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
#'      phi = c(0.3, 0.5, 1.0, 0.9, 0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0))
#'
#' # Random sampling
#' rLIMF(n = 100, phi = c(0.3, 0.5, 1.0, 0.9, 0.5, 0.5, 1.0, 0.5, 0.0, 0.0, 1.0))
#' @author Raphael Hartmann & Matthew Murrow
#' @name LIMF
NULL




########### PDF ###########



#' @rdname LIMF
#' @useDynLib "ream", .registration=TRUE
#' @export
dLIMF <- function(rt,
                  resp,
                  phi,
                  x_res = "default",
                  t_res = "default") {


  # constants
  modelname <- "LIMF"
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



#' @rdname LIMF
#' @useDynLib "ream", .registration=TRUE
#' @export
pLIMF <- function(rt,
                  resp,
                  phi,
                  x_res = "default",
                  t_res = "default") {


  # constants
  modelname <- "LIMF"
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



#' @rdname LIMF
#' @useDynLib "ream", .registration=TRUE
#' @export
rLIMF <- function(n,
                  phi,
                  dt = 0.00001) {

  # constants
  modelname <- "LIMF"
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



#' Generate Grid for PDF of the Leaky Integration Model With Flip
#'
#' Generate a grid of response-time values and the corresponding PDF values.
#'   For more details on the model see, for example, \code{\link{dLIMF}}.
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
#'     \item Stimulus strength 1 (\eqn{\mu_1}). Strength of the stimulus prior to \eqn{t_0}.
#'     \item Stimulus strength 2 (\eqn{\mu_2}). Strength of the stimulus after \eqn{t_0}.
#'     \item Log10-leakage (\eqn{log_{10}(L)}). Rate of leaky integration.
#'     \item Flip-time (\eqn{t_0}). Time when stimulus strength changes.
#'     \item Noise scale (\eqn{\sigma}). Model scaling parameter.
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
#' Evans, N. J., Trueblood, J. S., & Holmes, W. R. (2019). A parameter recovery assessment of
#'   time-variant models of decision-making. \emph{Behavior Research Methods, 52}(1), 193-206.
#'
#' Trueblood, J. S., Heathcote, A., Evans, N. J., & Holmes, W. R. (2021). Urgency, leakage,
#'   and the relative nature of information processing in decision-making.
#'   \emph{Psychological Review, 128}(1), 160-186.
#' @author Raphael Hartmann & Matthew Murrow
#' @useDynLib "ream", .registration=TRUE
#' @export
dLIMF_grid <- function(rt_max = 10.0,
                       phi,
                       x_res = "default",
                       t_res = "default") {


  # constants
  modelname <- "LIMF"
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
