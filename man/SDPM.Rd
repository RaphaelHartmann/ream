% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDPM.R
\name{SDPM}
\alias{SDPM}
\alias{dSDPM}
\alias{pSDPM}
\alias{rSDPM}
\title{Sequential Dual Process Model}
\usage{
dSDPM(rt, resp, phi, x_res = "default", t_res = "default")

pSDPM(rt, resp, phi, x_res = "default", t_res = "default")

rSDPM(n, phi, dt = 1e-05)
}
\arguments{
\item{rt}{vector of response times}

\item{resp}{vector of responses ("upper" and "lower")}

\item{phi}{parameter vector in the following order:
\enumerate{
\item Non-decision time (\eqn{t_{nd}}). Time for non-decision processes such as stimulus
encoding and response execution. Total decision time t is the sum of the decision
and non-decision times.
\item Relative start (\eqn{w}). Sets the start point of accumulation as a ratio of
the two decision thresholds. Related to the absolute start z point via equation
\eqn{z = b_l + w*(b_u - b_l)}.
\item Relative start of the target selection process (\eqn{w_{ts}}). Sets the start point
of accumulation for the target selection process as a ratio of the two decision
thresholds. Related to the absolute start \eqn{z_{ts}} point via equation
\eqn{z_{ts} = b_{lts} + w_ts*(b_{uts} – b_{lts})}.
    \item Stimulus strength (\eqn{\mu}).
\item Stimulus strength of process 2 (\eqn{\mu_2}).
\item Noise scale (\eqn{\sigma}). Model scaling parameter.
\item Effective noise scale of continuous approximation (\eqn{\sigma_{eff}}). See ream
publication for full description.
\item Decision thresholds (\eqn{b}). Sets the location of each decision threshold. The
upper threshold \eqn{b_u} is above 0 and the lower threshold \eqn{b_l} is below 0 such that
\eqn{b_u = -b_l = b}. The threshold separation \eqn{a = 2b}.
\item Target selection decision thresholds (\eqn{b_{ts}}). Sets the location of each decision
threshold for the target selection process. The upper threshold \eqn{b_{uts}} is above 0
and the lower threshold \eqn{b_{lts}} is below 0 such that \eqn{b_{uts} = -b_{lts} = b_{ts}}. The
threshold separation \eqn{a_{ts} = 2b_{ts}}.
\item Contamination (\eqn{g}). Sets the strength of the contamination process. Contamination
process is a uniform distribution \eqn{f_c(t)} where \eqn{f_c(t) = 1/(g_u-g_l)}
if \eqn{g_l <= t <= g_u} and \eqn{f_c(t) = 0} if \eqn{t < g_l} or \eqn{t > g_u}. It is
combined with PDF \eqn{f_i(t)} to give the final combined distribution
\eqn{f_{i,c}(t) = g*f_c(t) + (1-g)*f_i(t)}, which is then output by the program.
If \eqn{g = 0}, it just outputs \eqn{f_i(t)}.
\item Lower bound of contamination distribution (\eqn{g_l}). See parameter \eqn{g}.
\item Upper bound of contamination distribution (\eqn{g_u}). See parameter \eqn{g}.
}}

\item{x_res}{spatial/evidence resolution}

\item{t_res}{time resolution}

\item{n}{number of samples}

\item{dt}{step size of time. We recommend 0.00001 (1e-5)}
}
\value{
For the density a list of PDF values, log-PDF values, and the sum of the
log-PDFs, for the distribution function a list of of CDF values, log-CDF values,
and the sum of the log-CDFs, and for the random sampler a list of response
times (rt) and response thresholds (resp).
}
\description{
The Sequential Dual Process Model (SDPM) is similar in principle to the DSTP, but instead
of simultaneous accumulators, it contains sequential accumulator s. Its drift rate is given by
\eqn{v(x,t) = w(t)*\mu} where \eqn{w(t)} is 0 if the second process hasn't crossed a
threshold yet and 1 if it has. The noise scale has a similar structure \eqn{D(x,t) = w(t)*\sigma}.
}
\examples{
# Probability density function
dSDPM(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
      phi = c(0.3, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.0, 0.0, 1.0))

# Cumulative distribution function
pSDPM(rt = c(1.2, 0.6, 0.4), resp = c("upper", "lower", "lower"),
      phi = c(0.3, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.0, 0.0, 1.0))

# Random sampling
rSDPM(n = 100, phi = c(0.3, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.0, 0.0, 1.0),
      dt = 0.001)
}
\references{
Hübner, R., Steinhauser, M., & Lehle, C. (2010). A dual-stage two-phase model of
selective attention. \emph{Psychological Review, 117}(3), 759-784.
}
\author{
Raphael Hartmann & Matthew Murrow
}
