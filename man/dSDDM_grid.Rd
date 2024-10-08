% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDDM.R
\name{dSDDM_grid}
\alias{dSDDM_grid}
\title{Generate Grid for PDF of the Simple Drift Diffusion Model}
\usage{
dSDDM_grid(rt_max = 10, phi, x_res = "default", t_res = "default")
}
\arguments{
\item{rt_max}{maximal response time <- max(rt)}

\item{phi}{parameter vector in the following order:
\enumerate{
\item Non-decision time (\eqn{t_{nd}}). Time for non-decision processes such as stimulus
encoding and response execution. Total decision time t is the sum of the decision
and non-decision times.
\item Relative start (\eqn{w}). Sets the start point of accumulation as a ratio of
the two decision thresholds. Related to the absolute start z point via equation
\eqn{z = b_l + w*(b_u - b_l)}.
\item Stimulus strength (\eqn{\mu}). Strength of the stimulus and used to set the drift
rate. For the SDDM, \eqn{v(x,t) = \mu}.
\item Noise scale (\eqn{\sigma}). Model scaling parameter.
\item Decision thresholds (\eqn{b}). Sets the location of each decision threshold. The
upper threshold \eqn{b_u} is above 0 and the lower threshold \eqn{b_l} is below 0 such that
\eqn{b_u = -b_l = b}. The threshold separation \eqn{a = 2b}.
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
}
\value{
list of RTs and corresponding defective PDFs at lower and upper threshold
}
\description{
Generate a grid of response-time values and the corresponding PDF values.
For more details on the model see, for example, \code{\link{dSDDM}}.
}
\references{
Ratcliff, R. (1978). A theory of memory retrieval. \emph{Psychological Review, 85}(2), 59-108.

Ratcliff, R., & McKoon, G. (2008). The Diffusion Decision Model: Theory and Data
for Two-Choice Decision Tasks. \emph{Neural Computation, 20}(4), 873-922.
}
\author{
Raphael Hartmann & Matthew Murrow
}
