% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMC.R
\name{dDMC_grid}
\alias{dDMC_grid}
\title{Generate Grid for PDF of Diffusion Model of Conflict Tasks}
\usage{
dDMC_grid(rt_max = 10, phi, x_res = "default", t_res = "default")
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
\item Coherence parameter (\eqn{s}). Sets stimulus coherence. If \eqn{s = 1}, coherent condition;
if \eqn{s = 0}, neutral condition; if \eqn{s = -1}, incoherent condition.
\item Automatic process amplitude (\eqn{A}). Max value of automatic process.
\item Scale parameter (\eqn{\tau}). Contributes to time automatic process. Time to max
\eqn{t_{max} = (\alpha – 1)*\tau}.
    \item Shape parameter (\eqn{\alpha}). Indicates the shape of the automatic process. Must
have value more than 1 (\eqn{\alpha > 1}).
\item Drift rate of the controlled process (\eqn{\mu_c}).
\item Noise scale (\eqn{\sigma}). Model noise scale parameter.
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
For more details on the model see, for example, \code{\link{dDMC}}.
}
\references{
Ulrich, R., Schröter, H., Leuthold, H., & Birngruber, T. (2015). Automatic
and controlled stimulus processing in conflict tasks: Superimposed diffusion
processes and delta functions. \emph{Cognitive psychology, 78}, 148-174.
}
\author{
Raphael Hartmann & Matthew Murrow
}
