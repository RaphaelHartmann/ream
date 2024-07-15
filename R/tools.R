
dist_checks <- function(rt, resp, phi, Nphi, x_res, t_res, modelname) {

  char_res <- c("default", "high", "higher", "max")

  # checking input
  if (any(rt < 0)) stop("rt must be larger than 0.")
  if (!all(resp %in% c("lower", "upper"))) stop("resp must be either \"upper\" or \"lower\".")
  if (length(phi) != Nphi) stop(paste0("phi must be of length" , Nphi, " for the ", modelname))
  if (!x_res %in% char_res) stop("x_res has not a valid entry")
  if (!t_res %in% char_res) stop("t_res has not a valid entry")
  if (length(resp) != length(rt) & length(resp) != 1) stop("resp must be the same length as rt or of length one")
  if (length(resp) == 1) resp <- rep(resp, length(rt))

}

dist_options <- function(rt, x_res, t_res) {

  char_res <- c("default", "high", "higher", "max")

  # setting options
  x_ind <- which(char_res == x_res)
  t_ind <- which(char_res == t_res)

  out <- vector(mode = "list", length = 3L)
  out[[1]] <- max(rt)
  out[[2]] <- 151 + c(0, 100, 200, 300)[x_ind]
  out[[3]] <- 0.025 * c(1, 0.75, .5, 0.25)[t_ind]

  return(out)

}

sim_checks <- function(n, phi, Nphi, dt, modelname) {

  # checking input
  if (!is.numeric(n) | n %% 1 != 0) stop("n must be a whole number")
  if (length(phi) != Nphi) stop(paste0("phi must be of length" , Nphi, " for the ", modelname))
  if (!is.numeric(dt)) stop("dt must be a numeric value")

}

grid_checks <- function(rt_max, phi, Nphi, x_res, t_res, modelname) {

  char_res <- c("default", "high", "higher", "max")

  # checking input
  if (rt_max < 0) stop("rt must be larger than 0.")
  if (length(phi) != Nphi) stop(paste0("phi must be of length" , Nphi, " for the ", modelname))
  if (!x_res %in% char_res) stop("x_res has not a valid entry")
  if (!t_res %in% char_res) stop("t_res has not a valid entry")

}

grid_options <- function(x_res, t_res) {

  char_res <- c("default", "high", "higher", "max")

  # setting options
  x_ind <- which(char_res == x_res)
  t_ind <- which(char_res == t_res)

  out <- vector(mode = "list", length = 2L)
  out[[1]] <- 151 + c(0, 100, 200, 300)[x_ind]
  out[[2]] <- 0.025 * c(1, 0.75, .5, 0.25)[t_ind]

  return(out)

}
