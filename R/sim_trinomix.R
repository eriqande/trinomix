# Random generation of of trinomix data sets directly from the Dirichlet dsn
#'
#' This is done with no stick breaking.
#'
#' @param n_obs Number of observations (rows of data matrix to simulate). Defaults to 10
#' @param n_groups Number of categories for each observation (columns of data matrix). Defaults to 10
#' @param ess_fraction The effective sample size fraction, defaults to 1
#' @param tot_n The total sample size to simulate for each observation. This is approximate and the actual
#' simulated sample size will be slightly smaller. Defaults to 100
#' @param p The stock proportions to simulate from, as a vector. Optional, and when not included,
#' random draws from the dirichlet are used
#'
#' @export
#' @importFrom gtools rdirichlet
#' @importFrom stats qbeta rgamma
#'
sim_trinomix <- function(n_obs = 1000,
                         n_groups = 10,
                         ess_fraction = 1,
                         tot_n = 100,
                         p = NULL
) {

  if (is.null(p)) {
    p <- gtools::rdirichlet(1, rep(1, n_groups))
  } else {
    stopifnot(length(p) == n_groups)
  }

  # get some different variable names here
  N <- tot_n
  K <- n_groups
  ess <- N * ess_fraction
  a <- p * ess # these are the alphas
  abd <- sum(a)  # abd = \alpha_\bigdot
  b <- abd - a # the beta parameters for the corresponding marginals

  # simulate all the gammas you need. Set their scale to 1.
  # each COLUMN is a separate set of gammas
  G <- matrix(
    rgamma(n = K * n_obs, shape = a, scale = 1),
    nrow = K,
  )

  # normalize those to Dirichlet r.v.'s
  D <- sweep(G, 2, colSums(G), "/")

  # get the desired zero probabilities
  p0 <- (1 - p)^ess

  # determine the value of each marginal beta at which the
  # cdf for the marginal beta is p0
  qs <- qbeta(p0, a, b)

  # make an indicator matrix. 0 means it is 0, and 1
  # means it is not.
  Im <- matrix(as.integer(D > qs), nrow = K)

  # zero out the gammas that must be zero and normalize the rest
  Gprime <- G * Im
  Gret <- sweep(Gprime, 2, colSums(Gprime), "/")

  # return transposed version for consistency with broken_stick
  return(
    list(
      X_obs = N * t(Gret),
      p = p
    )
  )

}
