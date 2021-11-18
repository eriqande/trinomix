# Random generation of datasets to show that marginals can be satisfied by different distributions
#'
#' Random generation of datasets using the dirichlet broken stick method
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
#' @importFrom stats median pbeta qbeta quantile rbeta runif
#'
#' @examples
#' \donttest{
#' y <- broken_stick(n_obs = 3, n_groups = 5, tot_n = 100)
#'
#' # add custom proportions
#' y <- broken_stick(
#'   n_obs = 3, n_groups = 5, tot_n = 100,
#'   p = c(0.1, 0.2, 0.3, 0.2, 0.2)
#' )
#' }
broken_stick2 <- function(n_obs = 1000,
                         n_groups = 10,
                         ess_fraction = 1,
                         tot_n = 100,
                         p = NULL) {
  if (is.null(p)) {
    p <- gtools::rdirichlet(1, rep(1, n_groups))
  }

  ess <- tot_n * ess_fraction
  # first, determine the presence of zeros using the stick breaking algorithm
  # second, for instances where zeros and tot_n are not observed, rescale the parameters and
  #         use the stick breaking algorithm for the Dirichlet process a second time

  # ECA: so, what we have here is that p is going to define the means of the
  # Dirichlet distribution that these things will be coming from.  (i.e. p is
  # the normalized parameters of the Dirichlet distribution.)
  X_mean_prob <- matrix(p, n_obs, n_groups, byrow = T)

  # ECA: Down here, they are figuring out what they would like the marginal variance for
  # each component of a single replicate of the Dirichlet random vector to be, given
  # the effective sample size.  So, they are thinking about that here as the variance
  # of an actual Dirichlet random vector with parameters given by p * ess.  Where sum(p)
  # is one and ess is what is called phi in Douma and Weedon, box 1.
  X_var_prob <- matrix((p * (1 - p)) / (ess + 1),
                       n_obs, n_groups,
                       byrow = T
  )
  # ECA: I guess that they call these X_mean_prob and X_var_prob to indicate that they
  # are talking about the mean and variance of the actual Dirichlet random vector
  # that sums to 1 (i.e., looks like a probability.)


  # ECA: Here are the means and variances of the actual continuous variables they
  # are going to want in their model (ignoring values of 0 and N).  It is interesting
  # to note that this indicates that, essentially, X = tot_N * X_prob.  In other words,
  # they are casting their vectors as a scaled Dirichlet random vector, rather than
  # as a sum of tot_n iid Dirichlet random vectors. (This is clear because the variance
  # of X is tot_n^2, rather than tot_n of the X_var_prob).  So, I think they have to
  # be careful since their examples seem to be modeling sums of proportions. At least the
  # mixed fishery one.  (In practice, I think that it means their approach will
  # overestimate the phi parameter if the data are sums of Dirichlet random vectors
  # rather than sums of them. i.e. things will be estimated to be less dispersed than
  # they really are).
  X_mean <- X_mean_prob * tot_n
  X_var <- X_var_prob * tot_n^2

  # ECA: Presumably this is a matrix of values that will be compared
  # against a CDF later on. (But actually it really isn't used again).
  rand_unif <- matrix(runif(length(X_mean)), n_obs, n_groups)

  # ECA, start everything off as O?
  X_obs <- rand_unif * 0
  X_indicator <- X_obs


  #### Now, here we can create X_indicator a different way ####
  ## make a function to compute the realized probs of 0's from P0 and PN. ##
  realized_0probs <- function(P0, PN) {
    Bsum <- sum( (1 - P0) * prod(P0) / P0)
    Asum <- Bsum - (1 - P0) * prod(P0) / P0
    Exp0s <- (
      (1 - sum(PN)) *
        (P0 - prod(P0) - Asum) / (1 - prod(P0) - Bsum)
    ) +
      (sum(PN) * (1 - PN / sum(PN)))

    names(Exp0s) <- paste0("P", 1:length(Exp0s))
    Exp0s
  }

  ## Make a function based on the above for root-finding ##
  func_for_roots <- function(x, PN, D) {
    realized_0probs(x, PN) - D
  }

  ## Find the required values of P0 by solving the system of equations ##
  D <- (1 - p) ^ ess
  init <- D
  PN <- p ^ ess

  result <- rootSolve::multiroot(f = func_for_roots, start = init, PN = PN, D = D)

  ## Define a function for creating the Indicator matrix ##
  simulate_indicator_matrix <- function(n_obs, P0, PN) {
    NumNs <- rbinom(n = 1, size = n_obs, prob = sum(p ^ ess))
    NumNotNs <- n_obs - NumNs
    Ncounts <- t(rmultinom(n = NumNs, size = 1, prob = (p ^ ess) / sum(p ^ ess)))

    # and then we also simulate the ones without those Ns.  We will
    # oversimulate, and then just take the ones we need.
    XI <- matrix(
      runif(n = floor(2 * n_obs * n_groups)) > P0,
      ncol = n_groups,
      nrow = 2 * n_obs,
      byrow = TRUE
    )
    all0 <- rowSums(XI) == 0
    all0_but1 <- rowSums(XI) == 1
    XIf <- XI[!(all0 | all0_but1), ]

    # then merge those together:
    SimMat <- rbind(Ncounts, XIf[1:NumNotNs,])

    # permute the rows of the final result
    SimMat[sample(1:nrow(SimMat)), ]
  }


  ## Then, simulate that new Indicator matrix. ##
  X_indicator <- simulate_indicator_matrix(n_obs, P0 = result$root, PN = PN)

  #### Done with creating X_indicator a different way ####

  ### COMMENT OUT THE STUFF WE NO LONGER NEED with a shebang: #! ####
  # Scale out the sample size to move to proportion space.
  # calculate the probability of 0 and probability of 1 first
  # ECA: Since the still have X_mean_prob laying around, they can use that
  # without actually scaling anything.  This merely is giving the probability
  # that each component is identically zero or N (not sure why they call it
  # p_one.  This corresponds to the first two lines of equation 1.
  #! p_zero <- (1 - X_mean_prob)^ess
  #! p_one <- X_mean_prob^ess

  # Unconditional mean calculation
  # ECA: These don't get used subsequently, so I am not going to worry about
  # them.  It is pretty clear what they are, though.  The COND.MEAN is conditional
  # on everything being greater than 1.
  #! UNCOND.MEAN <- p_one[1, ] * tot_n + (1 - p_zero[1, ] - p_one[1, ]) * X_mean[1, ]
  #! COND.MEAN <- X_mean[1, ]

  # Method of moments to calculated alpha from the dirichlet
  # ECA: Here, following the notation of Douma and Weedon they are saying,
  # "we want alpha = mu_c * phi.  We know the variance we want, so we can solve
  # for phi using that, and then just multiply that by mu_c."
  # HOWEVER, the way they have formulated things here X_alpha is always
  # going to be the same as X_mean_prob * ess, so these are some
  # unnecessary machinations. (But might come in handy if you want to
  # define X_alpha a little differently in other applications.)
  X_alpha <- (X_mean_prob) * ((X_mean_prob * (1 - X_mean_prob) / X_var_prob) - 1)
  X_alpha_mod <- X_alpha * 0
  # Calculate the betas for the marginal Beta distribution for potential later use
  X_beta <- (1 - X_mean_prob) * ((X_mean_prob * (1 - X_mean_prob) / X_var_prob) - 1)


  # ECA: Got it. The mu_vals are the beta parts and he q_vals are the same things
  # but rescaled via stick breaking.
  mu_vals <- X_mean * 0 # These will be independent Beta draws, conditioned on being non-zero.
  q_vals <- mu_vals # These will be equivalent to dirichlet draws (mu_vals modified by stick-breaking algorithm)


  # ECA: Standard stick breaking...
  #! for (i in 1:n_obs) {
     BREAK <- "FALSE"
  #!   for (j in 1:(n_groups - 1)) { # Loop over stocks for the Multinomial component
      # ECA: each mu_val is a beta with a_j and sum_{a_j+1}^K.
  #!     mu_vals[i, j] <- rbeta(1, X_alpha[i, j], sum(X_alpha[i, (j + 1):n_groups]))
  #!     if (j == 1) {
  #!       q_vals[i, j] <- mu_vals[i, j]
  #!     } else if (j > 1) {
        # ECA: scale each such beta according to how much stick remains.
  #!       q_vals[i, j] <- prod(1 - mu_vals[i, (1:j - 1)]) * mu_vals[i, j]
  #!     }
  #!   }
    # ECA, the last component gets whatever remains.
  #!   q_vals[i, n_groups] <- 1 - sum(q_vals[i, (1:n_groups - 1)])
    # X_indicator[i,] <- rmultinom(1,tot_n,q_vals[i,])
  #! }

  # Calculate the quantile of each of the q_vals from their respective marginal Betas
  #! X_cdf <- pbeta(q_vals, X_alpha, X_beta)
  #! X_indicator <- X_cdf - p_zero
  #! X_indicator[X_indicator < 0] <- 0  # ECA: this must mean it is meant to be identically 0.
  #! X_indicator[X_indicator > 0] <- 1

  # ECA: At this juncture, it appears that they are getting ready to call things
  # zero if the probability of simulating the q_val (which, remember is a component of
  # a Dirichlet random vector), or something smaller, from a beta distribution (that
  # is in fact the marginal dsn for that component of the Dirichlet random vector),
  # is less than the pzero.  And now,
  # they go through the stick breaking again, but leave out those values that are
  # identically 0.
  q_vals <- q_vals * 0
  mu_vals <- mu_vals * 0

  # ECA: simply a second round of stick breaking, but the Alphas may be different
  # this time, because some of them might be zeroed out completely. Since those alphas
  # have changed, then I can see why they would not want to use the same mu_vals from
  # before.  I had wondered whether they could just be rescaled to sum to one, but if that
  # happens, they would be coming from a different distribution.  Of course, it is not
  # entirely clear what is wanted, because the distribution is only specified in terms
  # of its marginals.
  for (i in 1:n_obs) {
    if (BREAK == "FALSE") {
      # modify the alphas so that they are simply 0 for the components that are identically 0
      X_alpha_mod[i, ] <- X_alpha[i, ] * X_indicator[i, ]
      for (j in 1:(n_groups - 1)) { # Loop over stocks for dirichlet component
        if (j == 1) {
          if (X_alpha_mod[i, j] > 0) {
            mu_vals[i, j] <- rbeta(1, X_alpha_mod[i, j], sum(X_alpha_mod[i, (j + 1):n_groups]))
            q_vals[i, j] <- mu_vals[i, j]
          }
        } else if (j > 1) {
          if (X_alpha_mod[i, j] > 0) {
            mu_vals[i, j] <- rbeta(1, X_alpha_mod[i, j], sum(X_alpha_mod[i, (j + 1):n_groups]))
            q_vals[i, j] <- prod(1 - mu_vals[i, (1:j - 1)]) * mu_vals[i, j]
          }
        }
      }
      if (X_alpha_mod[i, n_groups] > 0) {
        q_vals[i, n_groups] <- 1 - sum(q_vals[i, (1:n_groups - 1)])
      }
      X_obs[i, ] <- q_vals[i, ] * (tot_n)
    }
  }
  return(list(X_obs = X_obs, p = p))
}
