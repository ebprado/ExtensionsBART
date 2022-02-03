# -------------------------------------------------------------------------#
# Description: this script contains 2 functions that are used to generate  #
#              the predictions, update variance and compute the tree prior #
#              and the marginalised likelihood                             #
# -------------------------------------------------------------------------#

# 1. simulate_mu: generate the predicted values (mu's)
# 2. update_sigma2: updates the parameters sigma2
# 3. update_z: updates the latent variables z. This is required for MOTR-BART for classification.
# 4. get_tree_prior: returns the tree log prior score
# 5. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree
# 6. get_number_distinct_cov: counts the number of distinct covariates that are used in a tree to create the splitting rules
# Compute the full conditionals -------------------------------------------------

tree_full_conditional = function(tree, R, sigma2, sigma2_mu, m, w) {

  # Function to compute log full conditional distirbution for an individual tree
  # R is a vector of partial residuals

  # Need to calculate log complete conditional, involves a sum over terminal nodes

  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  # Get node sizes for each terminal node
  nj = tree$tree_matrix[which_terminal,'node_size']

  # Get sum of residuals and sum of residuals squared within each terminal node
  sumRsq_j = aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[,2]
  S_j = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate the log posterior
  log_post = 0.5 * ( sum(log( (sigma2 * m^2) / (nj*sigma2_mu*w^2 + sigma2 * m^2))) +
              sum( (sigma2_mu * w^2 * S_j^2) / (sigma2 * (nj*sigma2_mu*w^2 + m^2 * sigma2))))
  return(log_post)
}


# Simulate_par -------------------------------------------------------------

simulate_mu = function(tree, R, sigma2, sigma2_mu, m, w) {

  # Simulate mu values for a given tree

  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  # Get node sizes for each terminal node
  nj = tree$tree_matrix[which_terminal,'node_size']

  # Get sum of residuals in each terminal node
  sumR = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate mu values
  mu = rnorm(length(nj),
             mean = ((sumR*w) / (sigma2 * m)) / ((nj*w^2)/(sigma2 * m^2) + 1/sigma2_mu),
             sd = sqrt(1/((nj*w^2)/(sigma2 * m^2) + 1/sigma2_mu)))

  # Wipe all the old mus out for other nodes
  tree$tree_matrix[,'mu'] = NA

  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal,'mu'] = mu

  return(tree)
}

# Update sigma2 -------------------------------------------------------------

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}

# Update weights ------------------

# tree = curr_trees[[j]]
# w = w[j]
# R = current_partial_residuals
# m = 1

log_sum_exp = function(x) {
  m = max(x);
  return(m + log(sum(exp(x - m))));
}

vec_tree_full_conditional = Vectorize(tree_full_conditional, "w")

update_w <- function(tree, candidate_w, R, ntrees, m, sigma2, sigma2_mu, b_w_prior = 2){
  a_w_prior = b_w_prior/(ntrees-1) # set the hyperparameter a (w ~ Beta(a,b)) so that E(w) = 1/number of trees (m).
  loglik = vec_tree_full_conditional(tree, R, sigma2, sigma2_mu, m, candidate_w) # log of the marginalised likelihood
  prior = (a_w_prior - 1)*candidate_w + (b_w_prior - 1)*(1-candidate_w) # log of the prior on w
  fc_w = loglik + prior # log of the full conditional of w
  # rescaled_fc_w = exp(fc_w)/ sum(exp(fc_w)) # rescale values so that they sum one
  rescaled_fc_w = exp(fc_w - log_sum_exp(fc_w)) # rescale values so that they sum one
  u = runif(1)
  new_w =candidate_w[min(which(cumsum(rescaled_fc_w) > u))] # inverse transform sampling
  return(new_w)
}

# Update the latent variable z ---------------

update_z = function(y, prediction){

  ny0 = sum(y==0)
  ny1 = sum(y==1)
  z = rep(NA, length(y))

  z[y==0] = rtruncnorm(ny0, a = -Inf, b=0,   mean = prediction[y==0], 1)
  z[y==1] = rtruncnorm(ny1, a = 0   , b=Inf, mean = prediction[y==1], 1)

  return(z)
}

# Get tree priors ---------------------------------------------------------
get_tree_prior = function(tree, alpha, beta) {

  # Need to work out the depth of the tree
  # First find the level of each node, then the depth is the maximum of the level
  level = rep(NA, nrow(tree$tree_matrix))
  level[1] = 0 # First row always level 0

  # Escpae quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0
  }

  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = as.numeric(tree$tree_matrix[i,'parent'])
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }

  # Only compute for the internal nodes
  internal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 0)
  log_prior = 0
  for(i in 1:length(internal_nodes)) {
    log_prior = log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]])
  }
  # Now add on terminal nodes
  terminal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 1)
  for(i in 1:length(terminal_nodes)) {
    log_prior = log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  }

  return(log_prior)

}