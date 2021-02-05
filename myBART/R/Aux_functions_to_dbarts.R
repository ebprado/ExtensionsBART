# Thanks to Dr. Vincent Dorie for kindly helping out with this!

dbarts_rebuild_tree <- function(tree.flat) {
  if (tree.flat$var[1L] == -1) return(tree.flat$value[1L])

  left <- dbarts_rebuild_tree(tree.flat[-1L,])
  if (!is.list(left)) {
    n.left <- 1L
  } else {
    n.left <- left$n
    left$n <- NULL
  }
  right <- dbarts_rebuild_tree(tree.flat[seq.int(2L + n.left, nrow(tree.flat)),])
  if (!is.list(right)) {
    n.right <- 1L
  } else {
    n.right <- right$n
    right$n <- NULL
  }

  list(var = tree.flat$var[1L], value = tree.flat$value[1L],
       left = left, right = right, n = 1L + n.left + n.right)
}

# Uses a list of lists tree to partition x and returns the
# vector of predictions.
dbarts_get_predictions <- function(tree, x) {
  if (!is.list(tree)) return(rep_len(tree, nrow(x)))

  goes_left <- x[, tree$var] <= tree$value
  predictions <- numeric(nrow(x))

  predictions[ goes_left] <- dbarts_get_predictions(tree$left,  x[ goes_left, , drop = FALSE])
  predictions[!goes_left] <- dbarts_get_predictions(tree$right, x[!goes_left, , drop = FALSE])

  predictions
}

# Set this to the variable you want to only require that trees contain. By
# setting it to -1L, we can double check that the method returns the correct
# result by including all trees.
filter_variable <- -1L
filter_variable <- 4

filter_trees_out <- function(trees.flat) {
  tree_predictions <- by(trees.flat[c("var", "value")], trees.flat$tree,
                         function(tree.flat)
                         {
                           # If tree doesn't contain filter variable, return 0.
                           if (all(tree.flat$var != filter_variable))
                             return(numeric(length(y)))

                           tree <- dbarts_rebuild_tree(tree.flat)

                           dbarts_get_predictions(tree, x)
                         })
  tree_predictions <- matrix(unlist(tree_predictions), nrow = length(y))

  rowSums(tree_predictions)
}

# Uses the flatted trees to get predictions and returns them as a 2-d list,
# n_chains x n_samples.

#' @export
dbarts_marginal_prediction <- function(object, filter_variable){

  trees <- object$fit$getTrees()

  samples_list <- by(trees[c("tree", "var", "value")],
                     trees[c("chain", "sample")],
                     filter_trees_out)

  # Convert list
  samples <- array(unlist(samples_list),
                   c(length(y), dim(samples_list)[1L], dim(samples_list)[2L]))

  # permute to match stored values in bart object
  samples <- aperm(samples, c(2L, 3L, 1L))

  # undo dbarts' internal scaling
  samples <- diff(range(y)) * (samples + 0.5) + min(y)

  # Take the mean of the predicted values over the chains
  samples = apply(apply(samples, c(1,3), mean),2 ,mean)

  return(samples)
}