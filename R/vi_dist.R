#' Variation of information distance
#'
#' (From 'mcclust' R package) Computes the 'variation of information'
#' distance of Meila (2007) between two clusterings/partitions of the same objects.
#'
#' @param cl1 Vector of cluster membership indices
#' @param cl2 Vector of cluster membership indices
#' @param parts logical; should the two conditional entropies also be returned?
#' @param base base of logarithm used for computation of entropy and mutual information
#'
#' @return The VI distance. If parts=TRUE the two conditional entropies are appended.
#' @export
#'
#' @details The variation of information distance is the sum of the two conditional
#' entropies of one clustering given the other. For details see Meila (2007).
#' @references
#' Meila, M. (2007) Comparing Clusterings - an Information Based Distance. Journal of Multivariate Analysis, 98, 873 – 895.
vi_dist = function (cl1, cl2, parts = FALSE, base = 2)
{
  if (length(cl1) != length(cl2))
    stop("cl1 and cl2 must have same length")
  ent <- function(cl) {
    n <- length(cl)
    p <- table(cl)/n
    -sum(p * log(p, base = base))
  }
  mi <- function(cl1, cl2) {
    p12 <- table(cl1, cl2)/length(cl1)
    p1p2 <- outer(table(cl1)/length(cl1), table(cl2)/length(cl2))
    sum(p12[p12 > 0] * log(p12[p12 > 0]/p1p2[p12 > 0], base = base))
  }
  if (!parts)
    return(ent(cl1) + ent(cl2) - 2 * mi(cl1, cl2))
  ent1 <- ent(cl1)
  ent2 <- ent(cl2)
  mi12 <- mi(cl1, cl2)
  c(vi = ent1 + ent2 - 2 * mi12, `H(1|2)` = ent1 - mi12, `H(2|1)` = ent2 -
      mi12)
}
