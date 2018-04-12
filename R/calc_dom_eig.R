#' Estimate dominant eigenvalue and vector
#'
#' @param A IPM kernel (IPM$K)
#' @param tol Tolerence level
#' @return list(lambda,w) where w is the eigen vector
#' @seealso \code{\link{domEig_speed}} a cut-down version
#' @export
calc_dom_eig <- function(MIPM,kernel = "K",tol=1e-8) {
  A <- MIPM[["megakernels"]][kernel][[1]]

  qmax <- 10 * tol
  lam <- 1
  x <- rep(1,nrow(A))

  while(qmax>tol) {
    x1 <- A %*% x
    qmax <- sum(abs(x1-lam*x))
    lam <- sum(x1)
    x <- x1/lam
  }

  # Having found w (within tol), get lambda
  x1 <- A %*% x
  lam <- sum(x1)
  x <- x1/lam
  return(list(lambda=lam,w=x/sum(x)))
}
