# Descripction:
#
# This function calculates the forces that model a single mutation using the LF-ENM 
# ("Linearly Forced - Elastic Network Model").
#
# Usage:
#
# CalculateForce(l, r, kij, fmax)
#
#  Args:
#    - l: index of the site to mutate.
#    - r: 3 x n.sites matrix containing the equilibrium coordinates of the protein.
#    - kij: n.sites x n.sites matrix containing kij between all sites i and j.
#    - fmax: Maximum value for the forces.
#
#  Returns:
#    f: 3 * n.sites vector containing the forces that model the mutation at site l.
#    sum.fij.square: sum of squares of fijs.
#
CalculateForce = function(l, r, kij, fmax) {
  n.sites = ncol(r)
  f = matrix(0, 3, n.sites)
  sum.fij.square = 0
  for (j in seq(n.sites)){
    if(kij[l, j] > 0.9) {
      rij = r[1:3, j] - r[1:3, l]
      eij = rij/sqrt(sum(rij ^ 2))
      fij = fmax * runif(1, -1, 1)
      fij.square = fij ^ 2
      f[, j] = fij * eij
      f[, l] = f[, l] - f[, j]
      sum.fij.square = sum.fij.square + fij.square
    }
  }
  dim(f) = c(3 * n.sites)
  output = list("f" = f, "sum.fij.square" = sum.fij.square)
  output
}




