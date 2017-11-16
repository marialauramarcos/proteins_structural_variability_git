# This function calculates ANM kij between 2 sites i and j.
#
#  Args:
#    dij: Distance between equilibrium coordinates of sites i and j.
#    R0: cut-off for ANM.
#
#  Returns:
#    kij: kij ANM for sites i and j.
CalculateKij = function(dij, R0) { 
  if (dij > tolerance) {
    kij = 1/(dij ^ 2)
  }
  kij
}


