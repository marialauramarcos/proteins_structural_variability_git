# This function calculates ANM kij between 2 sites i and j.
#
#  Args:
#    dij: Distance between equilibrium coordinates of sites i and j.
#    R0: cut-off for ANM.
#
#  Returns:
#    kij: kij ANM for sites i and j.
CalculateKij = function(dij, R0) { 
  if (dij <= R0) {
    kij = 1
  } else {
    kij = 0
  }
  kij
}


