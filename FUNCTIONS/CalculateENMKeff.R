# This function calculates K effective of r, its eigenvalues and eigenvectors.
#
#  Args:
#    r: 3 x n.sites matrix containing the equilibrium coordinates of each site.
#    aligned.index: aligned sites.
#    not.aligned.index: not aligned sites.
#    R0: cut-off for ANM.
#    tolerance: 0 tolerance.
#
#  Requires:
#    CalculateENMK
#
#  Returns:
#    KEFF: 3 * n.aligned x 3 * n.aligned KEFF matrix.
#    va: ordered eigenvalues > tolerance of KEFF.
#    ve: ordered eigenvectors > tolerance of KEFF.

CalculateENMKeff <- function(r, 
                             aligned.index, 
                             not.aligned.index, 
                             R0, 
                             tolerance, 
                             K.analysis) {
  
  n.aligned = length(aligned.index)
  n.not.aligned = length(not.aligned.index)
  n.sites = ncol(r)
  
  # Bind aligned and not aligned index
  aligned.not.aligned.index = c(aligned.index, not.aligned.index)
  
  # Calculate K
  ENMK = CalculateENMK(r, CalculateKij, R0, tolerance)
  K = ENMK$K
  
  # K.analysis = "K"
  if (K.analysis == "K") {
    KEFF = K
    vaKEFF = ENMK$va
    veKEFF = ENMK$ve
    n.modes = ENMK$n.modes
  }
  
  # K.analysis = "Keff"
  # Order K in order to get KPP, KQQ, KPQ and KQP and KEFF
  if (K.analysis == "Keff") {
    sortm = matrix(0, ncol = 3 * n.sites, nrow = 1)
    count = 1
    for (i in (1:n.sites)) {
      for (j in (1:n.aligned)) {
        if (aligned.not.aligned.index[j] == i) {
          sortm[1, (3 * count - 2)] = 3 * aligned.not.aligned.index[j] - 2
          sortm[1, (3 * count - 1)] = 3 * aligned.not.aligned.index[j] - 1
          sortm[1, (3 * count)] = 3 * aligned.not.aligned.index[j]
          count = count + 1
        }
      }
    }
    if (n.not.aligned > 0) {
      for (i in (1:n.sites)) {
        for (j in ((n.aligned + 1):n.sites)) {
          if (aligned.not.aligned.index[j] == i) {
            sortm[1, (3 * count - 2)] = 3 * aligned.not.aligned.index[j] - 2
            sortm[1, (3 * count - 1)] = 3 * aligned.not.aligned.index[j] - 1
            sortm[1, (3 * count)] = 3 * aligned.not.aligned.index[j]
            count = count + 1
          }
        }
      }
    }
    Kord = K[sortm, sortm]
    KPP = Kord[(1:(3 * n.aligned)), (1:(3 * n.aligned))]
  
    if(n.not.aligned > 0) {
      KQQ = Kord[((3 * n.aligned + 1):(3 * n.sites)), ((3 * n.aligned + 1):(3 * n.sites))]
      KPQ = Kord[(1:(3 * n.aligned)), ((3 * n.aligned + 1):(3 * n.sites))]
      KQP = t(KPQ)
    
      # Calculate KQQ^-1.
      eigQQ = eigen(KQQ, symmetric = T)
      veQQ = eigQQ$vectors
      vaQQ = eigQQ$values
      modes = vaQQ > tolerance #there aren´t negative evalues#
      veQQ = veQQ[, modes]
      vaQQ = vaQQ[modes]
      covQQ = veQQ %*% ((1 / vaQQ) * t(veQQ))
    
      # Calculate KEFF.
      KEFF <- KPP - (KPQ %*% covQQ %*% KQP)
    }
    if(n.not.aligned == 0) {
      KEFF = KPP
    }
  
    # Calculate and order eigenvalues and eigenvectors of KEFF.
    eigKEFF = eigen(KEFF, symmetric = T)
    veKEFF = eigKEFF$vectors
    vaKEFF = eigKEFF$values
    ord = sort.list(Mod(vaKEFF), decreasing = FALSE)
    vaKEFF = vaKEFF[ord]
    veKEFF = veKEFF[, ord]
    modes = vaKEFF > tolerance  # there aren´t negative evalues.
    veKEFF = veKEFF[, modes]
    vaKEFF = vaKEFF[modes]
    n.modes = length(vaKEFF)
    covKEFF = veKEFF %*% ((1 / vaKEFF) * t(veKEFF))
    
  }
  
  # Create a list for the output 
  output = list("KEFF" = KEFF, "ve" = veKEFF, "va" = vaKEFF, "n.modes" = n.modes, "covKEFF" = covKEFF)
  output
}


