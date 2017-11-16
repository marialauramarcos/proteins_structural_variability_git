# This function calculates measures of variability between two structures. 
#
#  Args:
#    r.p.1: vector of coordinates of p.1.
#    r.p.2: vector of coordinates of p.2.
#    aligned.p.1.index: aligned sites of p.1.
#    aligned.p.2.index: aligned sites of p.2.
#    not.aligned.p.1.index: not aligned sites of p.1.
#    not.aligned.p.2.index: not aligned sites of p.2.
#    R0: cut-off for the ANM.
#    rotate: it can be "TRUE" or "FALSE". If it is "TRUE", r.p.2 is rotated in order to minimize RMSD with r.p.1.
#    K.analysis: It can be "K" or "Keff". For "K" or "Keff", the analysis is based on normal modes of "K" or "Keff"
#    respectibly. 
#    tolerance: 0 tolerance.
# 
#  Required functions:
#    CalculateENMKeff
#    CalculateENMK
#    CalculateKij
#    WindowsRMSD
#    LocalEnviroment
#    CalculateDynamicalVariability
#
#  Returns:
#    Pn
#    dr.squarei
#    va

CalculateVariability <- function(r.p.1, 
                                 r.p.2, 
                                 n.aa.p.1,
                                 n.aa.p.2,
                                 aligned.p.1.index, 
                                 aligned.p.2.index, 
                                 not.aligned.p.1.index,
                                 not.aligned.p.2.index,
                                 R0,
                                 rotate,
                                 K.analysis,
                                 tolerance) {
  
  # calculate 3N indexes
  aligned.p.1.index.3N = sort(c(aligned.p.1.index * 3, aligned.p.1.index * 3 - 2, aligned.p.1.index * 3 - 1))
  aligned.p.2.index.3N = sort(c(aligned.p.2.index * 3, aligned.p.2.index * 3 - 2, aligned.p.2.index * 3 - 1))
  
  # rotate r.p.2 in order to minimize RMSD with r.p.1
  if (rotate == TRUE) {
    r.p.2 = fit.xyz(fixed = r.p.1,
                   mobile = r.p.2,
               fixed.inds = aligned.p.1.index.3N, 
              mobile.inds = aligned.p.2.index.3N) 
  }

  # calculate dr and dr.squarei
  dr = r.p.2[aligned.p.2.index.3N] - r.p.1[aligned.p.1.index.3N]
  dr.squarei = colSums(matrix(dr, nrow = 3) ^ 2) 
  
  # transform r.p.2 using an aa seq based windows and calculate dr.squarei.windows.rot
  dr.squarei.windows.rot = WindowsRMSD(10,
                                       r.p.1,
                                       r.p.2,
                                       aligned.p.1.index.3N,
                                       aligned.p.2.index.3N)
  
  # transform r.p.2 using contacts based windows and calculate a score (dr.squarei.windows.contacts.rot) in each window
  # the score of each window is mean(di^2), being i the sites inside the window
  
  ## calculate ENMK of aligned sites of p.1
  ENMK.p.1 = CalculateENMK(matrix(r.p.1, nrow = 3)[, aligned.p.1.index],
                           CalculateKij,
                           R0,
                           tolerance)
  
  ## trasform r.p.2 and get the score of the windows of each site 
  dr.squarei.windows.contacts.rot = WindowsRMSDcontacts(ENMK.p.1$kij,
                                                        r.p.1,
                                                        r.p.2,
                                                        aligned.p.1.index.3N,
                                                        aligned.p.2.index.3N)
  
  # similarity of local enviroment - local score without transformating coordinates
  ## set R0 for local score
  R0.local.score = R0
  
  ## calculate ENMK of aligned sites of p.1
  ENMK.p.1 = CalculateENMK(matrix(r.p.1, nrow = 3)[, aligned.p.1.index],
                           CalculateKij,
                           R0.local.score,
                           tolerance)
  
  ## get kij of aligned sites of p.1
  ENMkij.p.1 = ENMK.p.1$kij

  ## calculate the distance matrix of aligned sites
  dist.p.1 = as.matrix(dist(t(matrix(r.p.1, nrow = 3)[, aligned.p.1.index]), diag = F, upper = T))
  dist.p.2 = as.matrix(dist(t(matrix(r.p.2, nrow = 3)[, aligned.p.2.index]), diag = F, upper = T))
  
  ## create a vcetor to save the local score of each site
  local.score = c()
  
  ## start a loop for each aligned site
  for (i in (1:length(aligned.p.1.index))) {
    
    ## get distances of contacts
    dist.p.1.i = ENMkij.p.1[i, ] * dist.p.1[i, ]
    dist.p.2.i = ENMkij.p.1[i, ] * dist.p.2[i, ]
    
    ## get aligned contacts of i
    aligned.ENMkij.i = ENMkij.p.1[i, ]
    aligned.contacts.i = which(aligned.ENMkij.i == 1)
    
    ## get distances between contacts 
    dist.p.1.contacts.i = as.vector(upper.tri(dist.p.1[aligned.contacts.i, aligned.contacts.i]))
    dist.p.2.contacts.i = as.vector(upper.tri(dist.p.2[aligned.contacts.i, aligned.contacts.i]))
    
    ## get complete distances
    complete.dist.p.1.i = c(dist.p.1.i, dist.p.1.contacts.i)
    complete.dist.p.2.i = c(dist.p.2.i, dist.p.2.contacts.i)
    
    ## calculate the local score of the site i
    local.score.i = mean((complete.dist.p.1.i - complete.dist.p.2.i) ^ 2)
    local.score[aligned.p.1.index[i]] = local.score.i
  }

  # normal modes
  ## calculate complete ENM Keff of p.1
  ENMKeff.p.1 = CalculateENMKeff(matrix(r.p.1, nrow = 3),
                                 aligned.p.1.index, 
                                 not.aligned.p.1.index,
                                 R0,
                                 tolerance,
                                 K.analysis)
  
  ## get eigenvalues and eigenvectors of p.1
  va = ENMKeff.p.1$va
  ve = ENMKeff.p.1$ve
  
  ## correct ve for K.analysis == "K"
  if (K.analysis == "K") {
    ve = ve[aligned.p.1.index.3N,]
  }
  
  ## calculate measures of variability between p.1 and p.2
  Pn = ((t(ve) %*% dr) ^ 2) / (sum((t(ve) %*% dr) ^ 2))
  
  ## dynamical variability
  
  dynamical.variability = CalculateDynamicalVariability(r.p.1,
                                                        r.p.2,
                                                        aligned.p.1.index,
                                                        aligned.p.2.index,
                                                        not.aligned.p.1.index,
                                                        not.aligned.p.2.index,
                                                        R0, 
                                                        tolerance, 
                                                        K.analysis)
  square.dif.MSF = dynamical.variability$square.dif.MSF
  nH = dynamical.variability$nH
  nR = dynamical.variability$nR
  
  # create a list for the output
  output = list(             "dr" = dr,
                             "va" = va,
                             "ve" = ve,
                             "Pn" = Pn,
                     "dr.squarei" = dr.squarei,
         "dr.squarei.windows.rot" = dr.squarei.windows.rot,
"dr.squarei.windows.contacts.rot" = dr.squarei.windows.contacts.rot,
                    "local.score" = local.score,
                 "square.dif.MSF" = square.dif.MSF,
                             "nH" = nH,
                             "nR" = nR)
  output
}

