# Description:
#
# This function analyzes experimental and theoretical proteins calculating measures of structural variability in cartesian 
# coordinates and projected on the normal modes of a reference protein. 
# 
# Usage:
#
# AnalyzeExperimentalTheoretical <- function(family, p.ref, chain.p.ref, n.mut.p, R0, rotate = TRUE/FALSE,
# heme = TRUE/FALSE, K.analysis, data.dir, out.dir, mut.fname.id, analysis.fname.id, tolerance)
#
#  Args:
#    - family: the family of the protein to mutate. It can be "globins", "serinProteases", 
#    "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#    - p.ref: the reference protein of the family.
#    - chain.p.ref: the chain of p.ref in the pdb file obtained from Homstrad.
#    - n.mut.p: the number of mutants to generate for each member of the family. For example, if the family has 20 
#    members, the program generates n.mut.p x 20 mutants.
#    - R0: the cut-off for the ANM ("Anisotropic Network Model") that represents the proteins.
#    - rotate: it can be "TRUE" or "FALSE". If it is "TRUE", r.p.2 is rotaded in order to minimize RMSD with r.p.ref.
#    - heme: argument for globins. It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group. 
#    - K.analysis: It can be "K" or "Keff". For "K" or "Keff", the analysis is based on normal modes of "K" or "Keff"
#    respectibly.
#    - data.dir: directory of the data. It must contain the dataset ("data.dir/family_dataset.csv") and the pdb file 
#    obtained from Homstrad ("data.dir/family_coordinates.csv").
#    - out.dir: directory of the output. It must contain output files generated with AnalyzeFamily() and GenerateMutants().
#    The output of this function is also saved in out.dir.
#    - mut.fname.id: ID of filnames of mutant proteins.
#    - analysis.fname.id: ID of output filenames.
#    - tolerance: 0 tolerance.
#
#  Required libraries
#    {Bio3d}
#
#  Required functions
#    ReadCA()
#    ReadHeme()
#    CalculateSideChainCM()
#    CalculateVariability()
#    CalculateENMKeff()
#    CalculateENMK()
#    WindowsRMSD()

AnalyzeExperimentalTheoreticalCM <- function(family,
                                             p.ref,
                                             chain.p.ref,
                                             n.mut.p,
                                             R0, 
                                             rotate,
                                             heme,
                                             K.analysis,
                                             data.dir,
                                             out.dir,
                                             mut.fname.id, 
                                             analysis.fname.id,
                                             tolerance) {
  ### GENERAL ###
  
  # filenames
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))
  pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = "")) 
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  
  m.n.aligned.fname <- file.path(out.dir, paste(family, "_out_m.n.aligned.csv", sep = ""))
  m.aligned.p.ref.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.p.ref.index.csv", sep = ""))
  m.aligned.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.p.2.index.csv", sep = ""))
  m.not.aligned.p.ref.index.fname <- file.path(out.dir, paste(family, "_out_m.not.aligned.p.ref.index.csv", sep = ""))
  m.not.aligned.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.not.aligned.p.2.index.csv", sep = ""))
  
  # read the dataset
  dataset <- read.csv(dataset.fname)
  pdbid.dataset <- as.character(dataset$pdbid)
  chain <- as.character(dataset$chain)
  n.prot = length(pdbid.dataset)
  
  ### THEORETICAL ###
  print("analyzing theoretical...")
  
  # filenames
  theo.r.p.ref.fname <- file.path(out.dir, paste(mut.fname.id, "_out_r.p.ref.csv", sep = ""))
  m.r.mut.fname <- file.path(out.dir, paste(mut.fname.id, "_out_m.r.mut.csv", sep = ""))
  
  # read coordinates
  theo.r.p.ref = read.csv(theo.r.p.ref.fname)$x
  m.r.mut = read.csv(m.r.mut.fname)
  n.sites.p.ref = length(theo.r.p.ref)/3
  n.mut = ncol(m.r.mut)
  
  r.CM.p.ref = CalculateSideChainCM(pdb.fname, chain.p.ref)
  n.CM = ncol(r.CM.p.ref)
  
  theo.r.p.ref = c(theo.r.p.ref, as.vector(r.CM.p.ref))
  n.sites.tot = length(theo.r.p.ref)/3
  
  # create matrices to save measures of variability of each mutant
  m.theo.Pn = matrix(nrow = n.mut, ncol = 3 * n.sites.p.ref)
  m.theo.va = matrix(nrow = n.mut, ncol = 3 * n.sites.p.ref)
  m.theo.dr.squarei = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.dr.squarei.windows.rot = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.dr.squarei.windows.contacts.rot = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.norm.dr.squarei = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.norm.dr.squarei.windows.rot = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.norm.dr.squarei.windows.contacts.rot = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.local.score = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.square.dif.MSF = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.nH = matrix(nrow = n.mut, ncol = 3 * n.sites.p.ref)
  m.theo.nR = matrix(nrow = n.mut, ncol = 3 * n.sites.p.ref)
  
  # start a loop to analyze each mutant
  for (mut in (1:(n.mut))) {
    print(c(mut))
      
    # get theo.r.p.2
    theo.r.p.2 = m.r.mut[, mut]
    
    # calculate measures of variavility
    theo.variability = CalculateVariability(as.vector(theo.r.p.ref),
                                            as.vector(theo.r.p.2),
                                            n.sites.p.ref,
                                            n.sites.p.ref,
                                            seq(1, n.sites.p.ref),                                               
                                            seq(1, n.sites.p.ref), 
                                            seq((n.sites.p.ref + 1), n.sites.tot),
                                            seq((n.sites.p.ref + 1), n.sites.tot),
                                            R0,
                                            rotate,
                                            K.analysis,
                                            tolerance)
      
    m.theo.va[mut, 1:length(theo.variability$va)] = theo.variability$va
    m.theo.Pn[mut, 1:length(theo.variability$Pn)] = theo.variability$Pn
    m.theo.dr.squarei[mut, 1:length(theo.variability$dr.squarei)] = theo.variability$dr.squarei
    m.theo.dr.squarei.windows.rot[mut, 1:length(theo.variability$dr.squarei)] = theo.variability$dr.squarei.windows.rot
    m.theo.dr.squarei.windows.contacts.rot[mut, 1:length(theo.variability$dr.squarei)] = theo.variability$dr.squarei.windows.contacts.rot
    m.theo.local.score[mut, 1:length(theo.variability$local.score)] = theo.variability$local.score
    m.theo.square.dif.MSF[mut, 1:length(theo.variability$square.dif.MSF)] = theo.variability$square.dif.MSF
    m.theo.nH[mut, 1:length(theo.variability$nH)] = theo.variability$nH
    m.theo.nR[mut, 1:length(theo.variability$nR)] = theo.variability$nR
    
    # calculate norm.dr.squarei
    m.theo.norm.dr.squarei[mut, ] = m.theo.dr.squarei[mut, ]/ mean(m.theo.dr.squarei[mut, ], na.rm = T)
    m.theo.norm.dr.squarei.windows.rot[mut, ] = m.theo.dr.squarei.windows.rot[mut, ]/ mean(m.theo.dr.squarei.windows.rot[mut, ], na.rm = T)
    m.theo.norm.dr.squarei.windows.contacts.rot[mut, ] = m.theo.dr.squarei.windows.contacts.rot[mut, ]/ mean(m.theo.dr.squarei.windows.contacts.rot[mut, ], na.rm = T)
  }
  
  # create files to save the data
  write.csv(m.theo.va, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.va.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.Pn, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.Pn.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.norm.dr.squarei, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.norm.dr.squarei.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.norm.dr.squarei.windows.rot, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.norm.dr.squarei.windows.rot.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.norm.dr.squarei.windows.contacts.rot, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.norm.dr.squarei.windows.contacts.rot.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.local.score, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.local.score.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.square.dif.MSF, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.square.dif.MSF.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.nH, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.nH.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.nR, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.nR.csv", sep = "")), row.names = FALSE)
  
  ### EXPERIMENTAL ###
  print("analyzing experimental...")
  
  # read pdb of exp.p.ref
  exp.r.p.ref = theo.r.p.ref

  if(heme == "TRUE") {
    n.aa.p.ref = n.sites.p.ref - 5 
  } else {
    n.aa.p.ref = n.sites.p.ref
  }
  
  # read indexes files
  m.n.aligned = read.csv(m.n.aligned.fname)
  m.aligned.p.ref.index = read.csv(m.aligned.p.ref.index.fname)
  m.aligned.p.2.index = read.csv(m.aligned.p.2.index.fname)
  m.not.aligned.p.ref.index = read.csv(m.not.aligned.p.ref.index.fname)
  m.not.aligned.p.2.index = read.csv(m.not.aligned.p.2.index.fname)
  
  # create matrices to save measures of variability of each mutant
  m.exp.Pn = matrix(nrow = n.prot, ncol = 3 * (n.sites.p.ref)) 
  m.exp.va = matrix(nrow = n.prot, ncol = 3 * (n.sites.p.ref))
  m.exp.dr.squarei = matrix(nrow = n.prot, ncol = (n.sites.p.ref))
  m.exp.dr.squarei.windows.rot = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.exp.dr.squarei.windows.contacts.rot = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.exp.norm.dr.squarei = matrix(nrow = n.prot, ncol = (n.sites.p.ref))
  m.exp.norm.dr.squarei.windows.rot = matrix(nrow = n.prot, ncol = n.sites.p.ref)
  m.exp.norm.dr.squarei.windows.contacts.rot = matrix(nrow = n.prot, ncol = n.sites.p.ref)
  m.exp.local.score = matrix(nrow = n.prot, ncol = n.sites.p.ref)
  m.exp.square.dif.MSF = matrix(nrow = n.prot, ncol = n.sites.p.ref)
  m.exp.nH = matrix(nrow = n.prot, ncol = 3 * n.sites.p.ref)
  m.exp.nR = matrix(nrow = n.prot, ncol = 3 * n.sites.p.ref)
  
  # start a loop to evaluate each protein of the family
  for (P in (1:n.prot)) {
    print(P)
    
    # get aligned and not aligned indexes
    n.aligned = as.numeric(m.n.aligned[P, ])
    aligned.p.ref.index = as.numeric(m.aligned.p.ref.index[P, !is.na(m.aligned.p.ref.index[P, ])])
    aligned.p.2.index = as.numeric(m.aligned.p.2.index[P, !is.na(m.aligned.p.2.index[P, ])])
    not.aligned.p.ref.index = as.numeric(m.not.aligned.p.ref.index[P, !is.na(m.not.aligned.p.ref.index[P, ])])
    not.aligned.p.2.index = as.numeric(m.not.aligned.p.2.index[P, !is.na(m.not.aligned.p.2.index[P, ])])
    
    # read PDB of exp.p.2
    chain.p.2 <- chain[[P]]
    exp.pdb.p.2 = ReadCA(pdbs.fname, chain.p.2)
    exp.r.p.2 = exp.pdb.p.2$xyz.calpha
    exp.n.aa.p.2 = ncol(exp.r.p.2)
    
    # calculate heme coordinates, add them to CA´s coordinates and calculate the number of sites and not aligned indexes
    if (heme == "TRUE") {
      exp.r.heme.p.2 = ReadHeme(pdbs.fname, chain.p.2)
      exp.r.p.2 = cbind(exp.r.p.2, exp.r.heme.p.2)
      exp.n.sites.p.2 = ncol(exp.r.p.2)
      
      aligned.p.ref.index <- c(aligned.p.ref.index, t(seq((n.aa.p.ref + 1), n.sites.p.ref)))
      aligned.p.2.index <- c(aligned.p.2.index, t(seq((exp.n.aa.p.2 + 1), exp.n.sites.p.2)))
    } else {
      exp.n.sites.p.2 = exp.n.aa.p.2
    }
    
    # add CM to exp.p.2
    exp.r.CM.p.2 = CalculateSideChainCM(pdbs.fname, chain.p.2)

    exp.r.p.2 = c(exp.r.p.2, as.vector(exp.r.CM.p.2))
    n.sites.tot.p.2 = length(exp.r.p.2)/3
    
    # add indexes of CMs to not aligned
    not.aligned.p.ref.index <- c(not.aligned.p.ref.index, t(seq((n.sites.p.ref + 1), n.sites.tot)))
    not.aligned.p.2.index <- c(not.aligned.p.2.index, t(seq((exp.n.sites.p.2 + 1), n.sites.tot.p.2)))
    
    # calculate measures of variability
    exp.variability = CalculateVariability(as.vector(exp.r.p.ref), 
                                           as.vector(exp.r.p.2), 
                                           n.sites.p.ref,
                                           exp.n.sites.p.2,
                                           aligned.p.ref.index, 
                                           aligned.p.2.index, 
                                           not.aligned.p.ref.index,
                                           not.aligned.p.2.index,
                                           R0,
                                           rotate,
                                           K.analysis,
                                           tolerance)
    
    m.exp.va[P, 1:length(exp.variability$va)] = exp.variability$va
    m.exp.Pn[P, 1:length(exp.variability$Pn)] = exp.variability$Pn
    m.exp.local.score[P, 1:length(exp.variability$local.score)] = exp.variability$local.score
    m.exp.nH[P, 1:length(exp.variability$nH)] = exp.variability$nH
    m.exp.nR[P, 1:length(exp.variability$nR)] = exp.variability$nR
    
    exp.dr.squarei = exp.variability$dr.squarei
    exp.dr.squarei.windows.rot = exp.variability$dr.squarei.windows.rot
    exp.dr.squarei.windows.contacts.rot = exp.variability$dr.squarei.windows.contacts.rot
    exp.square.dif.MSF = exp.variability$square.dif.MSF
    
    for (i in (1:n.sites.p.ref)) {
      m.exp.dr.squarei[P, i] = matrix(exp.dr.squarei[aligned.p.ref.index == i], nrow = 1, ncol = 1)
      m.exp.dr.squarei.windows.rot[P, i] = matrix(exp.dr.squarei.windows.rot[aligned.p.ref.index == i], nrow = 1, ncol = 1)
      m.exp.dr.squarei.windows.contacts.rot[P, i] = matrix(exp.dr.squarei.windows.contacts.rot[aligned.p.ref.index == i], nrow = 1, ncol = 1)
      m.exp.square.dif.MSF[P, i] = matrix(exp.square.dif.MSF[aligned.p.ref.index == i], nrow = 1, ncol = 1)
    }
    
    # calculate norm.dr.squarei
    m.exp.norm.dr.squarei[P, ] = m.exp.dr.squarei[P, ]/ mean(m.exp.dr.squarei[P, ], na.rm = T)
    m.exp.norm.dr.squarei.windows.rot[P, ] = m.exp.dr.squarei.windows.rot[P, ]/ mean(m.exp.dr.squarei.windows.rot[P, ], na.rm = T)
    m.exp.norm.dr.squarei.windows.contacts.rot[P, ] = m.exp.dr.squarei.windows.contacts.rot[P, ]/ mean(m.exp.dr.squarei.windows.contacts.rot[P, ], na.rm = T)
  }  
  
  # create files to save the data
  write.csv(m.exp.va, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.va.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.Pn, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.Pn.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.norm.dr.squarei, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.norm.dr.squarei.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.norm.dr.squarei.windows.rot, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.norm.dr.squarei.windows.rot.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.norm.dr.squarei.windows.contacts.rot, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.norm.dr.squarei.windows.contacts.rot.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.local.score, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.local.score.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.square.dif.MSF, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.square.dif.MSF.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.nH, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.nH.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.nR, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.nR.csv", sep = "")), row.names = FALSE)
}

