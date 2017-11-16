# remove objects
rm(list = ls())

### INPUT ###

# set parameters
n.mut.p = 10 # set the numer of mutants per experimental protein

# set Elastic Network Model: "ANM" or "pfANM"
model <- "ANM"

# data dir
data.dir <- "DATA"

# output dir
out.dir <- "OUT/out_subset_CA_ANM"

# figures dir
fig.dir <- "OUT/mutants_score"

### PROGRAM ###

# load packages
library(bio3d) 
library(seqinr) 
library(ggplot2)

# general parameters
tolerance = 1e-10

# function filenames
ReadCA.fname <- "FUNCTIONS/ReadCA.R"
ReadHeme.fname <- "FUNCTIONS/ReadHeme.R"
CalculateENMKeff.fname <- "FUNCTIONS/CalculateENMKeff.R"
CalculateENMK.fname <- "FUNCTIONS/CalculateENMK.R"
if (model == "ANM") {
  CalculateKij.fname <- "FUNCTIONS/CalculateKij.R"
  CalculateForce.fname <- "FUNCTIONS/CalculateForce.R"
}
if (model == "pfANM") {
  CalculateKij.fname <- "FUNCTIONS/CalculateKijPFANM.R"
  CalculateForce.fname <- "FUNCTIONS/CalculateForcePFANM.R"
}

# source functions
source(ReadCA.fname) 
source(ReadHeme.fname)
source(CalculateENMKeff.fname)
source(CalculateENMK.fname)
source(CalculateKij.fname)
source(CalculateForce.fname)

# read input file
input.fname <- file.path("input_MainProgram.csv")
input <- read.csv(input.fname)

# start a loop for each row of the input
for (f in (1:nrow(input))) {
  print(f)

  family <- as.character(input$family)[f]
  p.ref <- as.character(input$p.ref)[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]
  fmax = input$fmax[f] 
  R0 = input$R0.CA[f]
  heme <- input$heme[f]
  
  ### READ EXPERIMENTAL DATA ###

  # Read score
  score.fname = file.path(out.dir, paste(p.ref, "_consurf.csv", sep = ""))
  score = read.csv(score.fname, sep =";")$SCORE
  
  # norm score [0, 1]
  rangos = range(score)
  maxMin = rangos[2] - rangos[1]
  norm.score = scale(score, rangos[1], ifelse(maxMin > 0, maxMin, 1))
  
  # get inverse scores
  norm.score.inv = 1 - norm.score
  
  # filenames
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 

  # read PDB of p.ref
  pdb = ReadCA(pdbs.fname, chain.p.ref)
  r.p.ref = pdb$xyz.calpha
  n.aa = pdb$n.sites
  
  # get heme coordinates, add them to CAB4s coordinates and calculate the new number of sites
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdbs.fname, chain.p.ref)
    r.p.ref = cbind(r.p.ref, r.heme)
    n.sites = ncol(r.p.ref)
  } else {
    n.sites = n.aa
  }
  
  # read identity
  m.identity = read.csv(file.path(out.dir, paste(family, "_out_m.identity.csv", sep = "")))$V1
  n.sites.mut.exp = n.sites - (m.identity * n.sites / 100)
  n.prot = length(m.identity)
  mean.n.sites.mut.exp = round(mean(n.sites.mut.exp))
  
  # read experimental data
  dri2.exp.fname = file.path(out.dir, paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.exp.norm.dr.squarei.csv", sep = ""))
  dri2.exp = read.csv(dri2.exp.fname, header = T)
  MSDi.exp = (colMeans(dri2.exp, na.rm = T))
  sd.exp = apply(dri2.exp, 2, sd, na.rm = T)
  
  ### CALCULATE ENM OF THE REFERENCE PROTEIN ###
  
  # calculate K of p.ref
  ENMK.p.ref = CalculateENMK(r.p.ref, CalculateKij, R0, tolerance) 
  
  ### GENERATE MUTANTS SCORE ###
  
  # create a matrix to save coordinates of each mutant
  m.dr.mut = matrix(0, nrow = 3 * n.sites, ncol = n.prot * n.mut.p)
  m.dr.squarei = matrix(0, nrow = n.sites, ncol = n.prot * n.mut.p)
  
  # create vectors to save simulated and accepted mutations of each mutant
  v.run = c()
  v.n.mut = c()
  v.site = c()
  v.p.accept = c()
  v.accept = c()
  v.n.accept.mut = c()
  
  # sart a counter for simulated proteins
  run = 0
  
  # start a loop for each protein of the family
  for (P in seq(n.prot * n.mut.p)) {
    print(P)
    # count the mut number
    run = run + 1
    
    # create a vector with possible sites to mutate
    sites.to.mutate = seq(n.aa) 
    
    # start a counter for simulated and accepted mutants
    n.mut = 0
    n.accept.mut = 0
    
    # start "dr.tot"
    dr.tot = 0
    
    # start a loop to generate mutations
    while (n.accept.mut < mean.n.sites.mut.exp) { 
      
      v.run = c(v.run, run)
      n.mut = n.mut + 1
      v.n.mut = c(v.n.mut, n.mut)
      
      # get the index to mutate between 1 and the number of aminoacids
      mut.index = sample(sites.to.mutate, 1) 
      v.site = c(v.site, mut.index)
      
      # calculate forces
      force = CalculateForce(mut.index, r.p.ref, ENMK.p.ref$kij, fmax) 
      f = force$f
      
      # calculate the acceptance probability of the mutation
      p.accept.mut = norm.score[mut.index]
      v.p.accept = c(v.p.accept, p.accept.mut)
      
      if (p.accept.mut >= runif(1, 0, 1)) {
        
        v.accept = c(v.accept, 1)
        
        # add to accepted mutants
        n.accept.mut = n.accept.mut + 1
        
        # calculate "dr.mut"
        dr.mut = ENMK.p.ref$cov %*% f
        dr.squarei = colSums(matrix(dr.mut, nrow = 3) ^ 2) 
        
        # keep "dr.mut" in a matrix
        m.dr.mut[, P] = dr.mut
        m.dr.squarei[, P] = dr.squarei
        
        # remove "mut.index" from "sites.to.mutate" in order to not to mutate the same site more than once
        sites.to.mutate = sites.to.mutate[sites.to.mutate != mut.index]
      }else{
        v.accept = c(v.accept, 0)
      }
    }
  }
  
  # save matrices
  m.dr.mut.1 = m.dr.mut
  m.dr.squarei.1 = m.dr.squarei
  
  ### GENERATE MUTANTS INV SCORE ###
  
  # create a matrix to save coordinates of each mutant
  m.dr.mut = matrix(0, nrow = 3 * n.sites, ncol = n.prot * n.mut.p)
  m.dr.squarei = matrix(0, nrow = n.sites, ncol = n.prot * n.mut.p)
  
  # create vectors to save simulated and accepted mutations of each mutant
  v.run = c()
  v.n.mut = c()
  v.site = c()
  v.p.accept = c()
  v.accept = c()
  v.n.accept.mut = c()
  
  # sart a counter for simulated proteins
  run = 0
  
  # start a loop for each protein of the family
  for (P in seq(n.prot * n.mut.p)) {
    print(P)
    # count the mut number
    run = run + 1
      
    # create a vector with possible sites to mutate
    sites.to.mutate = seq(n.aa) 
      
    # start a counter for simulated and accepted mutants
    n.mut = 0
    n.accept.mut = 0
      
    # start "dr.tot"
    dr.tot = 0
      
    # start a loop to generate mutations
    while (n.accept.mut < mean.n.sites.mut.exp) { 
        
      v.run = c(v.run, run)
      n.mut = n.mut + 1
      v.n.mut = c(v.n.mut, n.mut)
        
      # get the index to mutate between 1 and the number of aminoacids
      mut.index = sample(sites.to.mutate, 1) 
      v.site = c(v.site, mut.index)
        
      # calculate forces
      force = CalculateForce(mut.index, r.p.ref, ENMK.p.ref$kij, fmax) 
      f = force$f

      # calculate the acceptance probability of the mutation
      p.accept.mut = norm.score.inv[mut.index]
      v.p.accept = c(v.p.accept, p.accept.mut)
        
      if (p.accept.mut >= runif(1, 0, 1)) {
          
        v.accept = c(v.accept, 1)
          
        # add to accepted mutants
        n.accept.mut = n.accept.mut + 1
          
        # calculate "dr.mut"
        dr.mut = ENMK.p.ref$cov %*% f
        dr.squarei = colSums(matrix(dr.mut, nrow = 3) ^ 2) 
          
        # keep "dr.mut" in a matrix
        m.dr.mut[, P] = dr.mut
        m.dr.squarei[, P] = dr.squarei
          
        # remove "mut.index" from "sites.to.mutate" in order to not to mutate the same site more than once
        sites.to.mutate = sites.to.mutate[sites.to.mutate != mut.index]
      }else{
        v.accept = c(v.accept, 0)
      }
    }
  }

  # save matrices
  m.dr.mut.2 = m.dr.mut
  m.dr.squarei.2 = m.dr.squarei

  ### ANALYSIS ###
  
  # plot MSD of mutants
  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax, "_MSD_mutants.png", sep = "")))
  layout(matrix(1:2, 2, 1))
  MSD.P.1 = colMeans(m.dr.squarei.1)
  MSD.1 = mean(MSD.P.1)
  plot(MSD.P.1, xlab = "mutant", ylab = "MSD", main = "score", ylim = c(0, 0.04))
  legend("topright", paste ("mean = ", round(MSD.1, digits = 4)), bty = "n")
  abline(MSD.1, 0, col = "green")
  MSD.P.2 = colMeans(m.dr.squarei.2)
  MSD.2 = mean(MSD.P.2)
  plot(MSD.P.2, xlab = "mutant", ylab = "MSD", main = "inv score", ylim = c(0, 0.04))
  legend("topright", paste ("mean = ", round(MSD.2, digits = 4)), bty = "n")
  abline(MSD.2, 0, col = "red")
  dev.off()
  
  # norm matrices
  m.dr.squarei.1.norm = t(t(m.dr.squarei.1) / colMeans(m.dr.squarei.1))
  m.dr.squarei.2.norm = t(t(m.dr.squarei.2) / colMeans(m.dr.squarei.2))

  # calculate means and sd
  MSDi.1 = rowMeans(m.dr.squarei.1.norm)
  MSDi.2 = rowMeans(m.dr.squarei.2.norm)

  sd.1 = apply(m.dr.squarei.1.norm, 1, sd)
  sd.2 = apply(m.dr.squarei.2.norm, 1, sd)

  # plot MSDi, CN and SCORE
  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax, "_MSDi_CN_SCORE.png", sep = "")))
  layout(matrix(1:3, 3, 1))
  plot(MSDi.1, type = "l", col = "green", ylab = "MSDi +/- SD", xlab = "site")
  lines((MSDi.1 + (sd.1)/sqrt(n.mut.p * n.prot)), col = "green")
  lines((MSDi.1 - (sd.1)/sqrt(n.mut.p * n.prot)), col = "green")
  lines(MSDi.2, col = "red")
  lines((MSDi.2 + (sd.2)/sqrt(n.mut.p * n.prot)), col = "red")
  lines((MSDi.2 - (sd.2)/sqrt(n.mut.p * n.prot)), col = "red")
  plot((1/colSums(ENMK.p.ref$kij)), xlab = "site", ylab = "1/CN", type = "l")
  plot(norm.score, type = "l", col = "green", xlab = "site", ylab = "score")
  lines(norm.score.inv, col = "red")
  dev.off()
  
  # plot MSDi SCORE vs INVSCORE
  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax, "_MSDi_SCORE_INVSCORE.png", sep = "")))
  layout(matrix(1, 1, 1))
  plot(x = MSDi.1, y = MSDi.2, xlab = "MSDi score", ylab = "MSDi inv score", xlim = c(0.25, 3.5), ylim = c(0.25, 3.5))
  legend("topleft", paste ("r = ", round(cor(MSDi.1, MSDi.2, method = "spearman"), digits = 3)), bty = "n")
  abline(0, 1)
  dev.off()
  
  # plot MSDi exp vs MSDi SCORE of INVSCORE 
  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax, "_MSDi_EXP_SCORE_INVSCORE.png", sep = "")))
  layout(matrix(1:2, 2, 1))
  plot(x = sqrt(MSDi.1), y = sqrt(MSDi.exp), xlab = "sqrt MSDi score", ylab = "sqrt MSDi exp")
  legend("topleft", paste ("r = ", round(cor(sqrt(MSDi.1), sqrt(MSDi.exp), method = "spearman"), digits = 3)), bty = "n")
  plot(x = sqrt(MSDi.2), y = sqrt(MSDi.exp), xlab = "sqrt MSDi inv score", ylab = "sqrt MSDi exp")
  legend("topleft", paste ("r = ", round(cor(sqrt(MSDi.2), sqrt(MSDi.exp), method = "spearman"), digits = 3)), bty = "n")
  dev.off()
}  