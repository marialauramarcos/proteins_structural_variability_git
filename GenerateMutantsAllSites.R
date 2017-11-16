# remove objects
rm(list = ls())

### INPUT ###

# set parameters
n.mut.s = 50 # set number of mutants per site
fmax = 2 # set fmax for the forces

# set Elastic Network Model: "ANM" or "pfANM"
model <- "ANM"

# data dir
data.dir <- "DATA"

# figures dir
fig.dir <- "OUT/out_mutants_all_sites"

### PROGRAM ###

# load packages
library(bio3d) 
library(ggplot2)
library(plyr)

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

# load functions
source(ReadCA.fname) 
source(ReadHeme.fname)
source(CalculateENMKeff.fname)
source(CalculateENMK.fname)
source(CalculateKij.fname)
source(CalculateForce.fname)

# read input file
input.fname <- file.path("input_MainProgram.csv")
input <- read.csv(input.fname)

# build matrices to save data
s.rm.md.f = matrix(nrow = 300, ncol = nrow(input)) #no protein in the experimental dataset has more than 300 aa
s.cm.md.f = matrix(nrow = 300, ncol = nrow(input))

CN.f = matrix(nrow = 300, ncol = nrow(input))
MSF.f = matrix(nrow = 300, ncol = nrow(input))
RSA.f = matrix(nrow = 300, ncol = nrow(input))

# start a loop for each row of the input
for (p in (1:nrow(input))) {
  print(p)
  
  # read p
  family <- as.character(input$family)[p]
  p.ref <- as.character(input$p.ref)[p]
  chain.p.ref <- as.character(input$chain.p.ref)[p]
  R0 = input$R0.CA[p]
  heme <- input$heme[p]

  # pdbs filename
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 

  # read PDB of p.ref
  pdb = ReadCA(pdbs.fname, chain.p.ref)
  r.p.ref = pdb$xyz.calpha
  n.aa = pdb$n.sites
  
  # only for globins - get heme coordinates, add them to CAs coordinates and calculate the new number of sites
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdbs.fname, chain.p.ref)
    r.p.ref = cbind(r.p.ref, r.heme)
    n.sites = ncol(r.p.ref)
  } else {
    n.sites = n.aa
  }
  
  # calculate K of p.ref
  ENMK.p.ref = CalculateENMK(r.p.ref, CalculateKij, R0, tolerance) 

  # calculate CN
  CN = colSums(ENMK.p.ref$kij)
  CN.f[1:length(CN), p] = CN

  # calculate RSA
  ## filenames
  asa.fname = file.path(data.dir, paste(p.ref, "_dssp.csv", sep = ""))
  surfase.area.fname = file.path(data.dir, "surfase_area.csv")

  ## get asa or acc
  asa = read.csv(asa.fname, sep = ";")
  ACC = asa$ACC
  AA = asa$AA

  ## get the surfase area of AA
  surfase.area = read.csv(surfase.area.fname, sep = ";")
  AA.CODE = surfase.area$X1.LETTER.AA.CODE
  SURFASE.AREA = surfase.area$SURFASE.AREA

  ## calculate RSA
  RSA = c()
  for (i in (1:length(AA))) {
    RSA[i] = ACC[i]/SURFASE.AREA[which(AA.CODE == as.character(AA[i]))]
  }
  RSA.f[1:length(RSA), p] = RSA

  # calculate MSFi
  ## get cov matrix
  cov.p.ref = ENMK.p.ref$cov

  ## get the diagonal of the cov matrix
  diag.p.ref = diag(cov.p.ref)

  ## calculate the factor to split the diagonal
  factor = sort(rep(seq(1:n.sites), 3))

  ## split the diagonal
  s.diag.p.ref = split(diag.p.ref, factor)

  ## create a matrix to save the data
  MSF.p.ref = matrix(nrow = 1, ncol = n.sites)

  ## start a loop to calculate sums for each site
  for (i in (1:n.sites)) {
    MSF.p.ref[, i] = sum(unlist(s.diag.p.ref[i]), use.names = F)
  }

  ## transform to a vector and save
  MSF.p.ref = as.vector(MSF.p.ref)
  MSF.f[1:length(MSF.p.ref), p] = MSF.p.ref

  # create matrices to save data of each mutant
  m.dr.mut = matrix(nrow = 3 * n.sites, ncol = n.sites * n.mut.s)
  m.dr.squarei = matrix(nrow = n.sites, ncol = n.sites * n.mut.s)

  # start a loop for each site of the protein
  for (S in (1:n.sites)) {
  
    print(S)
  
    # start a loop to generate "n.mut.s" mutants for each S
    for(mut in seq(n.mut.s)) {
      
      print(mut)
    
      # calculate forces
      force = CalculateForce(S, r.p.ref, ENMK.p.ref$kij, fmax) 
      f = force$f

      # calculate "dr.mut"
      dr.mut = ENMK.p.ref$cov %*% f
      dr.squarei = colSums(matrix(dr.mut, nrow = 3) ^ 2) 
      
      # keep "dr.mut" and "dr.squarei"
      m.dr.mut[, n.mut.s * S - (n.mut.s - mut)] = dr.mut
      m.dr.squarei[, n.mut.s * S - (n.mut.s - mut)] = dr.squarei
    }
  }

  # set values for dr.squarei comparisons
  site.mut = sort(rep(seq(1:n.sites), (n.sites * n.mut.s))) # mutated site
  mut = sort(rep(seq(1:(n.mut.s * n.sites)), (n.sites))) # index of mutant
  site = rep(sort(rep(seq(1:(n.sites)))), (n.mut.s * n.sites)) # site of each mutant

  # bluid a data frame
  df.dr.mut = data.frame(site.mut, mut, site, SD = as.vector(m.dr.squarei)) 
  
  # norm SD
  df.norm = ddply(df.dr.mut, c("site.mut", "mut"), mutate,
              "norm.SD" = SD/mean(SD))

  # calculate means for each mutated site
  d = ddply(df.norm, c("site.mut", "site"),
            function(x) {
              mean.SD = mean(x$SD, na.rm = T)
              mean.norm.SD = mean(x$norm.SD, na.rm = T)
              data.frame(mean.SD, mean.norm.SD)})

  # write a file with d
  write.csv(d, file = paste(p.ref, "_fmax", fmax,"_sensitivity.csv", sep = ""))

  # transform d
  m.d = matrix(d$mean.SD, nrow = n.sites)

  # calculate means
  s.rm.md = scale(rowMeans(m.d))
  s.cm.md = scale(colMeans(m.d))

  s.rm.md.f[1:length(s.rm.md), p] = s.rm.md
  s.cm.md.f[1:length(s.cm.md), p] = s.cm.md

  # plots
  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax,"_dr2_site.png", sep = "")))
  plot(s.rm.md, type = "l", xlab = "site", ylab = "scale dr2", ylim = c(-2, 6), main = "black: dr2_i_x, red: dr2_x_j")
  lines(s.cm.md, col = "red")
  dev.off()

  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax,"_dr2_dr2.png", sep = "")))
  plot(x = s.rm.md, y = s.cm.md, xlab = "scale dr2_i_x", ylab =  "scale dr2_x_j", ylim = c(-2, 5.5), xlim = c(-2, 5.5))
  legend("topleft", paste ("r =", round(cor(s.rm.md, s.cm.md, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(s.cm.md ~ s.rm.md)
  #abline(lm)
  dev.off()
}

### PLOTS ALL FAMILIES ###

# plot dr2 vs site all families
png(filename = file.path(fig.dir, paste("fmax", fmax,"_dr2_site_all_families.png", sep = "")))
layout(matrix(1:8, 4, 2))
for (p in (1:nrow(input))){
  s.rm.md = s.rm.md.f[!is.na(s.rm.md.f[, p]), p]
  s.cm.md = s.cm.md.f[!is.na(s.cm.md.f[, p]), p]
  p.ref <- as.character(input$p.ref)[p]
  plot(s.rm.md, type = "l", xlab = "site", ylab = "scale dr2", ylim = c(-2, 6), main = "black: dr2_i_x, red: dr2_x_j")
  lines(s.cm.md, col = "red")
  legend("topleft", p.ref, bty = "n", col = "red")
}
dev.off()

# plot dr2 vs dr2 all families
png(filename = file.path(fig.dir, paste("fmax", fmax,"_dr2_dr2_all_families_1.png", sep = "")))
layout(matrix(1:4, 2, 2))
for (p in (1:4)) {
  s.rm.md = s.rm.md.f[!is.na(s.rm.md.f[, p]), p]
  s.cm.md = s.cm.md.f[!is.na(s.cm.md.f[, p]), p]
  p.ref <- as.character(input$p.ref)[p]
  plot(x = s.rm.md, y = s.cm.md, xlab = "scale dr2_i_x", ylab =  "scale dr2_x_j", ylim = c(-2, 5.5), xlim = c(-2, 5.5))
  legend("topleft", paste ("r =", round(cor(s.rm.md, s.cm.md, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(s.cm.md ~ s.rm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")
}
dev.off()

png(filename = file.path(fig.dir, paste("fmax", fmax,"_dr2_dr2_all_families_2.png", sep = "")))
layout(matrix(1:4, 2, 2))
for (p in (5:8)){
  s.rm.md = s.rm.md.f[!is.na(s.rm.md.f[, p]), p]
  s.cm.md = s.cm.md.f[!is.na(s.cm.md.f[, p]), p]
  p.ref <- as.character(input$p.ref)[p]
  plot(x = s.rm.md, y = s.cm.md, xlab = "scale dr2_i_x", ylab =  "scale dr2_x_j", ylim = c(-2, 5.5), xlim = c(-2, 5.5))
  legend("topleft", paste ("r =", round(cor(s.rm.md, s.cm.md, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(s.cm.md ~ s.rm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")
}
dev.off()

# plot dr vs CN, MSF and RSA all families

for (p in (1:nrow(input))){
  s.rm.md = s.rm.md.f[!is.na(s.rm.md.f[, p]), p]
  s.cm.md = s.cm.md.f[!is.na(s.cm.md.f[, p]), p]
  CN = CN.f[!is.na(CN.f[, p]), p]
  MSF = MSF.f[!is.na(MSF.f[, p]), p]
  RSA = RSA.f[!is.na(RSA.f[, p]), p]
  
  p.ref <- as.character(input$p.ref)[p]

  png(filename = file.path(fig.dir, paste(p.ref, "_fmax", fmax, "_dr2_CN_MSF_RSA.png", sep = "")))
  layout(matrix(1:6, 2, 3, byrow = TRUE))
  
  plot(x = s.rm.md, y = 1/CN, xlab = "scale dr2_i_x", ylab = "1/CN")
  legend("topleft", paste ("r =", round(cor(s.rm.md, CN, method = "spearman"), digits = 3)), bty = "n")
  lm = lm((1/CN) ~ s.rm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")
  
  plot(x = s.rm.md, y = RSA, xlab = "scale dr2_i_x", ylab =  "RSA")
  legend("topleft", paste ("r =", round(cor(s.rm.md, RSA, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(RSA ~ s.rm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")

  plot(x = s.rm.md, y = MSF, xlab = "scale dr2_i_x", ylab =  "MSF")
  legend("topleft", paste ("r =", round(cor(s.rm.md, MSF, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(MSF ~ s.rm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")

  plot(x = s.cm.md, y = 1/CN, xlab = "scale dr2_x_j", ylab =  "1/CN")
  legend("topleft", paste ("r =", round(cor(s.cm.md, CN, method = "spearman"), digits = 3)), bty = "n")
  lm = lm((1/CN) ~ s.cm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")

  plot(x = s.cm.md, y = RSA, xlab = "scale dr2_x_j", ylab =  "RSA")
  legend("topleft", paste ("r =", round(cor(s.cm.md, RSA, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(RSA ~ s.cm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")

  plot(x = s.cm.md, y = MSF, xlab = "scale dr2_x_j", ylab =  "MSF")
  legend("topleft", paste ("r =", round(cor(s.cm.md, MSF, method = "spearman"), digits = 3)), bty = "n")
  lm = lm(MSF ~ s.cm.md)
  #abline(lm)
  legend("topright", p.ref, bty = "n", col = "red")
  
  dev.off()
}