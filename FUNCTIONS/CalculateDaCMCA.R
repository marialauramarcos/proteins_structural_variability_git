# Description:
#
# This function calculates distances between functional sites of a protein and the other sites of the protein.
# It calculates distances between both CAs and CMs of the side chains of aminoacids.
#
# Usage:
#
# CalculateDaCMCA(family,
#                 p.ref,
#                 chain.p.ref,
#                 heme = TRUE/FALSE,
#                 data.dir,
#                 out.dir)
#
#  Args:
#    - family: The family of p.ref.
#    - p.ref: The pdb code (pdbid) of the reference protein (example: "1a6m"). The protein must be a member of
#    the selected family.
#    - chain.p.ref: The chain of p.ref in the pdb file obtained from Homstrad.
#    - heme: argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    - data.dir: directory of the data.
#    - out.dir: directory for the output. 
#
#  Required libraries:
#    {Bio3d}
#
#  Required functions:
#    ReadCA()
#    ReadHeme()
#    CalculateSideChainCM()

  CalculateDaCMCA <- function(family,
                              p.ref,
                              chain.p.ref,
                              heme,
                              data.dir,
                              out.dir) {
  
  # read functional sites
  active = read.csv(file.path(data.dir, paste(p.ref, "_functionalSites.csv", sep = "")), sep = ";")$index
  active = active[!is.na(active)]
  
  # build pdb filename
  pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = ""))

  # get coordinates of p.ref
  r.CA.p.ref = ReadCA(pdb.fname, chain.p.ref)$xyz
  r.CM.p.ref = CalculateSideChainCM(pdb.fname, chain.p.ref)
  n.aa = ncol(r.CA.p.ref)
  
  # read reference pdb index
  ref.pdb = read.pdb(pdb.fname)
  inds = atom.select(ref.pdb, elety = "CA")
  ref.pdb.ca = ref.pdb$atom[inds$atom,]
  resno = ref.pdb.ca$resno
  
  # correct active site index
  active.site.index.n = c()
  for (i in (1:length(active))) {active.site.index.n = c(active.site.index.n, which(resno == active[i]))}
  active = active.site.index.n
  
  # get heme coordinates
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdb.fname, chain.p.ref)
    r.CA.p.ref = cbind(r.CA.p.ref, r.heme)
    r.CM.p.ref = cbind(r.CM.p.ref, r.heme)
  }
  
  # calculate the number of sites
  n.sites = ncol(r.CA.p.ref)  # it is different than n.aa in the case of heme = "TRUE"

  # calculate the distances between each site and each functional site (CA or CM)
  m.da.CA = matrix(nrow = n.sites, ncol = length(active))
  m.da.CM = matrix(nrow = n.sites, ncol = length(active))
  
  for (i in (active)) {
    for (j in (1:n.sites)) {
      m.da.CA[j, which(active == i)] = sqrt(sum((r.CA.p.ref[, j] - r.CA.p.ref[, i]) ^ 2))
      m.da.CM[j, which(active == i)] = sqrt(sum((r.CM.p.ref[, j] - r.CM.p.ref[, i]) ^ 2))
    }
  }
  
  m.da.CA = cbind(seq(1:n.sites), m.da.CA)
  colnames(m.da.CA) = c("site", paste("da", seq(1:length(active)), sep = ""))
  
  m.da.CM = cbind(seq(1:n.sites), m.da.CM)
  colnames(m.da.CM) = c("site", paste("da", seq(1:length(active)), sep = ""))
  
  # get rid of heme
  if (family == "globins") {
    m.da.CA = m.da.CA[1:n.aa, ]
    m.da.CM = m.da.CM[1:n.aa, ]
  }
  
  # calculate the minimum distance of each site to the active sites
  m.min.da.CA = matrix(nrow = n.aa, ncol = 1)
  m.min.da.CM = matrix(nrow = n.aa, ncol = 1)
  
  for (i in (1:n.aa)) {
    m.min.da.CA[i, ] = min(m.da.CA[i, 2:ncol(m.da.CA)])
    m.min.da.CM[i, ] = min(m.da.CM[i, 2:ncol(m.da.CM)])
  }
  
  # calculate a measure of all distances to the active site
  m.inv.da.CA = 1/m.da.CA[, 2:ncol(m.da.CA)]
  m.inv.da.CA[is.infinite(m.inv.da.CA)] <- NA
  if (length(active) > 1) {
    m.sum.inv.da.CA = rowSums(m.inv.da.CA, na.rm = T)
  } else {
    m.sum.inv.da.CA = m.inv.da.CA
  }
  
  m.inv.da.CM = 1/m.da.CM[, 2:ncol(m.da.CM)]
  m.inv.da.CM[is.infinite(m.inv.da.CM)] <- NA
  if (length(active) > 1) {
    m.sum.inv.da.CM = rowSums(m.inv.da.CM, na.rm = T)
  } else {
    m.sum.inv.da.CM = m.inv.da.CM
  }
  
  # build a dataframe with minimum distances
  min.da.data = data.frame("site" = seq(1:n.aa), 
                      "min.da.CM" = as.vector(m.min.da.CM), 
                      "min.da.CA" = as.vector(m.min.da.CA))
  
  # build a dataframe with sum.inv distances
  sum.inv.da.data = data.frame("site" = seq(1:n.aa), 
                      "sum.inv.da.CM" = as.vector(m.sum.inv.da.CM), 
                      "sum.inv.da.CA" = as.vector(m.sum.inv.da.CA))
  
  # write output files
  write.csv(min.da.data, file.path(out.dir, paste(p.ref, "_min.da.CM.ca.csv", sep = "")))
  write.csv(sum.inv.da.data, file.path(out.dir, paste(p.ref, "_sum.inv.da.CM.ca.csv", sep = "")))
  write.csv(m.da.CA, file.path(out.dir, paste(p.ref, "_m.da.ca.csv", sep = "")))
  write.csv(m.da.CM, file.path(out.dir, paste(p.ref, "_m.da.CM.csv", sep = "")))
}
  