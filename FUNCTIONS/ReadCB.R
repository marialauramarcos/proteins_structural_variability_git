# This function reads a pdb file and returns CB equilibrium coordinates, sites and n.betas.
#
#  Args:
#    pdb.fname: pdb file name.
#    chain: chain to read.
#
#  Returns:
#    xyz.beta: Equilibrium coordinates of CBs.
#    site.beta: A vector with the sites.
#    n.betas: length(site.beta).
ReadCB <- function(pdb.fname, chain) {
  pdb <- read.pdb(file = pdb.fname)     
  sel <- atom.select(pdb, chain = chain, elety = "CB")    
  site.beta <- as.numeric(pdb$atom[sel$atom, c("resno")])    
  n.betas <- length(site.betas)    
  xyz.cbeta <- matrix(pdb$xyz[sel$xyz], ncol = n.betas, nrow = 3, byrow = F)    
  output <- list("xyz.cbeta" = xyz.cbeta, "site.beta" = site.beta, "n.betas" = n.betas)
  output
}

for (i in (1:n.aa)) { 
  sel <- atom.select(pdb, chain = chain, string = "protein", resno = i)    
  site.elety <- pdb$atom[sel$atom, c("elety")]     
  n.atoms = length(site.elety) 
  site.resid <- as.numeric(pdb$atom[sel$atom, c("resno")]) 
  xyz.site <- matrix(pdb$xyz[sel$xyz], nrow = 3, byrow = F)     
  
  weights = matrix(nrow = n.atoms)

  index.C = which(grepl("C", site.elety))
  weights[index.C, ] = 12

  for(j in (1:n.atoms)) {
    if (site.elety[j] == "CA") weights[j, ] = 0
    if (site.elety[j] == "N") weights[j, ] = 14
    if (site.elety[j] == "O") weights[j, ] = 16
    if (site.elety[j] == "S") weights[j, ] = 32 
  }

