# Description:
#
# This function reads a pdb file and returns the coordinates of the center of mass of side chains of each aminoacid.
#
#  Args:
#   - pdb.fname: pdb filename.
#   - chain: chain to read.
#
#  Returns:
#    xyz.side.chain.CM: coordinates of the centers of mass of side chains of each aminoacid.

CalculateSideChainCM <- function(pdb.fname, chain) {
  
  # Read pdb
  pdb <- read.pdb(file = pdb.fname)
  
  # get information from de pdb
  sel.ca <- atom.select(pdb, chain = chain, elety = "CA")
  site.ca = as.numeric(pdb$atom[sel.ca$atom, c("resno")])    
  n.aa = length(site.ca)  
  
  # create a matrix to the save the centers of mass
  xyz.side.chain.CM = matrix(nrow = 3, ncol = n.aa)
  
  # start to count sites
  site = 0
  
  # key for run
  key = 1
  # star a loop to analyze each residue
  for (i in (site.ca)) { 
    if (key == 1) {
      
      # count a site
      site = site + 1
      print(site)
  
      # Extract information of the resid
      sel <- atom.select(pdb, chain = chain, strain = "protein", resno = i)    
      site.elety <- pdb$atom[sel$atom, c("elety")] 
      site.resid <- pdb$atom[sel$atom, c("resid")]     
      n.atoms = length(site.elety) 
      xyz.site <- matrix(pdb$xyz[sel$xyz], nrow = 3, byrow = F)     
  
      # Get weights
      weights = matrix(nrow = n.atoms)
  
      index.C = which(grepl("C", site.elety))
      index.O = which(grepl("O", site.elety))
      index.N = which(grepl("N", site.elety))
      index.S = which(grepl("S", site.elety))
    
      weights[index.C, ] = 12
      weights[index.O, ] = 16
      weights[index.N, ] = 14
      weights[index.S, ] = 32
    
      # Eliminate the weights of CAs
      for(j in (1:n.atoms)) {
        if (site.elety[j] == "CA") weights[j, ] = 0
      }
    
      # calculate the center of mass of well enumerated residues
      xyz.side.chain.CM[1, site] = sum(xyz.site[1, ] * weights) / sum(weights)
      xyz.side.chain.CM[2, site] = sum(xyz.site[2, ] * weights) / sum(weights)
      xyz.side.chain.CM[3, site] = sum(xyz.site[3, ] * weights) / sum(weights)
      
      # correct for bad indexes in pdf
      for (j in (2:(n.atoms - 1))) {
        if(site.resid[j - 1] != site.resid[j]) {
        
          print(j)
          # print a warning message
          print("warning: duplicated indexes in the pdb file")
        
          # Calculate the center of mass of the first residue
          xyz.side.chain.CM[1, site] = sum(xyz.site[1, 1:(j - 1)] * weights[1:(j - 1), ]) / sum(weights[1:(j - 1), ])
          xyz.side.chain.CM[2, site] = sum(xyz.site[2, 1:(j - 1)] * weights[1:(j - 1), ]) / sum(weights[1:(j - 1), ])
          xyz.side.chain.CM[3, site] = sum(xyz.site[3, 1:(j - 1)] * weights[1:(j - 1), ]) / sum(weights[1:(j - 1), ])
        
          site = site + 1
          print(site)
        
          # calculate the center of mass of the second residue
          xyz.side.chain.CM[1, site] = sum(xyz.site[1, j:n.atoms] * weights[j:n.atoms, ]) / sum(weights[j:n.atoms, ])
          xyz.side.chain.CM[2, site] = sum(xyz.site[2, j:n.atoms] * weights[j:n.atoms, ]) / sum(weights[j:n.atoms, ])
          xyz.side.chain.CM[3, site] = sum(xyz.site[3, j:n.atoms] * weights[j:n.atoms, ]) / sum(weights[j:n.atoms, ])
        
          key = 0
        }
      }
    }else{
    key = 1
    }
  }
  xyz.side.chain.CM
}

  

