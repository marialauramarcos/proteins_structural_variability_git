# Description:
#
# This function reads a pdb file and returns the coordinates of the center of mass of side chains of each aminoacid.
#
#  Args:
#   - pdb.fname: pdb filename.
#   - chain: chain to read.
#
#  Returns:
#    xyz.side.chain.CM: coordinates of the center of mass of side chains of each aminoacid.

CalculateSideChainCM <- function(pdb.fname, chain) {
  
  # read pdb
  pdb <- read.pdb(file = pdb.fname)
  
  # get information from de pdb
  sel.ca <- atom.select(pdb, chain = chain, elety = "CA")
  site.ca = as.numeric(pdb$atom[sel.ca$atom, c("resno")])
  n.aa = length(site.ca)  
  
  # create a matrix to the save the centers of mass
  xyz.side.chain.CM = matrix(nrow = 3, ncol = n.aa)
  
  # eliminate duplicated values
  site.ca = site.ca[!duplicated(site.ca)] 
  
  # start to count sites
  site = 0
  
  # start a loop to analyze each residue
  for (i in (site.ca)) { 
    
    # set duplicate
    duplicate = "FALSE"
    
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
    
    # eliminate the weights of CAs
    for(j in (1:n.atoms)) {
      if (site.elety[j] == "CA") weights[j, ] = 0
    }
    
    # backbone
    bb = c(1, 2, 3, 4)
    
    # calculate the center of mass of well enumerated residues      
    xyz.side.chain.CM[1, site] = sum(xyz.site[1, ][-bb]  * weights[-bb] ) / sum(weights[-bb] )
    xyz.side.chain.CM[2, site] = sum(xyz.site[2, ][-bb]  * weights[-bb] ) / sum(weights[-bb] )
    xyz.side.chain.CM[3, site] = sum(xyz.site[3, ][-bb]  * weights[-bb] ) / sum(weights[-bb] )
      
    # correct for duplicated indexes in pdb
    if (length(which(site.elety == "CA")) > 1) {
      
      # set suplicate
      duplicate = "TRUE"
          
      # print a warning message
      print("warning: duplicated indexes in the pdb file")
          
      c.j = c((which(site.elety == "CA") - 1), (n.atoms + 1))
    }
    
    # operate with duplicates
    if (duplicate == "TRUE") {
      for (j in (2:length(c.j))) {
         
        # count a site
        if (j != 2) {site = site + 1}
        
        # set interval for j
        int = (c.j[j - 1]:(c.j[j] - 1)) # c.j[j] is the biginning of the interval, so the previous interval finishes in c.j[j] - 1 
        
        # calculate the center of mass of the interval
        xyz.side.chain.CM[1, site] = sum(xyz.site[1, int][-bb] * weights[int, ][-bb]) / sum(weights[int, ][-bb])
        xyz.side.chain.CM[2, site] = sum(xyz.site[2, int][-bb] * weights[int, ][-bb]) / sum(weights[int, ][-bb])
        xyz.side.chain.CM[3, site] = sum(xyz.site[3, int][-bb] * weights[int, ][-bb]) / sum(weights[int, ][-bb])
      }
    }
  }
  # build a side chain for sites with no information (exept for first and last residue)
  for (i in (2:(ncol(xyz.side.chain.CM) - 1))) {
    if (xyz.side.chain.CM[1, i] == "NaN") {
      xyz.side.chain.CM[, i] = (xyz.side.chain.CM[, (i - 1)] + xyz.side.chain.CM[, (i + 1)])/2
    }
  }
  xyz.side.chain.CM
}

  

