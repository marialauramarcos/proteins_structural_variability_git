# This function reads a pdb file and returns coordinates of NA, NB, NC, ND and Fe of the Heme group.
#
#  Args:
#    pdb.fname: pdb file name.
#    chain: chain to read.
#
#  Returns:
#    xyz.heme: Equilibrium coordinates of NA, NB, NC, ND and Fe of the Heme group.
ReadHeme <- function(pdb.fname, chain) {
  pdb <- read.pdb(file = pdb.fname, ATOM.only = TRUE)     
  selNA <- atom.select(pdb, chain = chain, elety = "NA") 
  selNB <- atom.select(pdb, chain = chain, elety = "NB")    
  selNC <- atom.select(pdb, chain = chain, elety = "NC")    
  selND <- atom.select(pdb, chain = chain, elety = "ND")    
  selFE <- atom.select(pdb, chain = chain, elety = "FE")    
  xyz.heme <- matrix(c(pdb$xyz[selNA$xyz], 
                       pdb$xyz[selNB$xyz], 
                       pdb$xyz[selNC$xyz], 
                       pdb$xyz[selND$xyz], 
                       pdb$xyz[selFE$xyz]), nrow=3)
  xyz.heme
}
