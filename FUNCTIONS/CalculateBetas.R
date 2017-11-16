# Description:
#
# This function caluclates 4 different values of the parameter "beta" for no selection, weak selection,
# medium selection and stong selection following the "Stress Model". 
#
# Usage:
#
# CalculateBetas(chain.p.ref, fmax, R0, heme = TRUE/FALSE, data.dir, out.dir,
# betas.fname.id, tolerance)
#
#  Args:
#    - chain.p.ref: The chain of p.ref in the pdb file obtained from Homstrad.
#    - fmax: The maximun value for the forces that model the mutations.
#    - R0: The cut-off for the ANM ("Anisotropic Network Model") that represents the proteins.
#    - heme: Argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    - data.dir: Directory of the data. It must contain the pdb file obtained from Homstrad ("data.dir/family_coordinates.csv").
#    - out.dir: Output directory.
#    - betas.fname.id: ID for output filenames.
#    - tolerance: 0 tolerance.
#
#  Required libraries:
#    {Bio3d}
#
#  Required functions:
#    ReadCA()
#    ReadHeme()
#    CalculateENMK()
#    CalculateKij()

CalculateBetas <- function(chain.p.ref,
                           fmax, 
                           R0,
                           heme,
                           data.dir,
                           out.dir,
                           betas.fname.id,
                           tolerance) {
  
  print("calculating betas...")
  
  ### READ EXPERIMENTAL DATA ###
  
  # filenames
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  
  # read PDB of p.ref
  pdb = ReadCA(pdbs.fname, chain.p.ref)
  r.p.ref = pdb$xyz.calpha
  n.aa = pdb$n.sites
  
  # get heme coordinates, add them to CA´s coordinates and calculate the new number of sites
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdbs.fname, chain.p.ref)
    r.p.ref = cbind(r.p.ref, r.heme)
    n.sites = ncol(r.p.ref)
  } else {
    n.sites = n.aa
  }
  
  ### CALCULATE ENM OF THE REFERENCE PROTEIN ###
  
  # calculate ENM of p.ref
  ENMK.p.ref = CalculateENMK(r.p.ref, CalculateKij, R0, tolerance) 
  
  ### CALCULATE PARAMETERS FOR THE SELECTION OF MUTANTS FOLLOWING THE "STRESS MODEL" ###
  
  # calculate <CN>
  mean.CN = mean(colSums(ENMK.p.ref$kij)) 
  
  # calculate "f2" (see the function CalculateForce())
  f2 = mean((fmax * runif(1e6, -1, 1)) ^ 2)
  
  # caculate betas for differente selection regimes
  beta.0.1 = round(- log(0.1) / (f2 * mean.CN), digits = 3) # Strong selection
  beta.0.5 = round(- log(0.5) / (f2 * mean.CN), digits = 3) # Medium selection
  beta.0.9 = round(- log(0.9) / (f2 * mean.CN), digits = 3) # Weak selection
  beta.1 = round(- log(1) / (f2 * mean.CN), digits = 3) # No selection
  
  # build a list with betas
  all.betas = list("strong.sel" = beta.0.1, "medium.sel" = beta.0.5, "weak.sel" = beta.0.9, "no.sel" = beta.1) 
  
  # create a file to save the data
  write.csv(all.betas, file = file.path(out.dir, paste(betas.fname.id, "_out_all.betas.csv", sep = "")), row.names = FALSE)
}
