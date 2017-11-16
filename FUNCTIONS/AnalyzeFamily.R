# Description:
# This function analyzes a multiple alignment of a family of proteins comparing a reference protein "p.ref" with the other proteins. 
#
# Usage:
# AnalyzeFamily(family, p.ref, data.dir, out.dir)
#
#  Args:
#   - family: the name of the family. It can be "globins", "serinProteases", "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#   - p.ref:  the pdb code (pdbid) of the protein to mutate (example: "1a6ma"). The protein must be a member of
#   the selected family. This pdbid must not be included in the dataset ("DATA/family_dataset.csv").
#   - data.dir: directory of the data. It must contain a file with the alignment of the family ("data.dir/family_alignment.txt") as obtanied
#   from Homstrad and a file with the dataset ("data.dir/family_dataset.csv"), containing in rows pdbids and chains in the multiple
#   pdb file obtained from Homstrad.
#   - out.dir: output directory.
#
#  Required functions:
#   ReadFasta()
#   AnalyzeAlignment()
#
#  Returns:
#   Files with the result of the analysis in out.dir.

AnalyzeFamily <- function(family, 
                          p.ref, 
                          data.dir, 
                          out.dir) {
  
  print("analyzing families...")
  
  # data filenames
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  
  alignment.fname <- file.path(data.dir, paste(family, "_alignment.txt", sep = ""))  
  
  # read the dataset
  dataset <- read.csv(dataset.fname)
  pdbid.dataset <- as.character(dataset$pdbid) 
  n.prot = length(pdbid.dataset)
  
  # read the alignment
  alignment.id <- ReadFasta(alignment.fname)
  alignment <- alignment.id$ali[, -ncol(alignment.id$ali)]  # the last column is "*"
  pdbid.alignment <- alignment.id$id
  l.alignment = ncol(alignment)
  
  # create matrices to save data
  m.n.sites.p.ref = matrix(nrow = n.prot, ncol = 1)
  m.n.sites.p.2 = matrix(nrow = n.prot, ncol = 1)
  m.n.aligned = matrix(nrow = n.prot, ncol = 1)
  m.aligned.p.ref.index = matrix(nrow = n.prot, ncol = l.alignment)
  m.aligned.p.2.index = matrix(nrow = n.prot, ncol = l.alignment)
  m.not.aligned.p.ref.index = matrix(nrow = n.prot, ncol = l.alignment)
  m.not.aligned.p.2.index = matrix(nrow = n.prot, ncol = l.alignment)
  m.n.aligned.mut.p.ref = matrix(nrow = n.prot, ncol = 1)
  m.aligned.mut.p.ref.index = matrix(nrow = n.prot, ncol = l.alignment)
  m.identity = matrix(nrow = n.prot, ncol = 1)
  
  # start a loop to evaluate each protein of the family
  for (P in (1:n.prot)) {
    print(P)
    
    # get pdbid 
    p.2 <- pdbid.dataset[P]
    
    # anylize the alignment
    analysis.alignment = AnalyzeAlignment(alignment, 
                                          pdbid.alignment,   
                                          p.ref, 
                                          p.2)
    
    # Save the analysis in the matrices
    m.n.sites.p.ref[P, ] = analysis.alignment$n.sites.p.ref
    m.n.sites.p.2[P, ] = analysis.alignment$n.sites.p.2
    m.n.aligned[P, ] = analysis.alignment$n.aligned
    m.aligned.p.ref.index[P, 1:length(analysis.alignment$aligned.p.ref.index)] = analysis.alignment$aligned.p.ref.index
    m.aligned.p.2.index[P, 1:length(analysis.alignment$aligned.p.2.index)] = analysis.alignment$aligned.p.2.index
    if (length(analysis.alignment$not.aligned.p.ref.index) > 0) {
      m.not.aligned.p.ref.index[P, 1:length(analysis.alignment$not.aligned.p.ref.index)] = analysis.alignment$not.aligned.p.ref.index
    }
    if (length(analysis.alignment$not.aligned.p.2.index) > 0) {
      m.not.aligned.p.2.index[P, 1:length(analysis.alignment$not.aligned.p.2.index)] = analysis.alignment$not.aligned.p.2.index
    }
    m.n.aligned.mut.p.ref[P, ] = analysis.alignment$n.aligned.mut.p.ref
    if (analysis.alignment$n.aligned.mut.p.ref > 0) {
      m.aligned.mut.p.ref.index[P, 1:length(analysis.alignment$aligned.mut.p.ref.index)] = analysis.alignment$aligned.mut.p.ref.index
    }
    m.identity[P, ] = analysis.alignment$identity
  }
  
  # Create files to save the data
  write.csv(m.n.sites.p.ref, file = file.path(out.dir, paste(family, "_out_m.n.sites.p.ref.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.sites.p.2, file = file.path(out.dir, paste(family, "_out_m.n.sites.p.2.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.aligned, file = file.path(out.dir, paste(family, "_out_m.n.aligned.csv", sep = "")), row.names = FALSE)
  write.csv(m.aligned.p.ref.index, file = file.path(out.dir, paste(family, "_out_m.aligned.p.ref.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.aligned.p.2.index, file = file.path(out.dir, paste(family, "_out_m.aligned.p.2.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.not.aligned.p.ref.index, file = file.path(out.dir, paste(family, "_out_m.not.aligned.p.ref.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.not.aligned.p.2.index, file = file.path(out.dir, paste(family, "_out_m.not.aligned.p.2.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.aligned.mut.p.ref, file = file.path(out.dir, paste(family, "_out_m.n.aligned.mut.p.ref.csv", sep = "")), row.names = FALSE)
  write.csv(m.aligned.mut.p.ref.index, file = file.path(out.dir, paste(family, "_out_m.aligned.mut.p.ref.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.identity, file = file.path(out.dir, paste(family, "_out_m.identity.csv", sep = "")), row.names = FALSE)
}
