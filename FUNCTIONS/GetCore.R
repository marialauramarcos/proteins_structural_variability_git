# Description:
# This function calculates the core of a multiple alignment. The core is composed of sites with no gaps (core)
# or of sites with no gaps and with 3 consecutive sites with no cagps (core.2). 
#
# Usage:
# AnalyzeAlignment(family, data.dir, p.ref)
#
#  Args:
#   - family: the name of the family. It can be "globins", "serinProteases", "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#   - data.dir: directory of the data. It must contain a file with the alignment of the family ("data.dir/family_alignment.txt") as obtanied
#   from Homstrad and a file with the dataset ("data.dir/family_dataset.csv"), containing in rows pdbids and chains in the multiple
#   pdb file obtained from Homstrad.
#   - p.ref: pdbid of the reference protein.
#
#  Required functions:
#   ReadFasta()


GetCore = function(family,
                   data.dir,
                   p.ref) {
  
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  
  alignment.fname <- file.path(data.dir, paste(family, "_alignment.txt", sep = ""))  
  
  # read the dataset
  dataset <- read.csv(dataset.fname)
  pdbid.dataset <- c(p.ref, as.character(dataset$pdbid))
  n.prot = length(pdbid.dataset)
  
  # read the alignment
  alignment.id <- ReadFasta(alignment.fname)
  pre.alignment <- alignment.id$ali[, -ncol(alignment.id$ali)]  # the last column is "*"
  l.alignment = ncol(pre.alignment)
  pdbid.alignment <- alignment.id$id
  
  # sort alignment
  alignment = matrix(nrow = n.prot, ncol = l.alignment)
  
  for (p in (1:n.prot)) {
    pdbid.p = pdbid.dataset[p]
    alignment[p, ] = pre.alignment[which(grepl(pdbid.p, pdbid.alignment)), ]
  }
  
  # Calculate core indexes
  gaps = matrix(nrow = l.alignment, ncol = 2)
  gaps[, 2] = seq(1:l.alignment)
  for (i in (1:l.alignment)) {
    ngaps = 0
    for (j in (1:n.prot)) {
      is.gap = is.gap(alignment[j, i])
      if (is.gap == "TRUE") {
        ngaps = ngaps + 1
      }
    }
    gaps[i, 1] = ngaps
  }
  core.index = gaps[gaps[, 1] == 0, 2]

  core.index.2 = c()
  for (i in (core.index)) {
    if ((i >= 3) & (i <= (l.alignment - 3))) {
      if (sum(gaps[(i - 3):(i + 3), 1]) == 0) {
        core.index.2 = c(core.index.2, i)
      }
    }
  }
  
  # Get aa indexes in the alignment
  p.index = matrix(nrow = n.prot, ncol = l.alignment)

  for (p in (1:n.prot)) {
    aa = 0
    for (i in (1:l.alignment)) {
      if (alignment[p, i] != "-") {
        aa = aa + 1
        p.index[p, i] = aa
      }
    }
  }

  # Get aa indexes in cores
  indexes.core = p.index[, core.index]
  indexes.core.2 = p.index[, core.index.2]

  write.csv(indexes.core, file.path(data.dir, paste("indexes.core", "_", family, ".csv", sep = ""))) 
  write.csv(indexes.core.2, file.path(data.dir, paste("indexes.core.2", "_", family, ".csv", sep = ""))) 
}

  
  