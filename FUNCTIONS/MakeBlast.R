# input
 n.seq = 10
 p.ref = "1a6m"
 chain = "A"

# libraries
library("bio3d")

# get the pdb file of p.ref
get.pdb(p.ref)

# read the file
pdb <- read.pdb(paste(p.ref, ".pdb", sep = "")) 
sel <- atom.select(pdb, chain = chain, elety = "CA")  
seq <- pdb$atom[sel$atom, c("resid")]

# code maps
code3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
           "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
           "Tyr", "Val")
code1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
           "M", "F", "P", "S", "T", "W", "Y", "V")

# for each code replace 3letter code by 1letter code:
for (i in 1:length(code3)) {
  seq <- gsub(code3[i], code1[i], seq, ignore.case = TRUE)
}

# make a blast
blast.info <- blast.pdb(seq, database = "pdb")
score = blast.info$hit.tbl$bitscore
pdbid = blast.info$hit.tbl$pdb.id
identity = blast.info$hit.tbl$identity

# filter too similar proteins
filter = (pdbid != paste(toupper(p.ref), "_", toupper(chain), sep = ""))
pdbid = pdbid[filter]
identity = identity[filter]

filter2 = (identity < 90)
pdbid = pdbid[filter2]
identity = identity[filter2]
