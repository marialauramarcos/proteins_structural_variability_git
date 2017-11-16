# Description:
#
# This is the main program of the project. The program simulates multiple mutants of a given protein using the "Linearly Forced - 
# Elastic Network Model" (LF-ENM) with different selection regimens according to the "Stress Model". The program also
# analyzes the multiple alignment of the family to which the protein belongs and calculates measures 
# of variabilty of theoretical and experimental structures.
# Experimental structures and the multiple alignment must be obtained from Homstrad.
#
# To run the program it is necessary to previously fill the input ("input_MainProgram.csv") with the following information:
#
#    - family: The family of the protein to mutate.
#    - p.ref: The pdb code (pdbid) of the protein to mutate (example: "1a6m"). The protein must be a member of
#    the selected family. This pdbid must not be included in the dataset ("DATA/family_dataset.csv").
#    - chain.p.ref: The chain of p.ref in the pdb file obtained from Homstrad.
#    - n.mut.p: The number of mutants to generate for each member of the family. For example, if the family has 20 
#    members, the program generates n.mut.p x 20 mutants.
#    - fmax: Argument for the "LF-ENM". It is the maximun value for the forces that model the mutations. 
#    - R0.CA: The Cut-off for the "ANM" (Anisotropic Network Model) that represents the proteins for the CA model.
#    - R0.CM: The Cut-off for the "ANM" (Anisotropic Network Model) that represents the proteins for the two nodes per site model. Not used in this program.
#    - rotate: It can be "TRUE" or "FALSE". If it is "TRUE", r.p.2 is rotated in order to minimize RMSD with r.p.ref.
#    - heme: Argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    - calculate.betas: It can be "TRUE" or "FALSE". If it is "TRUE", the program calculates betas of the "Stress Model".
#    - analyze.family: It can be "TRUE" or "FALSE". If it is "TRUE", the program analyzes the family.
#    - generate.mutants: It can be "TRUE" or "FALSE". If it is "TRUE", the program generates new mutants.
#    - analyze.experimental.theoretical: It can be "TRUE" or "FALSE". If it is "TRUE", the program analyzes the proteins.
#    - K.analysis: It can be "K" or "Keff". For "K" or "Keff", the analysis is based on normal modes of "K" or "Keff"
#    respectibly.

### PROGRAM ###

# remove objects from the workspace
rm(list = ls())

# load packages
library(bio3d) 
library(seqinr) 

# data dir
data.dir <- "DATA"

# set Elastic Network Model: "ANM" or "pfANM"
model <- "ANM"

# output dir
if (model == "ANM") out.dir <- "OUT/out_subset_CA_ANM"
if (model == "pfANM") out.dir <- "OUT/out_subset_CA_pfANM"

# general parameters
tolerance = 1e-10

# function filenames
AnalyzeExperimentalTheoretical.fname <- "FUNCTIONS/AnalyzeExperimentalTheoretical.R"
AnalyzeFamily.fname <- "FUNCTIONS/AnalyzeFamily.R"
AnalyzeAlignment.fname <- "FUNCTIONS/AnalyzeAlignment.R"
GenerateMutants.fname <- "FUNCTIONS/GenerateMutants.R"
ReadFasta.fname <- "FUNCTIONS/ReadFasta.R"
ReadCA.fname <- "FUNCTIONS/ReadCA.R"
ReadHeme.fname <- "FUNCTIONS/ReadHeme.R"
CalculateSideChainCM.fname <- "FUNCTIONS/CalculateSideChainCM.R"
CalculateENMKeff.fname <- "FUNCTIONS/CalculateENMKeff.R"
CalculateENMK.fname <- "FUNCTIONS/CalculateENMK.R"
CalculateVariability.fname <- "FUNCTIONS/CalculateVariability.R"
CalculateDynamicalVariability.fname <- "FUNCTIONS/CalculateDynamicalVariability.R"
CalculateDaCMCA.fname <- "FUNCTIONS/CalculateDaCMCA.R"
GetCore.fname <- "FUNCTIONS/GetCore.R"
WindowsRMSD.fname <- "FUNCTIONS/WindowsRMSD.R"
WindowsRMSDcontacts.fname <- "FUNCTIONS/WindowsRMSDcontacts.R"
if (model == "ANM") {
  CalculateBetas.fname <- "FUNCTIONS/CalculateBetas.R"
  CalculateKij.fname <- "FUNCTIONS/CalculateKij.R"
  CalculateForce.fname <- "FUNCTIONS/CalculateForce.R"
}
if (model == "pfANM") {
  CalculateBetas.fname <- "FUNCTIONS/CalculateBetasPFANM.R"
  CalculateKij.fname <- "FUNCTIONS/CalculateKijPFANM.R"
  CalculateForce.fname <- "FUNCTIONS/CalculateForcePFANM.R"
}

# source functions
source(AnalyzeExperimentalTheoretical.fname)
source(AnalyzeFamily.fname)
source(AnalyzeAlignment.fname)
source(GenerateMutants.fname)
source(ReadFasta.fname) 
source(ReadCA.fname) 
source(ReadHeme.fname)
source(CalculateSideChainCM.fname)
source(CalculateENMKeff.fname)
source(CalculateENMK.fname)
source(CalculateVariability.fname)
source(CalculateDynamicalVariability.fname)
source(CalculateDaCMCA.fname)
source(GetCore.fname)
source(WindowsRMSD.fname)
source(WindowsRMSDcontacts.fname)
source(CalculateBetas.fname)
source(CalculateKij.fname)
source(CalculateForce.fname)

# read input
input.fname <- file.path("input_MainProgram.csv")
input <- read.csv(input.fname)

# start a loop to analyze each family
for (f in (1:nrow(input))) { 
  family <- as.character(input$family)[f]
  p.ref <- as.character(input$p.ref)[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]
  n.mut.p = input$n.mut.p[f]
  fmax = input$fmax[f] 
  R0 = input$R0.CA[f]
  rotate <- input$rotate[f]
  heme <- input$heme[f]
  calculate.betas <- input$calculate.betas[f]
  analyze.family <- input$analyze.family[f]
  generate.mutants <- input$generate.mutants[f]
  analyze.experimental.theoretical <- input$analyze.experimental.theoretical[f]
  K.analysis <- input$K.analysis[f]
  
  print(family)
  
  # analyze the alignment of the family
  if (analyze.family == "TRUE") {
    AnalyzeFamily(family,
                  p.ref, 
                  data.dir,
                  out.dir)
  }
  
  # get the core of the alignment
  GetCore(family,
          data.dir,
          p.ref)
  
  # calculate distances to the active site
  CalculateDaCMCA(family,
                  p.ref,
                  chain.p.ref,
                  heme,
                  data.dir,
                  out.dir)
  
  # generate id for betas output filename
  betas.fname.id <- paste(family, "_", p.ref, "_R0_", R0, sep = "")
  
  # calculate betas of the "Stress Model"
  if (calculate.betas == "TRUE") {
    CalculateBetas(chain.p.ref,
                   fmax, 
                   R0,
                   heme,
                   data.dir,
                   out.dir,
                   betas.fname.id,
                   tolerance)
  }
  
  # read betas and stablish selection regimens
  all.betas <- read.csv(file.path(out.dir, paste(betas.fname.id, "_out_all.betas.csv", sep = "")))
  regimens <- c("strong.sel", "medium.sel", "weak.sel", "no.sel")
  
  # start a loop for each beta
  for (b in all.betas)  {
    
    # filter regimens
    if (regimens[all.betas == b] != "weak.sel") {
      if (regimens[all.betas == b] != "strong.sel") {
    
        # generate ids for output filenames
        mut.fname.id <- paste(family, "_R0_", R0, "_beta_", regimens[all.betas == b], sep = "")
        analysis.fname.id <- paste(mut.fname.id, "_K.analysis_", K.analysis, sep = "")
  
        # generate mutants
        if (generate.mutants == "TRUE") {
          GenerateMutants(family,
                          chain.p.ref, 
                          n.mut.p,
                          fmax, 
                          R0,
                          b,
                          heme, 
                          data.dir,
                          out.dir,
                          mut.fname.id,
                          tolerance)
        }
  
        # calculate measures of variability of theoretical and experimental proteins
        if (analyze.experimental.theoretical == "TRUE") {
          AnalyzeExperimentalTheoretical(family,
                                         p.ref,
                                         chain.p.ref,
                                         n.mut.p,
                                         R0, 
                                         rotate,
                                         heme,
                                         K.analysis,
                                         data.dir,
                                         out.dir,
                                         mut.fname.id, 
                                         analysis.fname.id,
                                         tolerance)
        }
      }
    }
  }
}