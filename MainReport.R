# Description: This program generates reports with the output of "MainProgram.R" and "MainProgramCM.R".
# To run the program it is neccessary to fill the input file "input_MainReport.csv".
# Files needed for each family in data.dir are:
#  "family_list.txt": with all of the proteins of the family (including p.ref).
#  "family_ref.txt": p.ref.
#  "p.ref.pdb": pdb file of p.ref.
#  "p.ref_consurf.csv": sequence divergence scores obtained from ConsurfDB.
# Files needed for each family in "DATA" are:
#  "p.ref_functionalSites.csv": functional sites of p.ref.
#  "p.ref_dssp.csv": ASA of sites of p.ref calculated by dssp.
#  "surfase_area.csv": surfase area of aa to calculate RSA.

### PROGRAM ###

# load packages
library(knitr)
library(markdown)

# read input
input.fname <- "input_MainReport.csv"
input <- read.csv(input.fname)

# satart a loop for each family
for (f in (1:nrow(input))) { 
  print(f)
  family <- as.character(input$family)[f]
  type <- as.character(input$type)[f]
  p.ref <- as.character(input$p.ref)[f]
  enm <- as.character(input$enm)[f]
  n.mut.p <- input$n.mut.p[f]
  R0.CA = input$R0.CA[f]
  R0.CM = input$R0.CM[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]
  K.analysis <- input$K.analysis[f]
  
  print(family)
  
  ## ANM CA
  data.dir <- "OUT/out_subset_CA_ANM"
  out.dir <- "FIGURES_REPORTS"
  
  R0 = R0.CA
  
  ### Structure - Normal modes analysis
  rmarkdown::render('analysis-normal-modes.Rmd', 
                    output_file =  paste(out.dir, "/report_normal_modes_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### Structure - Cartesian coordinates analysis
  rmarkdown::render('analysis-structure.Rmd', 
                    output_file =  paste(out.dir, "/report_structure_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### 1/LPD vs SD(EXP) - Comparison to determine the best R0 
  rmarkdown::render('analysis-structure-CN-WCN-exp.Rmd', 
                    output_file =  paste(out.dir, "/report_comparison_LPD_SDexp_", family, ".html", sep = ''))

  ### MSF
  #rmarkdown::render('analysis-dynamical-MSF-structure.Rmd', 
  #                  output_file =  paste(out.dir, "/report_dynamical_MSF_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nH
  #rmarkdown::render('analysis-dynamical-nH.Rmd', 
  #                  output_file =  paste(out.dir, "/report_dynamical_nH_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nR
  #rmarkdown::render('analysis-dynamical-nR.Rmd', 
  #                  output_file =  paste(out.dir, "/report_dynamical_nR_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD core
  #rmarkdown::render('analysis-structure-core.Rmd', 
  #                 output_file =  paste(out.dir, "/report_structure_CA_core_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD windows contacts rot
  #rmarkdown::render('analysis-structure-window-contacts.Rmd', 
  #                 output_file =  paste(out.dir, "/report_structure_window_contacts_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### local enviroment 
  #rmarkdown::render('analysis-structure-local-enviroment.Rmd', 
  #                 output_file =  paste(out.dir, "/report_structure_local_enviroment_CA_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
 
  ## CM
  #data.dir <- paste("OUT/out_subset_CM_ANM", sep = "")
  #R0 = R0.CM
  
  ### Pn
  #rmarkdown::render('analysis-normal-modes.Rmd', 
  #                   output_file =  paste(out.dir, "/report_normal_modes_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD
  #rmarkdown::render('analysis-structure.Rmd', 
  #                  output_file =  paste(out.dir, "/report_structure_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### MSF
  #rmarkdown::render('analysis-dynamical-MSF-structure.Rmd', 
  #                    output_file =  paste(out.dir, "/report_dynamical_MSF_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nH
  #rmarkdown::render('analysis-dynamical-nH.Rmd', 
  #                 output_file =  paste(out.dir, "/report_dynamical_nH_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### nR
  #rmarkdown::render('analysis-dynamical-nR.Rmd', 
  #                  output_file =  paste(out.dir, "/report_dynamical_nR_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD core
  #rmarkdown::render('analysis-structure-core.Rmd', 
  #                  output_file =  paste(out.dir, "/report_structure_CM_core_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### RMSD windows contacts rot
  #rmarkdown::render('analysis-structure-window-contacts.Rmd', 
  #                  output_file =  paste(out.dir, "/report_structure_window_contacts_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ### local enviroment 
  #rmarkdown::render('analysis-structure-local-enviroment.Rmd', 
  #                  output_file =  paste(out.dir, "/report_structure_local_enviroment_CM_", family, "_", enm, "_R0_", R0, ".html", sep = ''))
  
  ## pfANM CA
  #R0 = "NC"
  #data.dir <- paste("OUT/out_subset_CA_pfANM", sep = "")
  
  ### RMSD
  #rmarkdown::render('analysis-structure-pfANM.Rmd', 
  #                  output_file =  paste(out.dir, "/report_structure_pfANM_", family, ".html", sep = ''))
}
