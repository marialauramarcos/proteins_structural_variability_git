# Description: This program generates comparative reports with the output of "MainProgram.R" and "MainProgramCM.R".
# To run the program it is neccessary to fill the input file "input_MainReport.csv".
# Files needed for each family in data.dir are:
#  "family_list.txt": with all of the proteins of the family (including p.ref).
#  "family_ref.txt": p.ref.
#  "p.ref.pdb": pdb file of p.ref.
#  "p.ref_consurf.csv": sequence divergence scores obtained from ConsurfDB.

### PROGRAM ###

# load packages
library(knitr)
library(markdown)

# read input
input.fname <- "input_MainReport.csv"
input <- read.csv(input.fname)


## comparisons LPD vs SDexp all families
rmarkdown::render('comparison_families_LPD_vs_SDexp.Rmd', 
                  output_file =  paste("FIGURES_REPORTS/report_comparison_families_LPD_vs_SDexp", ".html", sep = ""))

## comparisons figures structure all families
R0 = 12.5
enm = "ANM"
data.dir = "OUT/out_subset_CA_ANM"

rmarkdown::render('figures2_structure_all_families.Rmd', 
                  output_file =  paste("FIGURES_REPORTS/report_comparison_figures_structure_all_families_CA","_", enm, "_R0_", R0, ".html", sep = ''))

## comparisons mutants all families
#rmarkdown::render('comparison_families.Rmd', 
#                  output_file =  paste("FIGURES_REPORTS/report_comparison_families", ".html", sep = ''))

# satart a loop for each family
for (f in (1:nrow(input))) { 
  print(f)
  family <- as.character(input$family)[f]
  type <- as.character(input$type)[f]
  p.ref <- as.character(input$p.ref)[f]
  enm <- as.character(input$enm)[f]
  n.mut.p <- input$n.mut.p[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]
  print(family)
  
  # generate reports - comparisons CA - CM and R0s
  
  ## comparisons CA CM
  
  #R0.CA = 12.5
  #R0.CM = 10
  #rmarkdown::render('comparison_CA_CM.Rmd', 
  #                  output_file =  paste("OUT/report_comparison_CA_CM_", family, "_", enm,"_R0_", R0.CA, "_", R0.CM, ".html", sep = ''))
  
  ## comparison R0s
  
  ### CA
  #data.dir <- paste("OUT/out_subset_CA_ANM", sep = "")
  
  #R0.1 = 10
  #R0.2 = 12.5
  
  #rmarkdown::render('comparison_R0.Rmd', 
  #                  output_file =  paste("OUT/report_comparison_CA_", family, "_", enm,"_R0_", R0.1, "_", R0.2, ".html", sep = ''))
  
  ### CM
  #data.dir <- paste("OUT/out_subset_CM_ANM", sep = "")
  
  #R0.1 = 7.5
  #R0.2 = 10
  
  #rmarkdown::render('comparison_R0.Rmd', 
  #                  output_file =  paste("OUT/report_comparison_CM_", family, "_", enm,"_R0_", R0.1, "_", R0.2, ".html", sep = ''))
  
}
