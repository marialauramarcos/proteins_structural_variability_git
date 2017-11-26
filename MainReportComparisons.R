# Description: This program generates comparative reports with the output of "MainProgram.R", "MainProgramCM.R" and "MainReport".
# To run the program it is neccessary to fill the input file "input_MainReport.csv".

### PROGRAM ###

# load packages
library(knitr)
library(markdown)

# read input
input.fname <- "input_MainReport.csv"
input <- read.csv(input.fname)

## all models

### 1/LPD vs SD(EXP) - Comparison to determine the best R0 - all families
rmarkdown::render('comparison_families_LPD_vs_SDexp.Rmd', 
                  output_file =  paste("FIGURES_REPORTS/report_comparison_families_LPD_vs_SDexp", ".html", sep = ""))
## ANM CA
data.dir <- "OUT/out_subset_CA_ANM"
enm = "ANM"
R0 = 12.5 # it might be changed

### comparisons figures structure all families
rmarkdown::render('figures2_structure_all_families.Rmd', 
                  output_file =  paste("FIGURES_REPORTS/report_comparison_figures_structure_all_families_CA","_", enm, "_R0_", R0, ".html", sep = ''))

## ANM CM
#data.dir <- "OUT/out_subset_CM_ANM"
#enm = "ANM"
#R0 = 10 # it might be changed

### comparisons figures structure all families
#rmarkdown::render('figures2_structure_all_families.Rmd', 
#                  output_file =  paste("FIGURES_REPORTS/report_comparison_figures_structure_all_families_CM","_", enm, "_R0_", R0, ".html", sep = ''))
## pfANM CA
#data.dir <- "OUT/out_subset_CA_pfANM"
#enm = "pfANM"
#R0 = NC #

### comparisons figures structure all families
#rmarkdown::render('figures2_structure_all_families.Rmd', 
#                  output_file =  paste("FIGURES_REPORTS/report_comparison_figures_structure_all_families_CA", "_", enm, "_R0_", R0, ".html", sep = ''))

## pfANM CM
#data.dir <- "OUT/out_subset_CM_pfANM"
#enm = "pfANM"
#R0 = NC #

### comparisons figures structure all families
#rmarkdown::render('figures2_structure_all_families.Rmd', 
#                  output_file =  paste("FIGURES_REPORTS/report_comparison_figures_structure_all_families_CM", "_", enm, "_R0_", R0, ".html", sep = ''))
