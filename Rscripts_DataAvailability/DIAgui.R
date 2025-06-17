################################################################################
######  DIAgui for DIA-NN output                                   #############
################################################################################
#DIAgui is an R package that contains a user-friendly shiny app to process output from DIA-nn 
# proteomics software. Since shiny can be quite slow to process big data file, a function that does 
# the complete workflow and saves all results in a directory is also present. The process is as follow: 
# filtering data according q.values, applying MaxLFQ algorithm for quantification, getting some other 
# informations, plot some graphs to check data quality.
# This package is based on diann-rpackage from Vadim Demichev. The c++ source code is exactly the same as 
# the one from V.Demichev, it is used for MaxLFQ algorithm. You can also use iq R package which is way faster 
# to run MaxLFQ algorithm.
# In addition to filter data and calculate MaxLFQ intensities, you can also get the Top3 and iBAQ quantification
# within the app.

# Install DIAgui from github:
if(!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools") 
}
devtools::install_github("mgerault/DIAgui")

# You can now load it and run the app with this commands:
library(DIAgui)
runDIAgui()

# Toy dataset
# If you want to test the application, you can download a toy dataset here which contain a 
# small report from saccharomyces cerevisiae with its corresponding FASTA file. For information 
# about the data, you access access its documentation in R via ?small_report

# Process data from DIAnn in a single workflow, save it and also save some quality check plots.
x <- diann_load("C:/Users/Laura/Downloads/Cardioprotection/5_speclib_z6/report.tsv")

report_process(
  data="C:/Users/Laura/Downloads/Cardioprotection/5_speclib_z6/report.tsv",
  header.id = "Protein.Group",
  sample.id = "File.Name",
  quantity.id = "Precursor.Normalised",
  secondary.id = "Precursor.Id",
  id_to_add = c("Protein.Names", "First.Protein.Description", "Genes"),
  qv = 0.01,
  pg.qv = 0.01,
  p.qv = 1,
  gg.qv = 1,
  only_proteotypic = TRUE,
  get_pep = TRUE,
  only_pepall = TRUE,
  get_Top3 = TRUE,
  get_iBAQ = TRUE,
  fasta = c("C:/Users/Laura/Downloads/Cardioprotection/Rattus_norvegicus_ref-proteome_UP000002494_10116_22456prot_2024-07-30.fasta",
            "C:/Users/Laura/Downloads/Cardioprotection/0602_Universal Contaminants.fasta",
            "C:/Users/Laura/Downloads/Cardioprotection/Aug2022_Rat Tissue Contaminants.fasta"),
  species = NULL,
  peptide_length = c(7, 50),
  format = c("xlsx", "csv", "txt")
)

diann_Protein_Group <- read_excel("C:/Users/Laura/Downloads/Cardioprotection/5_speclib_z6/DIAgui/240803_2121_report_Protein.Group/diann_Protein.Group.xlsx")
diann_Protein_Group <- read_excel("C:/Users/Laura/Downloads/Cardioprotection/5_speclib_z6/DIAgui/240803_1327_ProteinGroup_dia.xlsx")

heatmapDIA(
  diann_Protein_Group,
  transformation = c(#"none",
                     #"log2",
                     #"z.score_fraction",
                     "z.score_protein"),
  maxna = 0,
  print_val = F,
  nm_id = "Protein.Group",
  data_type = c(#"Top3", 
                #"all",
                #"iBAQ",
                "intensity"),
  gradient_color = c("#09009D", "#ffffff", "#BE0010"),
  static = F
)
