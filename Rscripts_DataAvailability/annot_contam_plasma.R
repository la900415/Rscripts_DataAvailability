library(tidyverse)
library(dplyr)
library(readr)
library(stringr)

# Define contaminant lists with modifications according to Gao et al. Systematic evaluation of blood contamination in nanoparticle-based plasma proteomics. bioRxiv 2025, 2025.2004.2026.650757. DOI: 10.1101/2025.04.26.650757.
platelet_markers <- c("PF4", "PF4V1", "PPBP", "SELP", "THBS1", "ITGA2B", "ITGB3", "GP1BA", "GP9", "TBXAS1")
coagulation_markers <- c("FGA", "FGB", "FGG", "F2", "F10", "F11", "F12", "F13A1", "F13B", "VWF", "SERPINC1","SERPINF2", "VTN", "HRG")
rbc_markers <- c("HBA1","HBA2","HBB","HBD", "HBG1","HBG2","HBZ", "ALAS2", "SLC4A1", "CA1", "ANK1", "EPB41", "EPB42", "GYPA", "GYPB","GYPC", "SPTA1", "SPTB", "ALDOA", "ENO1", "TPI1", "PGK1", "PKM", "LDHA")

# Function to annotate contaminants
annotate_contaminants <- function(df, gene_col, descript_col, id_col) {
  gene_vec <- pull(df, gene_col)
  descript_vec <- pull(df, descript_col)
  id_vec <- pull(df, id_col)
  
  df %>%
    mutate(
      Potential.contaminant = case_when(
        gene_vec %in% platelet_markers ~ "+",
        gene_vec %in% coagulation_markers ~ "+",
        gene_vec %in% rbc_markers ~ "+",
        str_detect(descript_vec, regex("Keratin, type", ignore_case = TRUE)) ~ "+",
        str_detect(descript_vec, regex("variable", ignore_case = TRUE)) &
          str_detect(descript_vec, regex("immunoglobulin", ignore_case = TRUE)) ~ "+",
        str_detect(id_vec, regex("Cont_", ignore_case = TRUE)) ~ "+",
        TRUE ~ ""
      ),
      Source.contaminant = case_when(
        gene_vec %in% platelet_markers ~ "platelet",
        gene_vec %in% coagulation_markers ~ "Coagulation",
        gene_vec %in% rbc_markers ~ "Rbd",
        str_detect(descript_vec, regex("Keratin, type", ignore_case = TRUE)) ~ "Keratin",
        str_detect(descript_vec, regex("variable", ignore_case = TRUE)) &
          str_detect(descript_vec, regex("immunoglobulin", ignore_case = TRUE)) ~ "Immunoglobulin",
        str_detect(id_vec, regex("Cont_", ignore_case = TRUE)) ~ "common.cont",
        TRUE ~ ""
      )
    )
}


# 🧪 Cargar archivo FragPipe manualmente
fragpipe_path <- file.choose()
fragpipe <- read_tsv(fragpipe_path)
fragpipe_annotated <- annotate_contaminants(fragpipe, gene_col = "Gene", descript_col = "Description", id_col = "Protein ID")
# 🧪 Cargar archivo AMICA-FragPipe manualmente
fragpipe_annotated <- annotate_contaminants(fragpipe, gene_col = "Gene.names", descript_col = "Description", id_col = "Protein.ID")
fragpipe_annotated <- fragpipe_annotated %>%
  mutate(quantified = if_else(Potential.contaminant == "+" & quantified == "+", NA_character_, quantified))

# Guardar archivo FragPipe anotado
fragpipe_save_path <- file.choose(new = TRUE)
write_tsv(fragpipe_annotated, fragpipe_save_path)
write_tsv(fragpipe_annotated, "combined_protein.annotated.tsv")
write_tsv(fragpipe_annotated, "2_amica_protein_groups_min1pept_imputed.annotated.tsv")

# 🧬 Cargar archivo DIA-NN manualmente
diann_path <- file.choose()
diann <- read_tsv(diann_path)
diann_annotated <- annotate_contaminants(diann, gene_col = "Genes", descript_col = "First.Protein.Description", id_col = "Protein.Group")

# Guardar archivo DIA-NN anotado
diann_save_path <- file.choose(new = TRUE)
write_tsv(diann_annotated, diann_save_path)

