CX.Dawn_Enriquecimiento <- read.delim2("C:/Users/Luis Ariel/Downloads/CX-Dawn_Enriquecimiento.txt")
CX.UP_Enriquecimiento <- read.delim2("C:/Users/Luis Ariel/Downloads/CX-UPEnriquecimiento.txt")

writexl::write_xlsx(CX.Dawn_Enriquecimiento, "C:/Users/Luis Ariel/Downloads/CX-Dawn_Enriquecimiento.xlsx")
writexl::write_xlsx(CX.UP_Enriquecimiento, "C:/Users/Luis Ariel/Downloads/CX-UP_Enriquecimiento.xlsx")
