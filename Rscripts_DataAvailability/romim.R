library(romim)

################### wrapper functions to query the OMIM API ##########################
#install if necessary
install.packages("XML")

#load library
library(XML)

#function to set API key
#the my_key object becomes a global object
set_key <- function(key){
  my_key <- paste('apiKey=', key, sep='')
}

#function to build the URL
#and to perform the query
get_omim <- function(omim_id = 100100, #default OMIM ID
                     text = FALSE, #Includes the text field sections with the entry.
                     existflags = FALSE, #Include the 'exists' flags with the entry (clinical synopsis, allelic variant, gene map & phenotype map).
                     allelicVariantList = FALSE, #Includes the allelic variant list with the entry.
                     clinicalSynopsis = FALSE, #Include the clinical synopsis with the entry.
                     seeAlso = FALSE, #Includes the 'see also' field with the entry.
                     referenceList = FALSE, #Include the reference list with the entry.
                     geneMap = FALSE, #Include the gene map/phenotype map data with the entry.
                     externalLinks = FALSE, #Include the external links with the entry.
                     contributors = FALSE, #Includes the 'contributors' field with the entry.
                     creationDate = FALSE, #Includes the 'creation date' field with the entry.
                     editHistory = FALSE, #Includes the 'edit history' field with the entry.
                     dates = FALSE, #Include the dates data with the entry.
                     all = FALSE #Include the above data with the entry.
){
  #get all the arguments of the function call
  a <- as.list(match.call())
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, my_key, sep = "&")
  #loop through all the arguments
  for (i in names(a)){
    #skip the omid_id and blank argument
    if(!i %in% '' && !i %in% 'omim_id'){
      my_include <- paste('&', 'include=', i, sep='')
      my_query <- paste(my_query, my_include, sep='')
    }
  }
  xmlParse(my_query)
}

#function to parse the XML
#for the title of an OMIM ID
get_title <- function(xml){
  xml_list <- xmlToList(xml)
  return(xml_list$entryList$entry$titles$preferredTitle)
}

###################  Here's a demonstration of the 3 functions ########################
set_key('NOT_A_REAL_KEY')

#I made the default OMIM ID 100100
#for the get_omim() function
get_omim()
<?xml version="1.0" encoding="UTF-8"?>
  <omim version="1.0">
  <entryList>
  <entry>
  <prefix>#</prefix>
  <mimNumber>100100</mimNumber>
  <status>live</status>
  <titles>
  <preferredTitle>ABDOMINAL MUSCLES, ABSENCE OF, WITH URINARY TRACT ABNORMALITY AND CRYPTORCHIDISM</preferredTitle>
  <alternativeTitles>PRUNE BELLY SYNDROME;;
EAGLE-BARRETT SYNDROME; EGBRS</alternativeTitles>
  </titles>
  </entry>
  </entryList>
  </omim>
  
  get_omim(geneMap=TRUE)
<?xml version="1.0" encoding="UTF-8"?>
  <omim version="1.0">
  <entryList>
  <entry>
  <prefix>#</prefix>
  <mimNumber>100100</mimNumber>
  <status>live</status>
  <titles>
  <preferredTitle>ABDOMINAL MUSCLES, ABSENCE OF, WITH URINARY TRACT ABNORMALITY AND CRYPTORCHIDISM</preferredTitle>
  <alternativeTitles>PRUNE BELLY SYNDROME;;
EAGLE-BARRETT SYNDROME; EGBRS</alternativeTitles>
  </titles>
  <phenotypeMapList>
  <phenotypeMap>
  <mimNumber>118494</mimNumber>
  <phenotype>?Eagle-Barrett syndrome</phenotype>
  <phenotypeMimNumber>100100</phenotypeMimNumber>
  <phenotypeMappingKey>3</phenotypeMappingKey>
  <phenotypeInheritance>Autosomal recessive</phenotypeInheritance>
  <sequenceID>1484</sequenceID>
  <chromosome>1</chromosome>
  <chromosomeSymbol>1</chromosomeSymbol>
  <chromosomeSort>1484</chromosomeSort>
  <chromosomeLocationStart>239549875</chromosomeLocationStart>
  <chromosomeLocationEnd>240078749</chromosomeLocationEnd>
  <transcript>uc001hyp.3</transcript>
  <cytoLocation>1q41-q44</cytoLocation>
  <computedCytoLocation>1q43</computedCytoLocation>
  <geneSymbols>CHRM3, EGBRS</geneSymbols>
  <geneInheritance/>
  </phenotypeMap>
  </phenotypeMapList>
  </entry>
  </entryList>
  </omim>
  
  get_omim(100200)
<?xml version="1.0" encoding="UTF-8"?>
  <omim version="1.0">
  <entryList>
  <entry>
  <mimNumber>100200</mimNumber>
  <status>live</status>
  <titles>
  <preferredTitle>ABDUCENS PALSY</preferredTitle>
  </titles>
  </entry>
  </entryList>
  </omim>
  
  get_title(get_omim(100200))
[1] "ABDUCENS PALSY"

################### Here's how to look up a list of 10 OMIM IDs ########################
my_list <- list(100070,100100,100300,100600,100800,100820,101000,101200,101400,101600)
my_list_title <- lapply(lapply(my_list, get_omim), get_title)

df <- data.frame(id=unlist(my_list),
                 title=unlist(my_list_title))
df
1  100070                                     AORTIC ANEURYSM, FAMILIAL ABDOMINAL, 1; AAA1
2  100100 ABDOMINAL MUSCLES, ABSENCE OF, WITH URINARY TRACT ABNORMALITY AND CRYPTORCHIDISM
3  100300                                                    ADAMS-OLIVER SYNDROME 1; AOS1
4  100600                                                             ACANTHOSIS NIGRICANS
5  100800                                                              ACHONDROPLASIA; ACH
6  100820                                                                   ACHOO SYNDROME
7  101000                                                  NEUROFIBROMATOSIS, TYPE II; NF2
8  101200                                                                   APERT SYNDROME
9  101400                                                    SAETHRE-CHOTZEN SYNDROME; SCS
10 101600                                                                PFEIFFER SYNDROME


################### Example of limits ###########################################
my_table <- read.table("mim2gene.txt",header=F,sep="\t")
dim(my_table)
head(my_table)

system.time(my_big_list_title <- lapply(lapply(my_table$V1, get_omim), get_title))
#almost 3 hours

df <- data.frame(id=my_table$V1,
                 title=unlist(my_big_list_title))

#name the columns to prepare for merging
colnames(my_table) <- c('id',
                        'type',
                        'entrez_gene_id',
                        'hgnc_symbol',
                        'ensembl_gene_id')

#merge the original table with the OMIM titles
new_table <- merge(x=my_table, y=df, by='id')
dim(new_table)
head(new_table)
tail(new_table)
write.table(new_table, file='mim2gene_title.txt', quote=F, sep="\t", row.names=F)

################### Using the romim package #####################################
library(romim)
setwd("C:/Users/Laura/Downloads/RENOMEII/OMIM_HPO_database")

my_table <- read.table("mim2gene.txt",header=F,sep="\t")
head(my_table)
