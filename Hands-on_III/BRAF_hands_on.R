#############################################################################
# BRAF practice                                                             #
# Bojos per la supercomputació - Life Sciences module                       #
# Barcelona Supercomputing Center                                           #
# Author : Victoria Ruiz                                                    #
# Contact: victoria.ruizserra@bsc.es                                        #
# Date : 05/10/2019                                                         #
#############################################################################

# Here you'll find all the necessary code to do perform a bioinformatics 
# analysis. More specifically, you'll adquire the knowledge to extract
# information from public data of real patients with cancer. 

# The final objective is to understand, process and interpret the given 
# data in such way of a normal routine of any researcher of our unit.

##########################
# Load necesary packages #
##########################
library(data.table)
library(ggplot2)
library(Biostrings)
library(geno2proteo)
library(stringr)
library(ggplot2)

##########################
# Set working directory  #
##########################
setwd("./")

###############################################################################
#                                                                             #
#                                                                             #
#                          PART 1- GENOMIC DATA                               #
#                                                                             #
#                                                                             #
###############################################################################

###############################################################################
# Vamos a utilizar el paquete de R "geno2proteo" para                         #
# localizar posiciones genómicas de mutaciones recurrentes                    #
# de pacientes reales en el gen BRAF                                          #
# Los siguientes pasos ya están precomputados (NO ES NECESARIO EJECUTARLO):   #
#  1) Descargar el genoma en formato .gtf                                     #
#     Para más información sobre este formato consultar el                    #
#     siguiente link                                                          #
#     (https://www.ensembl.org/info/website/upload/gff.html)                  #
#   ftp://ftp.ensembl.org/pub/grch37/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf
#
#  1.b) split the downloaded file into different chromosomes. In this way,the #
#     following code will run faster. To split it, type in terminal:          #
###############################################################################

#   awk -F ' ' 'NR==0 {h=$0; next} {f="./data/Homo_sapiens.GRChr38.90.Chr_"$1".gtf"}
#   !($1 in p) {p[$1]; print h > f} {print >>f ; close(f)}' ./data/Homo_sapiens.GRCh38.90.gtf

################################################################################
# 2) Descargar el genoma humano en formato .fasta (sólo el chr 7)              #
################################################################################

# download.file("ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz",
#               './scripts/data/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz')

################################################################################
# 3) Set the necessary parameters to execute the function                      #
# generatingCDSaaFile() from the geno2proteo package.                          #
#                                                                              #
# generatingCDSaaFile(): This function will find the DNA and protein           #
# sequences for each CDS region listed in the ENSEMBLgene annotation file      #
# (gtf file) provided, and store the CDS regions and the corresponding DNA     #
# andprotein sequences in an output data file. The output data file will be    #
# needed by some of the functionsin this package                               #
################################################################################

# dataFolder = system.file("extdata", package="geno2proteo")

# geneticCodeFile_line = file.path(dataFolder,"geneticCode_standardTable_lines.txt")

# gtfFile = "./data/Homo_sapiens.GRChr38.90.Chr_7.gtf"

# DNAfastaFile =  "./data/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz"

# outputFolder = "./output" # using the current folder as output folder

# calling the function.
# generatingCDSaaFile(geneticCodeFile_line=geneticCodeFile_line,
#                     gtfFile=gtfFile,
#                     DNAfastaFile=DNAfastaFile,
#                     outputFolder=outputFolder)

#filename00 = sub(".*/", "", gtfFile)

# get the output file's name.
# outputFile = paste(outputFolder, "/", filename00, "_AAseq.txt.gz", sep="")


###############################################################################
# CODE TO EXECUTE                                                             #
###############################################################################

# load PanCancer genomic data
#pancan_var_data = fread("./data/GDC-PANCAN.mutect2_snv.tsv")

# select data realted to the BRAF gene
#BRAF_pancan_var_data = subset(pancan_var_data,
#                              gene =="BRAF" &
#                              effect == "missense_variant",
#                              select = c("Sample_ID",
#                                         "gene", 
#                                         "chrom",
#                                         "start",
#                                         "end",
#                                         "ref",
#                                         "alt"))


# visualize the data 
BRAF_pancan_var_data = fread("./data/BRAF_pancan_var_data.txt")
BRAF_pancan_var_data

###############################################################################
# prepare data to execute genomicLocsToProteinSequence() function             #
#                                                                             #
# This function takes a list of genomic loci given in the input and tries to  #
# findthe protein sequences and DNA sequences of the coding                   #
# regions of genome which are within thosegenomic loci                        #
###############################################################################

# set the inputLoci data
inputLoci = BRAF_pancan_var_data[,c("chrom", "start", "end")]
inputLoci$strand = "-"

# set the CDSaaFile
CDSaaFile= "./output/Homo_sapiens.GRChr38.90.Chr_7.gtf_AAseq.txt.gz"

# get the codon
proteinSeq = genomicLocsToProteinSequence(inputLoci=inputLoci, 
                                          CDSaaFile=CDSaaFile)

# get the codons 
codons = paste(proteinSeq$dnaBefore, proteinSeq$dnaSeq, proteinSeq$dnaAfter, sep ="")

# alternative codons
alt_codons = paste(proteinSeq$dnaBefore, BRAF_pancan_var_data$ref, proteinSeq$dnaAfter, sep ="")

# translate codons
aa     = as.data.frame(translate(DNAStringSet(codons)))$x
alt_aa = as.data.frame(translate(DNAStringSet(alt_codons)))$x

# update data frame
BRAF_pancan_var_data$codons     = codons
BRAF_pancan_var_data$alt_codons = alt_codons
BRAF_pancan_var_data$aa         = aa
BRAF_pancan_var_data$alt_aa     = alt_aa

###############################################################################
#                                                                             #
#                                                                             #
#                          PART 2 - PROTEIN DATA                              #
#                                                                             #
#                                                                             #
###############################################################################

# Load TCGA data
pancan_data = fread("./data/Curated_PanCanAtlas_clinical_data.txt")

# visualize the nature of the loaded data
pancan_data

# put the patients ids into the rigth format
BRAF_pancan_var_data$Sample_ID = strtrim(BRAF_pancan_var_data$Sample_ID, 12)

# add data regarding cancer type to the BRAF mutations data frame 
pancan_type = pancan_data[, c("bcr_patient_barcode", "type")] #subset data
names(pancan_type)[1] = "Sample_ID"  #change colID

# merge both datasets taking as reference the patient ID
BRAF_pancan_var_data = merge(pancan_type,
                             BRAF_pancan_var_data,
                             by ="Sample_ID") 

# add the data of mutations
BRAF_pancan_var_data$mut = paste(BRAF_pancan_var_data$aa
                                 , BRAF_pancan_var_data$alt_aa, sep = "/")

BRAF_pancan_var_data # take a look of the data now

################################################################################
# How many unique mutations do we find of this gene among the patients?        #
################################################################################

# extract from the column "mut" the uniqued mutations
BRAF_unique_mutations = unique(BRAF_pancan_var_data$mut)

# take a look of the most common amino acids mutations
BRAF_unique_mutations 

# compute the total number of diferent mutations
length(BRAF_unique_mutations) 

# Answer: ---- mutations in total

###############################################################################
# among all the TCGA patients, how many of them do have mutations on BRAF     #
###############################################################################
# compute the total number of patients with patients
length(unique(BRAF_pancan_var_data$Sample_ID)) 
# Answer: --- patients do have mutations 

# The data frame 'BRAF_pancan_var_dara' stores more rows than the 
# number of patients, how is that possible?
# Answer : 

###############################################################################
# What is the distribution of this mutations?                                 #
###############################################################################
ggplot(BRAF_pancan_var_data, aes(x =mut))+
  geom_bar(fill="white", position="dodge", stat = "count", color ="orange")+
  ggtitle("Most common amino acid change")

# the most recurrent mutation is ------

###############################################################################
# For the most common mutation, what is the most common cancer type           #
###############################################################################
# Please, below, the most common amino acid change between quotation marks
most_common_mutation = "" #write the most common amino acid change, e.g.: "L/A"
patients_most_common_mutation = subset(BRAF_pancan_var_data,
                                mut == most_common_mutation)

ggplot(patients_most_common_mutation, aes(x =type, fill = type))+
  geom_bar()+
  ggtitle("Distribution of cancer types of most common amino acid change")

# what is/are the most common cancer types containing the 
# mutation?
# Answer: -------------
# Note: if you don't know what are the meaning of the acronysms, look for
# them on google!


###############################################################################
#                                                                             #
#                                                                             #
#                          PART 3 - MACHINE LEARNING                          #
#                                                                             #
#                                                                             #
###############################################################################
# The code below is the actual code used to generate our machine 
# learning example. Due to lack of time, we will not comment the method applied
# here. However, if you are courious, you can inspect the code and test it 
# yourself! :)

# Source of data "https://clincancerres.aacrjournals.org/content/25/11/3239"

##########################
# Load necesary packages #
##########################
library(caret)
library(foreach)
library(doParallel)

# Load data
braf_data = read.csv("./data/198021_2_supp_5455023_ppszzf.csv",
                     sep = ",")
# variables to train the different models
opts = c("BRAF.V600.Mut","Rx", "BMI", "Age", "Sex", "Stage", "LDH")

# generate a list of all posible combinations of all the variables 
# to generate the diferent models through iterations. 127 in total. 
l = list()
for(i in 1:7){
  l[[i]] =combn(opts, i) 
}

# create an input data frame to feed the machine learning alforithm. 
input_ML = braf_data[, c("Patient.ID","BORR", opts)]
# change column ID
names(input_ML)[2]= "phenotype"
input_ML = data.table(input_ML, stringsAsFactors = T)

# change working directory
setwd("./models")

# iterate over all the possible combinations of all the variables
for ( i in 1:length(l) ){
  for (j in 1:ncol(l[[i]])){
    
    set.seed(0) # for the sake of reproducibility
    
    # train and test data
    input_ML_subset = subset(input_ML, select =  c("Patient.ID","phenotype", l[[i]][,j]))
    # create 80% train data 20% test data
    partition = createDataPartition(input_ML_subset$phenotype, p=0.8, list=FALSE)
    # remove patient id column from input data
    input_ML_subset = input_ML_subset[,-"Patient.ID"]
    
    # create test and train datasets and remove NA values
    # since it will give errors messages 
    # !!!!Warning: in real life, the removal of NA values 
    # is not necessarily a good practice. Since this is 
    # an example, we do not really care. However, when we decide
    # to remove data from the original source we have to justify it
    # and do it carefully.  
    test_data = input_ML_subset[-partition,]
    train_data = input_ML_subset[partition,]
    train_data = na.omit(train_data)
    

    # select parameters to train the model
    # 10-fold crossvalidation is a recommended practice. 
    cctrl1 <- trainControl(method = "cv",
                           number = 10,
                           allowParallel = TRUE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE
    )
    
    # Option to run in parallel: 
    #cl <- makePSOCKcluster(4)
    #registerDoParallel(cl)
    
    # build and train the model
    # we have selected Random Forest model because is 
    # the "easiest" to understand and to compute!. 
    my_model <- train(phenotype ~ .,
                      data = train_data, 
                      method = 'parRF',
                      trControl = cctrl1,
                      na.action = na.pass)
    #stopCluster(cl)
    
    # save the model to use it in the shiny app
    saveRDS(my_model , file = paste('./models/model_predictors_', paste(l[[i]][,j],collapse = "_"),'.rds', sep = ""))
    
    # save the confusion matrix to use it in the shiny app
    cm <- confusionMatrix(predict(my_model, test_data), test_data[[1]]) 
    saveRDS(cm , file = paste('./models/confusionMatrix_predictors_', paste(l[[i]][,j],collapse = "_"),'.rds', sep = ""))
  }
}
