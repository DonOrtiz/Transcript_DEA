#Transcriptomics Pipeline for Differential Expression Analysis and other transcriptomical visualizations
# Created by Emilio Fabian Ortiz

### Step 1: Load package libraries required for data import and analysis
library(DESeq2)
library(tximport)
library(ggplot2)
library(dplyr)

### Step 2: Import Salmon derived quantification files and sample information
#For this we first need to establish a directory that will contain the files
dir <- "D:\\Emili\\Desktop\\Tec BASE\\Transcriptomics Pipeline"
#you will need to set your working directory on the 'Session' tab
samples <-read.table(file.path(dir,"samples.txt"), header=TRUE) 
#NOTE: this samples.txt must contain the file names as refrence and must be the same as the files on the folder
#Gather the files given on you txt file
files <- file.path(dir, samples$Sample)
names(files) <- paste(samples$Sample)

#Step 3: create experimental design objects for DEA
#The lengthScaledTPM parameter will also perform a normalization of the transcripts given its lenght
txi <- tximport(files, type="salmon", TRUE, countsFromAbundance = "lengthScaledTPM",  dropInfReps= TRUE)
dds <- DESeqDataSetFromTximport(txi, colData= samples, design = ~Condition)

#Step 4: We run properly the analysis
ddsDE <- DESeq(dds) #Main differential express analysis
res <- results(ddsDE) #Extracts results table with log2 fold changes, p values and adjusted p values.
#Adjusting the results with an specific P value
res05 <- results(ddsDE, alpha=0.05)

#Step 5: Manage the GFF information
GFF_Treat <- read.delim("FSN2QEPIA.gff", header=F, comment.char="#") 
GFF_Cntrl <- read.delim("control_gff_export.gff3", header=F, comment.char="#") 
GFF_Cntrl$V9<-gsub("ID=","",as.character(GFF_Cntrl$V9))
GFF_Cntrl$V9<-gsub(";Description=RecName:","",as.character(GFF_Cntrl$V9))
GFF_Cntrl$V9<-gsub("---NA---","",as.character(GFF_Cntrl$V9))
GFF_Cntrl$V9<-gsub(";Description=","",as.character(GFF_Cntrl$V9))





