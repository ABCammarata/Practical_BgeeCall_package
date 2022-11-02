setwd("./")

#load BgeeCall package
library(BgeeCall)
library(R.utils)
library(tidyverse)

#Download annotation and transcriptome files for Drosophila melanogaster
annotation_file_path =  "./Data/Drosophila_melanogaster_Annotation.chr.gtf.gz"
download.file("http://ftp.ensemblgenomes.org/pub/release-51/metazoa/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.51.chr.gtf.gz", annotation_file_path)
gunzip(annotation_file_path)


transcriptome_file_path = ".//Data/Drosophila_melanogaster_Transcriptome.cdna.all.fa.gz"
download.file("http://ftp.ensemblgenomes.org/pub/release-51/metazoa/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz", transcriptome_file_path)
gunzip(transcriptome_file_path)


###Retrieve intergenic information
##1) List all intergenic releases available in BgeeCall. How many exist?
list_intergenic_release()
#There are 5 intergenic releases that are available in BgeeCall, 3 provided by Bgee, one containing the ones done by the community and the possibility to create a custom one.

##2)Verify which species are available for the current Bgee intergenic release. How many exist?
bgee <- new("BgeeMetadata", intergenic_release="1.0")
species_list <- list_bgee_ref_intergenic_species(myBgeeMetadata = bgee)
#List of the available species
species_available <- species_list$speciesName
#There are 52 available species in the current Bgee intergenic release.

##3) Verify which species belong to the community. How many exist?
list_community_ref_intergenic_species()
#There are two species that belong to the community. We can see that they are Mesocricetus auratus
#and Solenopsis invicta (using the speciesId and NCBI taxonomy Browser).

##Use BgeeCall to download the pseudo-alignment software
#1)Create an object of the KallistoMetadata class.
#creation of the kallisto object from the KallistoMetadata class
kallisto <- new("KallistoMetadata", download_kallisto = TRUE)

###Run analysis: Drosophila melanogaster 1 sample
##1)Create a userMetadata object (note that you have to specify in the argument species_id the Taxonomy ID, you can verify that in https://bgee.org/ in the See species information). 
#Check for Drosophila melanogaster availability and speciesId
searched_speciesId <- species_list$speciesId[species_list$speciesName == "Drosophila melanogaster"]
#7227 %in% speciesID_available
user_BgeeCall <- new("UserMetadata", species_id = "7227", reads_size = 100)

##2)What happens if the argument reads_size is not specified by you when you create the new userMetadata object? What can be the impact in general? Reads size of RNA-Seq libraries can be found in SRA (e.g https://www.ncbi.nlm.nih.gov/sra/?term=SRX109278)
#The read_size will be set to 51 by default. Specifying the wrong read_size may interfer withe the alignment process.

##3)Specify by using the following functions setRNASeqLibPath(), setTranscriptomeFromFile(), setAnnotationFromFile(), setOutputDir() and setWorkingPath() the path to your library SRX109278, transcriptome file, annotation file as well as the output and working directory. 
#unzip RNA_seq libraries

user_BgeeCall <- setRNASeqLibPath(user_BgeeCall, "./RNA_seq/SRX109278/")
user_BgeeCall <-setTranscriptomeFromFile(user_BgeeCall, "./Data/Drosophila_melanogaster_Transcriptome.cdna.all.fa")
user_BgeeCall <- setAnnotationFromFile(user_BgeeCall, "./Data/Drosophila_melanogaster_Annotation.chr.gtf")
user_BgeeCall <- setWorkingPath(user_BgeeCall, "./Results/")
user_BgeeCall <- setOutputDir(user_BgeeCall, "./Output/SRX109278/cutoff_0.05")

##4)Generate the present and absent calls for the library SRX109278 by using generate_calls_workflow(). Which type of information is provided in the output files?

calls_output <- generate_calls_workflow(abundanceMetadata = kallisto, userMetadata = user_BgeeCall)

head(read.table(calls_output$calls_tsv_path, header = TRUE), n = 5)
#In this file, we can see a multitude of information such as the presence/absence call, the counts number, the pvalues, the type of the region (genic, intergenic).

read.table(calls_output$cutoff_info_file_path)
#Here, we can read information summarizing the full data, such as the proportion of different regions and the cutoffTPM.

head(read.table(calls_output$abundance_tsv, header = TRUE), n = 5)
#In this file, we find the output of kallisto's execution.

read.table(calls_output$S4_slots_summary, header = TRUE, sep = "\t")
#We can see the values used for our call workflow. 

#We can also see a pdf showing the density distributions of the TPM values.

##5)Plot the frequency of p-values for the correspondent library.
call_data <- read.table(calls_output$calls_tsv_path, header = TRUE)

jpeg("./Plots/p-value_frequencies.jpg", width = 700, height = 700)
par(mfrow=c(1,2))
hist(call_data$pValue,
     col = "lightblue",
     breaks= 20,
     main = "Frequency of p-values for SRX109278 library",
     xlab = "p-Value",
     ylab = "Frequency")

hist(call_data$pValue,
     col = "lightblue",
     xlim=c(0,0.05),
     breaks= 200,
     main = "Frequency of p-values for SRX109278 library",
     xlab = "p-Value",
     ylab = "Frequency")
dev.off()

###Run analysis: multiple Drosophila melanogaster samples
##1)Create a user input file describing all RNA-Seq libraries previously downloaded, see https://github.com/BgeeDB/BgeeCall/blob/develop/inst/userMetadataTemplate.tsv and the vignette of the package for more information
#File Path: ./UserMetadata/Default_values.tsv

##2)Run the generation of present and absent calls from the user file with default values for all .
user_Metadata_file = read.table("./UserMetadata/Default_values.tsv", header = TRUE, sep = "\t")
calls_output_multiple_samples <- generate_calls_workflow(abundanceMetadata = kallisto, userFile = "./UserMetadata/Default_values.tsv")

#We can see the first few lines of the results for one of the RNA_seq_libraries
head(read.table(calls_output_multiple_samples[[1]][1]$calls_tsv_path, header = TRUE), n = 5)


##3)Combine multiple libraries per species using the merging_libraries() function. What is the proportion of genes present?
mergedLibraries <- merging_libraries(userFile = "./UserMetadata/Default_values.tsv", approach = "BH", condition = "species_id", cutoff = 0.05, outDir = "./Output/CallbySpecies/cutoff_0.05")

mergedlib <- read.table("./Output/CallbySpecies/cutoff_0.05/Calls_merging_BH_cutoff=0.05_species_id=7227.tsv", header = TRUE)
proportion <- count(mergedlib$call == "present")/(count(mergedlib$call == 'absent') + count(mergedlib$call == 'present'))
print(proportion)
#We can see that 80.5 percent of the genes are present.

##4)Modify the input file to combine libraries per species (species_id) and developmental stage (devStage)
#Modified input file path: ./UserMetadata/BySpecies+devStage.tsv

##5)Generate the present and absent calls with a more restrictive p-value = 0.01
user_Metadata_merged_file = read.table("./UserMetadata/BySpecies+devStage.tsv", header = TRUE, sep = "\t")

#Generate calls for all libraries with the more restrictive p-value (0.01)
kallisto <- new("KallistoMetadata", download_kallisto = TRUE, cutoff = 0.01)
calls_output_multiple_samples_0.01 <- generate_calls_workflow(abundanceMetadata = kallisto, userFile = "./UserMetadata/More_Restrictive_Default_0.01.tsv")

#Get minimal p-value by grouping by species and developmental Stage
mergedSpeciesDevStage <- merging_libraries(userFile = "./UserMetadata/BySpecies+devStage.tsv", approach = "BH", condition = c("species_id", "devStage"), cutoff = 0.01, outDir = "./Output/CallbySpecies+devStage/cutoff_0.01")

#Open the minimal p-values grouped in tables
Uberon <- read.table("./Output/CallbySpecies+devStage/cutoff_0.01/Calls_merging_BH_cutoff=0.01_species_id=7227_devStage=UBERON:0000066.tsv", header = TRUE)
FBdv <- read.table("./Output/CallbySpecies+devStage/cutoff_0.01/Calls_merging_BH_cutoff=0.01_species_id=7227_devStage=FBdv:00007079.tsv", header = TRUE)

##6)Get summary stats of all libraries by using get_summary_stats() function.

summaryStats <- get_summary_stats(userFile = "./UserMetadata/Default_all_cutoffs.tsv", outDir = "./Output/Summary")


all_libraries_summary <- read.table("./Output/Summary/summary_Stats_All_Libraries.tsv", header = TRUE)

##7)Plot the proportion of protein coding genes of all libraries for each p-value cut-off.
jpeg("./Plots/protein_coding_genes_proportion1.jpg", width = 700, height = 700)
ggplot(all_libraries_summary, aes(libraryId, proportionCodingPresent, fill=libraryId))+
  geom_bar(stat='identity')+
  facet_wrap(~cutoff)+
  ylab("Proportion of protein coding genes present")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()  

jpeg("./Plots/protein_coding_genes_proportion2.jpg", width = 700, height = 700)
ggplot(all_libraries_summary, aes(cutoff, proportionCodingPresent, fill=as.factor(cutoff)))+
  geom_bar(stat='identity')+
  facet_wrap(~libraryId)+
  ylab("Proportion of protein coding genes present")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


###Downstream analysis
##1)Perform a differential expression analysis between different developmental stage conditions.
#Redoing calculation with same number of samples from the two developmental stages
user_Metadata_balanced_file = read.table("./UserMetadata/Balanced_devStage.0.05.tsv", header = TRUE, sep = "\t")
mergedDevStage_balanced <- merging_libraries(userFile = "./UserMetadata/Balanced_devStage.0.05.tsv", approach = "BH", condition = c("species_id", "devStage"), cutoff = 0.05, outDir = "./Output/CallbySpecies+devStage/Balanced/cutoff_0.05")

Uberon <- read.table("./Output/CallbySpecies+devStage/Balanced/cutoff_0.05/Calls_merging_BH_cutoff=0.05_species_id=7227_devStage=UBERON:0000066.tsv", header = TRUE)
FBdv <- read.table("./Output/CallbySpecies+devStage/Balanced/cutoff_0.05/Calls_merging_BH_cutoff=0.05_species_id=7227_devStage=FBdv:00007079.tsv", header = TRUE)

read_expressed_in_uberon <- Uberon$id[Uberon$call == "present"]
read_expressed_in_FBdv <- FBdv$id[FBdv$call == "present"]

difference_uberon <- setdiff(read_expressed_in_uberon, read_expressed_in_FBdv)
difference_FBdv <- setdiff( read_expressed_in_FBdv, read_expressed_in_uberon)
#We can see that there are a variety of different genes that are expressed in one stage but not the other.

##2)Filter the results by providing just genes with FDR < 0.01. Provide a visualization graphic as MA plot. 
#Make presence absence calls using the qvalue approach with a cutoff of 0.01
kallisto_qvalue <- new("KallistoMetadata", download_kallisto = TRUE, cutoff = 0.01, cutoff_type = "qValue")
calls_output_multiple_samples_qvalue0.01 <- generate_calls_workflow(abundanceMetadata = kallisto_qvalue, userFile = "./UserMetadata/Q_value_approach_0.01.tsv")

#Merging based on developmental stage
user_Metadata_2_2_qvalue_file = read.table("./UserMetadata/Balanced_devStage_qvalue_0.01.tsv", header = TRUE, sep = "\t")
mergedDevStage_2_2_qvalue <- merging_libraries(userFile = "./UserMetadata/Balanced_devStage_qvalue_0.01.tsv", approach = "fdr_inverse", condition = c("species_id", "devStage"), cutoff = 0.01, outDir = "./Output/CallbySpecies+devStage/Balanced/qvalue/cutoff_0.01")

Uberon_q <- read.table("./Output/CallbySpecies+devStage/Balanced/qvalue/cutoff_0.01/Calls_merging_fdr_inverse_cutoff=0.01_species_id=7227_devStage=UBERON:0000066.tsv", header = TRUE)
FBdv_q <- read.table("./Output/CallbySpecies+devStage/Balanced/qvalue/cutoff_0.01/Calls_merging_fdr_inverse_cutoff=0.01_species_id=7227_devStage=FBdv:00007079.tsv", header = TRUE)

genes_uberon_q <- Uberon_q$id[Uberon_q$call == "present"]
genes_FBdv_q <- FBdv_q$id[FBdv_q$call == "present"]

difference_uberon_q <- setdiff(genes_uberon_q, genes_FBdv_q)
difference_FBdv_q <- setdiff( genes_FBdv_q, genes_uberon_q)

#We can make the MA plot
#We start by reading the tables containing the reads counts.
abundance_FBdv1 <- read.table(calls_output_multiple_samples_qvalue0.01[[1]][1]$calls_tsv_path, header = TRUE)
abundance_FBdv2 <- read.table(calls_output_multiple_samples_qvalue0.01[[2]][1]$calls_tsv_path, header = TRUE)
abundance_uberon1 <- read.table(calls_output_multiple_samples_qvalue0.01[[5]][1]$calls_tsv_path, header = TRUE)
abundance_uberon2 <- read.table(calls_output_multiple_samples_qvalue0.01[[6]][1]$calls_tsv_path, header = TRUE)

#We calculate the mean for the RNA libraries involved with both developmental stages.
abundance_FBdv1$counts <- log(abundance_FBdv1$counts)
abundance_FBdv1$counts2 <- log(abundance_FBdv2$counts)
abundance_FBdv1$mean <- rowMeans(abundance_FBdv1[, c('counts', 'counts2')], na.rm=TRUE)

abundance_uberon1$counts <- log(abundance_uberon1$counts)
abundance_uberon1$counts2 <- log(abundance_uberon2$counts)
abundance_uberon1$mean <- rowMeans(abundance_uberon1[, c('counts', 'counts2')], na.rm=TRUE)

#Create MA dataframe
#Calculate the average expression and the difference of expression between the develpmental stages.
A = (abundance_FBdv1$mean + abundance_uberon1$mean)/2
M = abundance_FBdv1$mean - abundance_uberon1$mean
MA <- data.frame(M, A)

#Remove -infinity and NA values
MA$M[!is.finite(MA$M)] <- NA
MA$A[!is.finite(MA$A)] <- NA
MA <- na.omit(MA)


jpeg("./Plots/MA_plot.jpg", width = 700, height = 700)
ggplot(MA, aes(A, M))+
  geom_point(colour = "lightblue")+
  xlab("log2(Average count number)")+
  ylab("log2(count number difference)")
dev.off()

jpeg("./Plots/MA_plot_zoomed.jpg", width = 700, height = 700)
ggplot(MA, aes(A, M))+
  geom_point(colour = "lightblue")+
  xlim(0,12)+
  ylim(-10, 10)+
  xlab("log2(Average count number)")+
  ylab("log2(count number difference)")
dev.off()

#Average count lower than 0
MA_neg <- MA[MA$A < 0, ]

##3)Make a GO analysis. Which type of information do you retrieve? 
#We can look back at the genes that were expressed differently in the two developmental stages (difference_uberon_q vs difference_FBdv_q and difference_Uberon vs difference_FBdv)
#We will check if the gene expressed differently in the 2 developmental stages have particular function. This will be done using the gene ontology website (http://geneontology.org/)

write.table(difference_uberon_q, "./Data/GO/differencial_gene_expression/Uberon_genes_qvalue.txt", sep = "\n",row.names = FALSE, col.names = FALSE, quote = FALSE)
#We see that that genes affecting the catalytic activity on nucleic acids are underrepresented compared to the FBdv stage

write.table(difference_FBdv_q, "./Data/GO/differencial_gene_expression/FBdv_genes_qvalue.txt", sep = "\n",row.names = FALSE, col.names = FALSE, quote = FALSE)
#We can see that there is a multitude of function that are over or under-represented in the gene_set, such as mating function that  are over-represented and developmental function that are under-represented.

write.table(difference_uberon, "./Data/GO/differencial_gene_expression/Uberon_genes.txt", sep = "\n",row.names = FALSE, col.names = FALSE, quote = FALSE)
#We see that only genes associated with cellular anatomical entitiy seem to be under-represented.

write.table(difference_FBdv, "./Data/GO/differencial_gene_expression/FBdv_genes.txt", sep = "\n",row.names = FALSE, col.names = FALSE, quote = FALSE)
#Again we see an over-representation of genes related to mating for females.
