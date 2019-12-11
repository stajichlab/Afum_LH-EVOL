#This script looks at specific mutations across the Pbs2 gene across the whole A.fum tree (V8 containing 356 isolates)

###first in bash, here's what I did:
#location of PBS2 / Afu1g15950
#I put this into a bed file pbs2.bed
#(with tabs as the seperator):

#Chr1_A_fumigatus_Af293	4327904	4330995   

#the version you want is V.8.
#A_fumigiatus_Af293.Popgen8.filtered.SNP.vcf


#module load bcftools

#search in the form 
#bcftools view -R bedfile.bed vcffile (remember to use the zipped version).

#bcftools view -R pbs2.bed A_fumigiatus_Af293.Popgen8.SNP.vcf.gz > SNPS_in_PBS2.vcf

#remake with filtered 
#bcftools view -R pbs2.bed A_fumigiatus_Af293.Popgen8.filtered.SNP.vcf.gz > SNPS_in_PBS2.vcf

#print to a usable format
#bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%AF[\t%TGT]\n' SNPS_in_PBS2.vcf> SNPS_in_PBS2_b.tab

#NOTE - you can grep to PASS here - but lets do it in R to see how many failed.


#load packages
require(data.table)
require(splitstackshape)
library(plyr)
library(dplyr)
library(tidyverse)

#set dir
setwd("~/Desktop/Afum_AF100_Clades/data/")
options(stringsAsFactors = FALSE)

SNPs_in_PBS2<-read.delim("SNPS_in_PBS2_filter.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
#subset to only PASS
SNPs_in_PBS2_Pass<- SNPs_in_PBS2[SNPs_in_PBS2$X.6.FILTER == "PASS", ]

###Return only collums that do not match the reference

#subset 
samples<- SNPs_in_PBS2_Pass[,c(3,7:ncol(SNPs_in_PBS2_Pass))]

#convert to T/F dataframe
are_any_true<- as.data.frame(lapply(samples[], `%in%`, samples$X.3.REF))
are_any_true

#Get only collums where some are true 
only_valid_varriants<-SNPs_in_PBS2_Pass[,apply(are_any_true,2,function(x) !all(x == TRUE))]

dim(SNPs_in_PBS2_Pass)
dim(only_valid_varriants)
#324 out of 362 strains have mutations that differ from the reference SOMEWHERE in pbs2. 

#look at the specific mutations that Brandon is interested in for PBS2
#positions 4329162 (Clade 1) and 4329159 (Clade 2)
pos4329162_C1<- SNPs_in_PBS2[SNPs_in_PBS2$X.2.POS == 4329162,]
pos4329159_C2<- SNPs_in_PBS2[SNPs_in_PBS2$X.2.POS == 4329159,]

#subset 
#Reference is G
pos4329162_C1_sub<- pos4329162_C1[,7:ncol(pos4329162_C1)]
#Reference is C
pos4329159_C2_sub<- pos4329159_C2[,7:ncol(pos4329159_C2)]

#get totals and unique numbers 
pos4329162_C1_sub_b<-t(pos4329162_C1_sub)
C1_table<-table(pos4329162_C1_sub_b)
View(C1_table)
C1_table

pos4329159_C2_sub_b<-t(pos4329159_C2_sub)
C2_table<-table(pos4329159_C2_sub_b)
View(C2_table)

pos4329162_C2_sub_c<- as.data.frame(t(pos4329159_C2_sub_b))
#pull out the names of the 11 isolates for Clade 2 
G_to_C_mutations<- pos4329162_C2_sub_c[,pos4329162_C2_sub_c == "T"]








######test code with HRdata#######

#make test df
REF<- c("A", "A", "T", "A")
a<- c("A", "A", "T", "A")
b<- c("G", "G", "G", ".") 
c<- c("A", "A", "T", ".")
d<- c("A", "A", "T", "A")
dftest<- data.frame(t(rbind(REF,a,b,c,d)))


#convert to T/F dataframe
are_any_true<- as.data.frame(lapply(dftest[], `%in%`, dftest$REF))

#Get only collums where some are true 
test<-dftest[,apply(are_any_true,2,function(x) !all(x == TRUE))]
test



