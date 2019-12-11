##This is the main R parcing script for the AF100 project 
##It preforms basic data cleaning, visulization, and statistics.
#The main input files are the output of snpEff (default output = snpEff_genes.txt) 
#and manually input calde designators and a time line of days since isolation 

#load modules
require(data.table)
require(splitstackshape)
library(plyr)
library(dplyr)
library(tidyverse)
library(rlist)

#set dir
setwd("~/Desktop/Afum_AF100_Clades/data/")
options(stringsAsFactors = FALSE)

#load main dat files from snpEff 
snpEff_clade_1<-read.delim("Clade_1.snpEff.matrix.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
snpEff_clade_2<-read.delim("Clade_2.snpEff.matrix.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

###clean the data for easy processing###
#remove the brackets form ALT
snpEff_clade_1$ALT<-gsub("\\[|\\]", "", snpEff_clade_1$ALT)
snpEff_clade_2$ALT<-gsub("\\[|\\]", "", snpEff_clade_2$ALT)
#remove the "/" from all the the SNP calls 
snpEff_clade_1[] <- lapply(snpEff_clade_1, gsub, pattern='/', replacement="")
snpEff_clade_2[] <- lapply(snpEff_clade_2, gsub, pattern='/', replacement="")

#subset: if the ALT = the outgroup (most of the calls), then it's probably not interesting 
clade1_to_OG<-snpEff_clade_1[snpEff_clade_1$ALT != snpEff_clade_1$B8783_CDC.17,]
clade2_to_OG<-snpEff_clade_2[snpEff_clade_2$ALT != snpEff_clade_2$AZTEC_15_181,]

##also probably not interesting if the mutation is only in the outgroup
#subset based on time points (excluding the reference and the outgroup)
timepts_clade1<- as.numeric(c(3,4,5,6,10, 12)) 
timepts_clade2<- as.numeric(c(8,9,10,11,12)) 
#add AF100, _ and * to potential time points for grep search 
timepts_clade1<-paste("AF100.", timepts_clade1, "_*", sep = "", collapse = NULL)
timepts_clade2<-paste("AF100.", timepts_clade2, "_*", sep = "", collapse = NULL)

subsetdf_C1<- data.frame(clade1_to_OG[, grepl(paste(timepts_clade1, collapse="|"),names(clade1_to_OG)) ])
subsetdf_C2<- data.frame(clade2_to_OG[, grepl(paste(timepts_clade2, collapse="|"),names(clade2_to_OG)) ])

#calculate unique calls (excluding the outgroup)
#C1
unique_calls_per_row_1<- t(apply(subsetdf_C1, 1, function(x) unique(x)))
new_ALT_1<- lapply(unique_calls_per_row_1, function(x) paste(x))
new_ALT2_1<- lapply(new_ALT_1, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_1 <- data.frame(matrix(unlist(new_ALT2_1), nrow=length(new_ALT2_1), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_1)<- "NEW_ALT"

#C2
unique_calls_per_row_2<- t(apply(subsetdf_C2, 1, function(x) unique(x)))
new_ALT_2<- lapply(unique_calls_per_row_2, function(x) paste(x))
new_ALT2_2<- lapply(new_ALT_2, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_2 <- data.frame(matrix(unlist(new_ALT2_2), nrow=length(new_ALT2_2), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_2)<- "NEW_ALT"

#add back in the reference and the type
#C1
new_df_1<- cbind(subsetdf_C1, new_ALT2_df_1)
new_df_w_ref_1<- cbind(new_df_1, "REF" = clade1_to_OG$REF, "TYPE" = clade1_to_OG$TYPE, "GENE" = clade1_to_OG$GENE)
#C2
new_df_2<- cbind(subsetdf_C2, new_ALT2_df_2)
new_df_w_ref_2<- cbind(new_df_2, "REF" = clade2_to_OG$REF, "TYPE" = clade2_to_OG$TYPE, "GENE" = clade2_to_OG$GENE)

#remove any rows where NEW_ALT == "." (insertion in OG)
no_miss_calls_in_OG_1<- new_df_w_ref_1[! with(new_df_w_ref_1, NEW_ALT == "."),]
no_miss_calls_in_OG_2<- new_df_w_ref_2[! with(new_df_w_ref_2, NEW_ALT == "."),]

#strip off all of the "." calls from the NEW_ALT collumn so we can compare new alt to the outgroup
no_miss_calls_in_OG_1$NEW_ALT <- sapply(no_miss_calls_in_OG_1$NEW_ALT, function(x) gsub("[.]", "", x))
no_miss_calls_in_OG_2$NEW_ALT <- sapply(no_miss_calls_in_OG_2$NEW_ALT, function(x) gsub("[.]", "", x))

#if the REF col matches the new alt, don't return that line:
cleaned_df_C1<- no_miss_calls_in_OG_1[! with(no_miss_calls_in_OG_1, REF == no_miss_calls_in_OG_1$NEW_ALT),]
cleaned_df_C2<- no_miss_calls_in_OG_2[! with(no_miss_calls_in_OG_2, REF == no_miss_calls_in_OG_2$NEW_ALT),]



###Investigate the data###
##combine to better visualize the differnces 
#make new df's to fill in missing values
df_c1<- data.frame(cbind(rep("Clade_1", nrow(cleaned_df_C1)), cleaned_df_C1$TYPE))
colnames(df_c1)<- "Clade"
df_c2<- data.frame(cbind(rep("Clade_2", nrow(cleaned_df_C2)), cleaned_df_C2$TYPE))
colnames(df_c2)<- "Clade"
#combine them
df_c1_2<- rbind(df_c1, df_c2)
df_c1_2_table<- data.frame(table(df_c1_2))

#get totals for each type of mutation 
mut_type_clade1<- table(cleaned_df_C1$TYPE)
mut_type_clade2<- table(cleaned_df_C2$TYPE)

#look at the clades individually 
barplot(mut_type_clade1[order(mut_type_clade1, decreasing=T)], cex.names = .5)
barplot(mut_type_clade2[order(mut_type_clade2, decreasing=T)], cex.names = .5)

#look at the clades together
barplot(df_c1_2_table$Freq,beside=T,ylim=c(0,11000),
        main="mutation type",
        ylab="frequency",
        axis.lty="solid", 
        col = c("grey", "maroon"), 
        names.arg = df_c1_2_table$NA.,
        axisnames = TRUE, 
        cex.names = .4, 
        las =2)

legend("topright", 
       legend = c("Clade 1", "Clade 2"), 
       fill = c("grey", "maroon"))

#lots of intergenic mutations, and we probably arn't interested in these, likewise synonamous variants
#subset to remove intergenic SNPs and INDELS
#C1
clade1_no_intergenic_a<- cleaned_df_C1[!cleaned_df_C1$TYPE == "intergenic",]
clade1_no_intergenic_b<- clade1_no_intergenic_a[!clade1_no_intergenic_a$TYPE == "intron_variant",]
clade1_no_intergenic<- clade1_no_intergenic_b[!clade1_no_intergenic_b$TYPE == "synonymous_variant",]
#C2
clade2_no_intergenic_a<- cleaned_df_C2[!cleaned_df_C2$TYPE == "intergenic",]
clade2_no_intergenic_b<- clade2_no_intergenic_a[!clade2_no_intergenic_a$TYPE == "intron_variant",]
clade2_no_intergenic<- clade2_no_intergenic_b[!clade2_no_intergenic_b$TYPE == "synonymous_variant",]

#combine the subset clades to better visualize the differnces 
#first malke new df's to fill in missing values. 
df_c1b<- data.frame(cbind(rep("Clade_1", nrow(clade1_no_intergenic)), clade1_no_intergenic$TYPE))
df_c2b<- data.frame(cbind(rep("Clade_2", nrow(clade2_no_intergenic)), clade2_no_intergenic$TYPE))
#combine them
df_c1_2b<- rbind(df_c1b, df_c2b)
df_c1_2_table_b<- data.frame(table(df_c1_2b))

#grouped bar plot without intergenic mutations 
barplot(df_c1_2_table_b$Freq,beside=T,ylim=c(0,1800),
        main="mutation type - no intergenic or syn. variants",
        ylab="frequency",
        axis.lty="solid", 
        col = c("grey", "maroon"), 
        names.arg = df_c1_2_table_b$X2,
        axisnames = TRUE, 
        cex.names = .4, 
        las =2)

legend("topright", 
       legend = c("Clade 1", "Clade 2"), 
       fill = c("grey", "maroon"))


##Missessence mutations likely have the most impact. 
#look only at missence mutations 
missence_clade1<- clade1_no_intergenic[clade1_no_intergenic$TYPE == "missense_variant",]
missence_clade2<- clade2_no_intergenic[clade2_no_intergenic$TYPE == "missense_variant",]

#get a list of the genes inpacted by the varriants and look at overlaps between the two clades
clade1_MGs<- data.frame(missence_clade1$GENE)
clade2_MGs<- data.frame(missence_clade2$GENE)

#remove redundancies
clade1_MGs_unique<- data.frame(unique(clade1_MGs))
clade2_MGs_unique<-data.frame(unique(clade2_MGs))

#get overlap in gene names
genes_in_common<- Reduce(intersect, list(clade1_MGs_unique$missence_clade1.GENE, clade2_MGs_unique$missence_clade2.GENE))

length(genes_in_common)

#print them to txt file - my c/p is not working ... 
write.table(clade1_MGs, "MISS_MUTS_Clade1.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(clade2_MGs, "MISS_MUTS_Clade2.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genes_in_common, "genes_in_common.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



###INVESTIGATE MUTATIONAL ACUMULATION###
#this scrip looks for mutational acumulation over time, while
#subtracting out all variants that apeared in the previous time point
#to report only novel mutations that appeared at that time point

#step 1:
for (i in 1:20){
        T_F_Clade2<- data.frame(synanomous_only_clade2$REF == synanomous_only_clade2[,14: ncol(synanomous_only_clade2) -2])
}

names(synanomous_only_clade1)

###step 1 functionlized: subset df to only include AF100 strains + the outgroup 
#parcing by time point
convert_to_t_f <- function(df, timepts){
        subsetdf<- data.frame(df[, grepl(paste(timepts, collapse="|"),names(df)) ])
        T_F_df <- data.frame(subsetdf$REF == subsetdf[,2: ncol(subsetdf)])
        return(T_F_df)
}

#works:
#subsetdf<- data.frame(clade1_no_intergenic[, grepl(paste(timepts_clade1, collapse="|"),names(clade1_no_intergenic)) ])
T_F_df <- data.frame(subsetdf$REF == subsetdf[,2: ncol(subsetdf)])

dim(subsetdf)
clade1_no_intergenic

##investigate all non-synanamous mutations (whould should be increasing if they're sellected for)
T_F_df_clade1<- convert_to_t_f(df = clade1_no_intergenic, timepts = timepts_clade1)
T_F_df_clade2<- convert_to_t_f(df = clade2_no_intergenic, timepts = timepts_clade2)

##investigate all synanamous mutations (these should be flat if the data is not bised)
T_F_df_clade1<- convert_to_t_f(df = synanomous_only_clade1, timepts = timepts_clade1)
T_F_df_clade2<- convert_to_t_f(df = synanomous_only_clade2, timepts = timepts_clade2)


##step 2 functionalize: #(do any of the time points in that time point match the reference?)
#split the dataframes based on timepoints
#list_of_dfs_by_timept9<- lapply(setNames(timepts, timepts), function(x) data.frame(T_F_Clade2[, grep(x, colnames(T_F_Clade2))]))

are_any_true <- function(df, timepts){
        T_F_df_list<- 0
        #remove reference col
        timepts2<- timepts[-1]
        #split dfs
        split_dfs<- lapply(setNames(timepts2, timepts2), function(x) data.frame(df[, grep(x, colnames(df))]))
        #for each position, are there any for that time pt taht are true?
        for (i in 1:length(split_dfs)){
                T_F_df_list[i]<- data.frame(apply(split_dfs[[i]], 1, any))
        }
        names(T_F_df_list)<- timepts2
        T_F_df_list_df<- as.data.frame(lapply(T_F_df_list, cbind))
        return(T_F_df_list_df)
}
#run function
clade_1_TF_by_time<- are_any_true(df = T_F_df_clade1, timepts = timepts_clade1)
View(clade_1_TF_by_time)



##step 3: functionalize: is it the first instance of a new mutation?
are_any_true <- function(df){
        status<- 0
        #attach a null collumn for the first time point calculation 
        df_w_null<-cbind(logical(nrow(df)), df)
        #run loop
        for (i in 2:length(df_w_null)){
                status[i]<- data.frame(apply((!df_w_null[,i, drop=FALSE] == df_w_null[,i-1, drop=FALSE]), 1, isTRUE))
                bstatuschangedf<- data.frame(lapply(status[2:length(status)], cbind))
                return(bstatuschangedf)
                #rename
        }
        
}



#test function
status_change<- are_any_true(df = clade_1_TF_by_time)
View(status_change)

#try piecemeal:
status<- 0
#attach a null collumn for the first time point calculation 
df_w_null<-cbind(logical(nrow(clade_1_TF_by_time)), clade_1_TF_by_time)
#run loop
bstatuschangedf<- for (i in 2:length(df_w_null)){
        status[i]<- data.frame(apply((!df_w_null[,i, drop=FALSE] == df_w_null[,i-1, drop=FALSE]), 1, isTRUE))
        #bstatuschangedf<- data.frame(lapply(status[2:length(status)], cbind))
        #return(status)
        #rename
}

for (i in 2:length(are_any_true_df_df_df2)){
        status[i]<- data.frame(apply((!are_any_true_df_df_df2[,i, drop=FALSE] == are_any_true_df_df_df2[,i-1, drop=FALSE]), 1, isTRUE))
        bstatuschangedf<- data.frame(lapply(status[2:length(status)], cbind))
        #rename
}        

View(status)

##you are here - fix from here. .. . . 
View(bstatuschangedf)





dim(status_change)



#attach a null collumn for the first time point calculation 
are_any_true_df_df_df2<-cbind(logical(nrow(are_any_true_df_df2)), are_any_true_df_df2)







#plot the new changes over time
totals_per_timept2<- colSums(bstatuschangedf)
#attach time pts. 
dsi2<- c(968, 1145, 1221, 1346, 1498)
totals_per_timept2<- data.frame(cbind(dsi2,totals_per_timept2))

#NO reapply time pts. 
#graph here Clade2 novel mutations here
Clade2_new_muts_over_time<- plot(totals_per_timept2$dsi, totals_per_timept2$totals_per_timept,
                                 col = "red", 
                                 xlab = "days after isolation", 
                                 ylab = "total new mutations",
                                 main = "Clade 1 Synanonamous Mutations", 
                                 pch = 16)

clade2_new_lm<- lm(totals_per_timept2$totals_per_timept ~ totals_per_timept2$dsi)
sum1<- summary(clade2_new_lm)
abline(clade2_new_lm)
rsq_clade2_new <- bquote(italic(R)^2 == .(format(summary(clade2_new_lm)$adj.r.squared, digits=2)))
text(x = 1100, y = 10, labels = rsq_clade2_new)
pval_calde2_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde2_new <- bquote(italic(P) == .(pval_calde2_new))
text(x = 1100, y = 11, labels = pval_calde2_new)






















#######OLD SCRIPTS - NON FUNCTIONALIZED BELOW:








#make summary of how many mutations there are per strain (different from the reference)
#pseudo code
#for each column that contains a speices in clad<>_no_intergenic, give me the sum of calls that are different than the reference and not "."

#this works on a single entry
#total_muts_per_strain<- sum(clad2_no_intergenic$REF == clad2_no_intergenic$X1F1SW_F4)
#total_muts_per_strain

#loop over all of them
for (i in 1:20){
        muts<- (clad1_no_intergenic$REF == clad1_no_intergenic[,13: ncol(clad1_no_intergenic) -1])
        total_muts_clade1<- data.frame((colSums(muts)))
}

for (i in 1:13){
        muts<- (clad2_no_intergenic$REF == clad2_no_intergenic[,14: ncol(clad2_no_intergenic) -2])
        total_muts_clade2<- data.frame(print(colSums(muts)))
}

#attach days since isolation 
#read in date file
Clade1_mut_totals<-read.csv("Clade1_mut_totals.csv", header = TRUE, fill = TRUE, strip.white = TRUE)
Clade2_mut_totals<-read.csv("Clade2_mut_totals.csv", header = TRUE, fill = TRUE, strip.white = TRUE)

#direct graph output to pdf
pdf("mut_over_time_AF100.pdf", 5, 7)
par(mfrow=c(2,1))

#take a look at how mutations are acumulation over time
Clade1plot<- plot(Clade1_mut_totals$dsi, Clade1_mut_totals$total_mut,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total mutations",
                  main = "Clade 1", 
                  pch = 16)

clade1_lm<- lm(Clade1_mut_totals$total_mut ~ Clade1_mut_totals$dsi)
sum1<- summary(clade1_lm)
abline(clade1_lm)
#legend("bottomleft", bty="n", legend=paste("R2=", format(summary(clade1_lm)$adj.r.squared, digits=4)), "\t yes")

rsq_clade1 <- bquote(italic(R)^2 == .(format(summary(clade1_lm)$adj.r.squared, digits=3)))
text(x = 1200, y = 208, labels = rsq_clade1)
pval_calde1<- round(sum1$coefficients[2,4], digits = 5)
pval_clade1 <- bquote(italic(P) == .(pval_calde1))
text(x = 1200, y = 215, labels = pval_clade1)



#is it normal?
#qqnorm(residuals(clade1_lm), main="", datax=TRUE)
#qqline(residuals(clade1_lm), datax=TRUE)
#qqplot(Clade1_mut_totals$total_mut, Clade1_mut_totals$dsi)


#take a look at how mutations are acumulation over time
Clade2plot<- plot(Clade2_mut_totals$dsi, Clade2_mut_totals$total_mut,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total mutations",
                  main = "Clade 2", 
                  pch = 16)

clade2_lm<- lm(Clade2_mut_totals$total_mut ~ Clade2_mut_totals$dsi)
sum2<- summary(clade2_lm)
abline(clade2_lm)
#legend("bottomleft", bty="n", legend=paste("R2=", format(summary(clade1_lm)$adj.r.squared, digits=4)), "\t yes")

rsq_clade2 <- bquote(italic(R)^2 == .(format(summary(clade2_lm)$adj.r.squared, digits=3)))
text(x = 1400, y = 450, labels = rsq_clade2)
pval_calde2<- round(sum2$coefficients[2,4], digits = 4)
pval_clade2 <- bquote(italic(P) == .(pval_calde2))
text(x = 1400, y = 475, labels = pval_clade2)

dev.off()




###remake with just synanamous mutations 
synanomous_only_clade1<- clade1_to_OGn[clade1_to_OGn$TYPE == "synonymous_variant",]
synanomous_only_clade2<- clade2_to_OGn[clade2_to_OGn$TYPE == "synonymous_variant",]

#get total n mutations for ea strain
#loop over all of them
for (i in 1:20){
        muts<- (synanomous_only_clade1$REF == synanomous_only_clade1[,13: ncol(synanomous_only_clade1) -1])
        total_syn_muts_clade1<- data.frame((colSums(muts)))
}

for (i in 1:13){
        muts<- (synanomous_only_clade2$REF == synanomous_only_clade2[,14: ncol(synanomous_only_clade2) -2])
        total_syn_muts_clade2<- data.frame(print(colSums(muts)))
}

#print it to attach dates
#read in date file
Clade1_syn_mut_totals<-read.csv("Clade1_syn_mut_totals.csv", header = TRUE, fill = TRUE, strip.white = TRUE)
Clade2_syn_mut_totals<-read.csv("Clade2_syn_mut_totals.csv", header = TRUE, fill = TRUE, strip.white = TRUE)

##plot syn. mutations for both clades
#send to pdf output
pdf("syn_mut_over_time_AF100.pdf", 5, 8)
par(mfrow=c(2,1))
#take a look at how mutations are acumulation over time
Clade1plot<- plot(Clade1_syn_mut_totals$dsi, Clade1_syn_mut_totals$total_mut,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total syn mutations",
                  main = "Clade 1 synonymous", 
                  pch = 16)

clade1_lm<- lm(Clade1_syn_mut_totals$total_mut ~ Clade1_syn_mut_totals$dsi)
sum1<- summary(clade1_lm)
abline(clade1_lm)
rsq_clade1 <- bquote(italic(R)^2 == .(format(summary(clade1_lm)$adj.r.squared, digits=2)))
text(x = 1220, y = 8.5, labels = rsq_clade1)
pval_calde1<- round(sum1$coefficients[2,4], digits = 3)
pval_clade1 <- bquote(italic(P) == .(pval_calde1))
text(x = 1200, y = 9, labels = pval_clade1)

#for clade 2
Clade2plot<- plot(Clade2_syn_mut_totals$dsi, Clade2_syn_mut_totals$total_mut,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total syn mutations",
                  main = "Clade 2 synonymous", 
                  pch = 16)

clade2_lm<- lm(Clade2_syn_mut_totals$total_mut ~ Clade2_syn_mut_totals$dsi)
sum2<- summary(clade2_lm)
abline(clade2_lm)
rsq_clade1 <- bquote(italic(R)^2 == .(format(summary(clade2_lm)$adj.r.squared, digits=3)))
text(x = 1400, y = 85, labels = rsq_clade1)
pval_calde2<- round(sum1$coefficients[2,4], digits = 3)
pval_clade2 <- bquote(italic(P) == .(pval_calde2))
text(x = 1400, y = 92, labels = pval_clade2)

dev.off()


### plots with all mutations remake with all mutations
#get total n mutations for ea strain
#loop over all of them
for (i in 1:20){
        muts<- (clade1_to_OGn$REF == clade1_to_OGn[,13: ncol(clade1_to_OGn) -1])
        total_muts_clade1_all<- data.frame((colSums(muts)))
}

for (i in 1:13){
        muts<- (clade2_to_OGn$REF == clade2_to_OGn[,14: ncol(clade2_to_OGn) -2])
        total_muts_clade2_all<- data.frame(print(colSums(muts)))
}

#read int he data files:
Clade1_mut_totals_all<-read.csv("Clade1_mut_totals_all.csv", header = TRUE, fill = TRUE, strip.white = TRUE)
Clade2_mut_totals_all<-read.csv("Clade2_mut_totals_all.csv", header = TRUE, fill = TRUE, strip.white = TRUE)

##plot syn. mutations for both clades
#send to pdf output
pdf("all_mut_over_time_AF100.pdf", 5, 8)
par(mfrow=c(2,1))
#take a look at how mutations are acumulation over time
Clade1plot<- plot(Clade1_mut_totals_all$dsi, Clade1_mut_totals_all$total_mut,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total n mutations",
                  main = "Clade 1 all mutations", 
                  pch = 16)

clade1_lm<- lm(Clade1_mut_totals_all$total_mut ~ Clade1_mut_totals_all$dsi)
sum1<- summary(clade1_lm)
abline(clade1_lm)
rsq_clade1 <- bquote(italic(R)^2 == .(format(summary(clade1_lm)$adj.r.squared, digits=2)))
text(x = 1220, y = 620, labels = rsq_clade1)
pval_calde1<- round(sum1$coefficients[2,4], digits = 7)
pval_clade1 <- bquote(italic(P) == .(pval_calde1))
text(x = 1230, y = 700, labels = pval_clade1)

#and clade 2
Clade2plot<- plot(Clade2_mut_totals_all$dsi, Clade2_mut_totals_all$total_mut,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total n mutations",
                  main = "Clade 2 all mutations", 
                  pch = 16)

clade2_lm<- lm(Clade2_mut_totals_all$total_mut ~ Clade2_mut_totals_all$dsi)
sum1<- summary(clade2_lm)
abline(clade2_lm)
rsq_clade2 <- bquote(italic(R)^2 == .(format(summary(clade2_lm)$adj.r.squared, digits=2)))
text(x = 1400, y = 2000, labels = rsq_clade2)
pval_calde2<- round(sum1$coefficients[2,4], digits = 7)
pval_clade2 <- bquote(italic(P) == .(pval_calde2))
text(x = 1420, y = 2800, labels = pval_clade2)

dev.off()





###calculate only new mutations over time###

###one thing to think about here is whether we need to change the "." into TRUE first (convert to NA or to same as reference call before running t/f matrix)

###to identify new mutations: (for each row in are_any_true_df), for each column, if FLASE pass
                                                                #if TRUE, is the previous col also TRUE?
                                                                #if previous collumnn is TRUE print FLASE
                                                                #if previous column is FALSE, print TRUE



#for clade 1:
##REAL: step 1 (does it match the reference?) 
for (i in 1:20){
        T_F_Clade1<- data.frame(clade1_to_OGn$REF == clade1_to_OGn[,13: ncol(clade1_to_OGn) -1])
}

##REAL step 2 (do and any of the time points in that time point match the reference?)
#split the dataframes by their timepoint 
#create a list of potential time points
#timepts<- as.numeric(2:12) - all - set specific for the time pts present in clade 1
timepts<- as.numeric(c(3,4,5,6,10,12)) 
#must change this for calde2

#add AF100, _ and * to potential time points for grep search 
timepts<-paste("AF100.", timepts, "_*", sep = "", collapse = NULL)
#split the dataframes based on timepoints
list_of_dfs_by_timept<- lapply(setNames(timepts, timepts), function(x) T_F_Clade1[, grep(x, colnames(T_F_Clade1))])

#for each position, are any true in that time pt? 
#obstantiate
are_any_true_df<- 1
#run the loop
for (i in 1:length(list_of_dfs_by_timept)){
        are_any_true_df[i]<- as.data.frame(apply(list_of_dfs_by_timept[[i]], 1, any))
}
#attach names
names(are_any_true_df)<- timepts

#reformat as dataframe
are_any_true_df_df<- as.data.frame(lapply(are_any_true_df, cbind))

##step 3 (is it the first instance of a new mutation?)

#attach a null collumn for the first time point calculation 
are_any_true_df_df_df<-cbind(logical(nrow(are_any_true_df_df)), are_any_true_df_df)
#obstantiate:
status<- 1
#run loop
for (i in 2:length(are_any_true_df_df_df)){
        status[i]<- as.data.frame(apply((!are_any_true_df_df_df[,i, drop=FALSE] == are_any_true_df_df_df[,i-1, drop=FALSE]), 1, isTRUE))
        bstatuschangedf<- data.frame(lapply(status[2:length(status)], cbind))

}
names(bstatuschangedf)<- timepts



#plot the new changes over time
totals_per_timept<- colSums(bstatuschangedf)
#attach time pts. 
dsi<- c(135, 249, 471, 624, 1221, 1498)
totals_per_timept<- data.frame(cbind(dsi,totals_per_timept))

Clade1_new_muts_over_time<- plot(totals_per_timept$dsi, totals_per_timept$totals_per_timept,
                  col = "red", 
                  xlab = "days after isolation", 
                  ylab = "total new mutations",
                  main = "Clade 1 New Mutations", 
                  pch = 16)

clade1_new_lm<- lm(totals_per_timept$totals_per_timept ~ totals_per_timept$dsi)
sum1<- summary(clade1_new_lm)
abline(clade1_new_lm)
rsq_clade1_new <- bquote(italic(R)^2 == .(format(summary(clade1_new_lm)$adj.r.squared, digits=2)))
text(x = 1300, y = 1100, labels = rsq_clade1_new)
pval_calde1_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde1_new <- bquote(italic(P) == .(pval_calde1_new))
text(x = 1300, y = 1000, labels = pval_calde1_new)

#ave ave mut per year ###GO BAKC HERE
perday<-  totals_per_timept$dsi / totals_per_timept$totals_per_timept
aveperday<- ave(perday)


#####CALDE 2
for (i in 1:20){
        T_F_Clade2<- data.frame(clade2_to_OGn$REF == clade2_to_OGn[,13: ncol(clade2_to_OGn) -1])
}


##step 2 (do and any of the time points in that time point match the reference?)
#split the dataframes by their timepoint 
#create a list of potential time points
#timepts<- as.numeric(2:12) - all - set specific for the time pts present in clade 1
timepts<- as.numeric(c(8,9,10,11,12)) 
#must change this for calde2

#add AF100, _ and * to potential time points for grep search 
timepts<-paste("AF100.", timepts, "_*", sep = "", collapse = NULL)
#split the dataframes based on timepoints
list_of_dfs_by_timept2<- lapply(setNames(timepts, timepts), function(x) data.frame(T_F_Clade2[, grep(x, colnames(T_F_Clade2))]))

#for each position, are any true in that time pt? 

#obstantiate
are_any_true_df<- 1

#run the loop
for (i in 1:length(list_of_dfs_by_timept2)){
        are_any_true_df2[i]<- data.frame(apply(list_of_dfs_by_timept2[[i]], 1, any))
}

#attach names
names(are_any_true_df2)<- timepts

#reformat as dataframe
are_any_true_df_df2<- as.data.frame(lapply(are_any_true_df2, cbind))

##step 3 (is it the first instance of a new mutation?)

#attach a null collumn for the first time point calculation 
are_any_true_df_df_df2<-cbind(logical(nrow(are_any_true_df_df2)), are_any_true_df_df2)

#obstantiate:
status<- 1
#run loop
for (i in 2:length(are_any_true_df_df_df2)){
        status[i]<- as.data.frame(apply((!are_any_true_df_df_df2[,i, drop=FALSE] == are_any_true_df_df_df2[,i-1, drop=FALSE]), 1, isTRUE))
        bstatuschangedf<- as.data.frame(lapply(status[2:length(status)], cbind))
        #rename
        }
        names(bstatuschangedf)<- timepts

#plot the new changes over time
totals_per_timept2<- colSums(bstatuschangedf)
#attach time pts. 
dsi2<- c(968, 1145, 1221, 1346, 1498)
totals_per_timept2<- data.frame(cbind(dsi2,totals_per_timept2))

#NO reapply time pts. 
#graph here Clade2 novel mutations here
Clade2_new_muts_over_time<- plot(totals_per_timept2$dsi, totals_per_timept2$totals_per_timept,
                                 col = "red", 
                                 xlab = "days after isolation", 
                                 ylab = "total new mutations",
                                 main = "Clade 2 New Mutations", 
                                 pch = 16)

clade2_new_lm<- lm(totals_per_timept2$totals_per_timept ~ totals_per_timept2$dsi)
sum1<- summary(clade2_new_lm)
abline(clade2_new_lm)
rsq_clade2_new <- bquote(italic(R)^2 == .(format(summary(clade2_new_lm)$adj.r.squared, digits=2)))
text(x = 1400, y = 2000, labels = rsq_clade2_new)
pval_calde2_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde2_new <- bquote(italic(P) == .(pval_calde2_new))
text(x = 1400, y = 2300, labels = pval_calde2_new)




####remake with synanomous mutations only
##clade 1: 
#####CALDE 1

for (i in 1:20){
        T_F_Clade2<- data.frame(synanomous_only_clade1$REF == synanomous_only_clade1[,13: ncol(synanomous_only_clade1) -1])
}


##step 2 (do and any of the time points in that time point match the reference?)
#split the dataframes by their timepoint 
#create a list of potential time points
#timepts<- as.numeric(2:12) - all - set specific for the time pts present in clade 1
timepts<- as.numeric(c(3,4,5,6,10,12)) 
#must change this for calde2

#add AF100, _ and * to potential time points for grep search 
timepts<-paste("AF100.", timepts, "_*", sep = "", collapse = NULL)
#split the dataframes based on timepoints
list_of_dfs_by_timept2<- lapply(setNames(timepts, timepts), function(x) data.frame(T_F_Clade2[, grep(x, colnames(T_F_Clade2))]))

#for each position, are any true in that time pt? 

#obstantiate
are_any_true_df<- 1

#run the loop
for (i in 1:length(list_of_dfs_by_timept2)){
        are_any_true_df2[i]<- data.frame(apply(list_of_dfs_by_timept2[[i]], 1, any))
}

#attach names
names(are_any_true_df2)<- timepts

#reformat as dataframe
are_any_true_df_df2<- as.data.frame(lapply(are_any_true_df2, cbind))

##step 3 (is it the first instance of a new mutation?)

#attach a null collumn for the first time point calculation 
are_any_true_df_df_df2<-cbind(logical(nrow(are_any_true_df_df2)), are_any_true_df_df2)

#obstantiate:
status<- 1
#run loop
for (i in 2:length(are_any_true_df_df_df2)){
        status[i]<- as.data.frame(apply((!are_any_true_df_df_df2[,i, drop=FALSE] == are_any_true_df_df_df2[,i-1, drop=FALSE]), 1, isTRUE))
        bstatuschangedf<- as.data.frame(lapply(status[2:length(status)], cbind))
        #rename
}
names(bstatuschangedf)<- timepts

#plot the new changes over time
totals_per_timept2<- colSums(bstatuschangedf)
#attach time pts. 
dsi2<- c(135, 249, 471, 624, 1221, 1498)
totals_per_timept2<- data.frame(cbind(dsi2,totals_per_timept2))

#NO reapply time pts. 
#graph here Clade2 novel mutations here
Clade2_new_muts_over_time<- plot(totals_per_timept2$dsi, totals_per_timept2$totals_per_timept,
                                 col = "red", 
                                 xlab = "days after isolation", 
                                 ylab = "total new mutations",
                                 main = "Clade 1 Synanonamous Mutations", 
                                 pch = 16)

clade2_new_lm<- lm(totals_per_timept2$totals_per_timept ~ totals_per_timept2$dsi)
sum1<- summary(clade2_new_lm)
abline(clade2_new_lm)
rsq_clade2_new <- bquote(italic(R)^2 == .(format(summary(clade2_new_lm)$adj.r.squared, digits=2)))
text(x = 1100, y = 10, labels = rsq_clade2_new)
pval_calde2_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde2_new <- bquote(italic(P) == .(pval_calde2_new))
text(x = 1100, y = 11, labels = pval_calde2_new)

#####CALDE 2

for (i in 1:20){
        T_F_Clade2<- data.frame(synanomous_only_clade2$REF == synanomous_only_clade2[,14: ncol(synanomous_only_clade2) -2])
}

##step 2 (do and any of the time points in that time point match the reference?)
#split the dataframes by their timepoint 
#create a list of potential time points
#timepts<- as.numeric(2:12) - all - set specific for the time pts present in clade 1
timepts<- as.numeric(c(8,9,10,11,12)) 
#must change this for calde2

#add AF100, _ and * to potential time points for grep search 
timepts<-paste("AF100.", timepts, "_*", sep = "", collapse = NULL)
#split the dataframes based on timepoints
list_of_dfs_by_timept9<- lapply(setNames(timepts, timepts), function(x) data.frame(T_F_Clade2[, grep(x, colnames(T_F_Clade2))]))
#for each position, are any true in that time pt? 


#obstantiate
are_any_true_df6<- 0

#run the loop
for (i in 1:length(list_of_dfs_by_timept9)){
        are_any_true_df6[i]<- data.frame(apply(list_of_dfs_by_timept9[[i]], 1, any))
}

View(are_any_true_df6)
#attach names
names(are_any_true_df6)<- timepts

#reformat as dataframe
are_any_true_df_df2<- as.data.frame(lapply(are_any_true_df6, cbind))

##step 3 (is it the first instance of a new mutation?)

#attach a null collumn for the first time point calculation 
are_any_true_df_df_df2<-cbind(logical(nrow(are_any_true_df_df2)), are_any_true_df_df2)

#obstantiate:
status<- 1
#run loop
for (i in 2:length(are_any_true_df_df_df2)){
        status[i]<- as.data.frame(apply((!are_any_true_df_df_df2[,i, drop=FALSE] == are_any_true_df_df_df2[,i-1, drop=FALSE]), 1, isTRUE))
        bstatuschangedf<- as.data.frame(lapply(status[2:length(status)], cbind))
        #rename
}
names(bstatuschangedf)<- timepts

#plot the new changes over time
totals_per_timept2<- colSums(bstatuschangedf)
#attach time pts. 
dsi2<- c(968, 1145, 1221, 1346, 1498)
totals_per_timept2<- data.frame(cbind(dsi2,totals_per_timept2))

#NO reapply time pts. 
#graph here Clade2 novel mutations here
Clade2_new_muts_over_time<- plot(totals_per_timept2$dsi, totals_per_timept2$totals_per_timept,
                                 col = "red", 
                                 xlab = "days after isolation", 
                                 ylab = "total new mutations",
                                 main = "Clade 1 Synanonamous Mutations", 
                                 pch = 16)

clade2_new_lm<- lm(totals_per_timept2$totals_per_timept ~ totals_per_timept2$dsi)
sum1<- summary(clade2_new_lm)
abline(clade2_new_lm)
rsq_clade2_new <- bquote(italic(R)^2 == .(format(summary(clade2_new_lm)$adj.r.squared, digits=2)))
text(x = 1100, y = 10, labels = rsq_clade2_new)
pval_calde2_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde2_new <- bquote(italic(P) == .(pval_calde2_new))
text(x = 1100, y = 11, labels = pval_calde2_new)


###rewrite as functions









######
##test step 1 on mock data: 
#make mock data
a<- c("A", "A", "T", "A")
b<- c("G", "G", "G", ".") 
c<- c("T", "A", "A", "T")
dftest<- data.frame(rbind(a,b,c))
colnames(dftest) <- c("REF", "OG", "SOI1", "SOI2")
timepts_test<- c("REF", "OG", "SOI1", "SOI2")  

#run the function
T_F_df_test<- convert_to_t_f(df = dftest, timepts = timepts_test)
T_F_df_test
#OK that works. . . 



##test step 2 on mock data:
#function
are_any_true <- function(df, timepts){
        T_F_df_list<- 0
        #remove reference col
        timepts2<- timepts[-1]
        #split dfs
        split_dfs<- lapply(setNames(timepts2, timepts2), function(x) data.frame(df[, grep(x, colnames(df))]))
        #for each position, are there any for that time pt taht are true?
        for (i in 1:length(split_dfs)){
                T_F_df_list[i]<- data.frame(apply(split_dfs[[i]], 1, any))
        }
        names(T_F_df_list)<- timepts2
        T_F_df_list_df<- as.data.frame(lapply(T_F_df_list, cbind))
        return(T_F_df_list_df)
}

#run function
test_are_any_true<- are_any_true(df = T_F_df_test, timepts = timepts_test)
View(test_are_any_true)
