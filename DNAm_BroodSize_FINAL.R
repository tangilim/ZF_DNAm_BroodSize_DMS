---
title: "Long-term effects of early-life adversity on DNA methylation in zebra finches"
author: "Marianthi Tangili"
date: "30/03/2025"
---

#script for the identification of differentially methylated CpG sites (DMSs) from whole genome EMSeq zebra finch data

#raw data can be found in https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1108628

#bioinformatics pipeline from raw sequencing data to .cov.gz methylation files can be found in https://zenodo.org/records/12074859



#load packages
library(data.table)
library(methylKit)
library(ggplot2)
library(readr)
library(purrr)
library(janitor)
library(dplyr)
library(tidyr)
library(arsenal)
library(readxl)

#load graph theme
graph_theme<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour =NA ,fill=NA,size=0.5),
  axis.title.x=element_text(size=15,hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=15,hjust=0.5,vjust=0.5,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=10),
  axis.text.y=element_text(colour="black",angle=0,size=10),
  axis.ticks=element_line(colour="black",size=0.5),
  axis.line.x=element_line(size=0.5),
  axis.line.y=element_line(size=0.5))



setwd("SOMEPATH")


#differential methylation calculation

#2 samples from each individual are merged and methylation is calculated

#define BirdID and Sample pairs
#BirdID= individual
#Sample= sample ID
bird_data <- data.frame(
  BirdID = c(1810, 1810, 1825, 1825, 1897, 1897, 1904, 1904, 1907, 1907, 1920, 1920, 1923, 1923, 1942, 1942, 1944, 1944, 1953, 1953, 1983, 1983, 2035, 2035, 2038, 2038, 2045, 2045, 2053, 2053, 2058, 2058, 2069, 2069, 2087, 2087, 2108, 2108, 2116, 2116, 2124, 2124, 2142, 2142, 2153, 2153, 2155, 2155, 2163, 2163, 2172, 2172, 2175, 2175, 2229, 2229, 2236, 2236, 2247, 2247, 2257, 2257, 2259, 2259, 2260, 2260, 2269, 2269, 2277, 2277, 2288, 2288, 2291, 2291, 2354, 2354, 2390, 2390, 2392, 2392, 2415, 2415, 2418, 2418, 2423, 2423, 2429, 2429, 2463, 2463, 2470, 2470, 2473, 2473, 2487, 2487, 2508, 2508, 2511, 2511),
  Sample = c(1, 52, 2, 53, 63, 76, 3, 83, 4, 70, 5, 50, 29, 78, 30, 84, 6, 75, 7, 45, 8, 26, 74, 88, 65, 82, 9, 38, 10, 54, 11, 31, 12, 55, 13, 33, 14, 77, 59, 89, 15, 69, 16, 86, 17, 56, 18, 51, 19, 46, 34, 79, 20, 90, 21, 66, 24, 97, 47, 92, 39, 62, 25, 60, 22, 48, 32, 94, 40, 67, 35, 96, 23, 41, 68, 98, 43, 95, 81, 100, 80, 93, 42, 99, 36, 61, 27, 72, 49, 64, 85, 91, 28, 73, 37, 87, 57, 58, 44, 71)
)

#function to process and merge samples by BirdID
process_and_merge_samples <- function(bird_id, samples) {
  #read and process each sample
  meth_data_list <- lapply(samples, function(sample_id) {
    meth_data <- fread(paste0("S", sample_id, ".deduplicated.bismark.cov.gz"), header = FALSE)
    setnames(meth_data, c("chr", "start", "end", "methylation", "numCs", "numTs"))
    return(meth_data)
  })
  
  #merge the data
  merged_data <- Reduce(function(x, y) {
    merge(x, y, by = c("chr", "start", "end"), suffixes = c(".1", ".2"))
  }, meth_data_list)
  
  #calculate combined numCs and coverage
  merged_data[, `:=`(
    numCs = numCs.1 + numCs.2,
    coverage = (numCs.1 + numTs.1) + (numCs.2 + numTs.2)
  )]
  
  #create methylRaw object
  methyl_obj <- new("methylRaw", 
                    data.frame(chr = merged_data$chr,
                               start = merged_data$start,
                               end = merged_data$end,
                               strand = "*",
                               coverage = merged_data$coverage,
                               numCs = merged_data$numCs,
                               numTs = merged_data$coverage - merged_data$numCs),
                    sample.id = paste0("Bird", bird_id),
                    assembly = "tgut1.4",
                    context = "CpG",
                    resolution = "base")
  
  return(methyl_obj)
}

#process and merge samples for each BirdID
bird_ids <- unique(bird_data$BirdID)
methyl_list <- lapply(bird_ids, function(bird_id) {
  samples <- bird_data$Sample[bird_data$BirdID == bird_id]
  process_and_merge_samples(bird_id, samples)
})

#name the list elements
names(methyl_list) <- paste0("Bird", bird_ids)

#create the treatment vector (1==broods of 2, 2==broods of 6)
treatment <- c(1, 1, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 1)

#create a methylRawList object with the treatment vector
methyl_raw_list <- new("methylRawList", methyl_list, treatment = treatment)

#define coverage filtering thresholds
min_coverage <- 12
max_coverage <- 99.9

#filter by coverage
methyl_filtered <- filterByCoverage(methyl_raw_list,
                                    lo.count = min_coverage,
                                    lo.perc = NULL,
                                    hi.count = max_coverage ,
                                    hi.perc = NULL)

#normalize coverage
methyl_normalized <- normalizeCoverage(methyl_filtered)

#unite the samples (positions present in >=70% of the samples)
meth_united <- unite(methyl_normalized)

meth_percentages <- percMethylation(meth_united)

#Identify positions to keep (not 0% or 100% methylation in all samples)
positions_to_keep <- apply(meth_percentages, 1, function(x) !all(x %in% c(0, 100)))

#filter
meth_united_filtered <- meth_united[positions_to_keep,]

#differential methylation 
myDiff <- calculateDiffMeth(meth_united_filtered, overdispersion=c("none"), adjust=c("bonferroni"))

head(myDiff)

# Figure 1 ----------------------------------------------------------------
#Volcano plot

#Bonferroni-corrected significance threshold
-log10(0.05/nrow(myDiff))
#7.613336

#create new data frame
meth.diff=myDiff$meth.diff
chr=myDiff$chr
pos=myDiff$start
sig=-log10(myDiff$pvalue)
df<-data.frame(fc,sig, chr, pos)
df$thre<-as.factor(abs(fc) >= 25 & sig >  7.613336)

#find significant threshold
significance_threshold <- -log10(0.05 / nrow(myDiff))

#define significant points based on criteria
significant <- myDiff %>%
  filter(abs(meth.diff) >= 25 & -log10(pvalue) >= significance_threshold) 

#number of DMS
nrow(significant)
#149

VolcanoPlot<- ggplot(data=df, aes(x=meth.diff, y=sig, color=thre)) +
  geom_point(alpha=0.8, size=1) +
  theme(legend.position = "none") +
  labs(y=expression("-log"^"10"*"(p-value)"), x = expression(Delta*" % DNA methylation level")) +
  scale_x_continuous(limits=c(-40,40), breaks=c(-35,-25,-15,0,15, 25,35)) +
  theme(panel.background = element_rect(fill="white", colour="white", size=0.5, linetype="solid",
                                        color="grey70")) +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75))+
  scale_color_manual(values = c("azure3","blue")) +
  geom_hline(yintercept =  7.613336,color="black",linetype="dashed",alpha=0.7) +
  geom_vline(xintercept = -25,color="black",linetype="dashed",alpha=0.7) +
  geom_vline(xintercept = 25,color="black",linetype="dashed",alpha=0.7)

png("Fig1_ZF_volcano_treatment_FINAL.png",  res = 300, width = 2000, height = 1200)
VolcanoPlot
dev.off()


# Figure 2 ----------------------------------------------------------------
#Manhattan plot
myDiff$Chromosome <- factor(myDiff$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32','33', '34', '35', '37','Z'))

# Adds CHR_State in order to have alternate colours ~
myDiff_CHR_IDs <- as.data.frame(unique(myDiff$Chromosome)); colnames(myDiff_CHR_IDs) <- c("Chromosome")
myDiff_CHR_IDs$CHR_IDs <- seq.int(nrow(myDiff_CHR_IDs))
myDiffUp <- merge(myDiff, myDiff_CHR_IDs, by = "Chromosome")
myDiffUp <- myDiffUp %>% arrange(CHR_IDs)
myDiffUp$CHR_State <- ifelse(myDiffUp$CHR_IDs %% 2 == 0, "Even", "Odd")


#calculates the max start position for each chromosome ~
cum_start <- myDiffUp %>% 
  group_by(Chromosome) %>% 
  summarise(max_pos = max(start)) %>% 
  mutate(cum_pos = cumsum(max_pos) - max_pos)


#merges the cumulative start positions back to the main data frame ~
myDiffUp <- myDiffUp %>% 
  left_join(cum_start, by = "Chromosome") %>%
  mutate(cum_start = start + cum_pos)


#calculates chromosome midpoints for x-axis labeling ~
chr_midpoints <- cum_start %>% 
  mutate(midpoint = cum_pos + max_pos / 2) %>% 
  select(Chromosome, midpoint)

significance= -log10(0.05/nrow(myDiff))

significant <- myDiffUp$pvalue < 0.05 & abs(myDiffUp$meth.diff) >= 25


#label Chromosomes in custom order
custom_chr_scale <- function(data) {
  # Get unique chromosomes and their maximum positions
  chr_summary <- data.frame(Chromosome = levels(data$Chromosome))
  chr_summary$end <- sapply(chr_summary$Chromosome, function(chr) {
    max_val <- max(data$end[data$Chromosome == chr], na.rm = TRUE)
    if (is.finite(max_val)) max_val else NA
  })
  
  #remove chromosomes with no valid data
  chr_summary <- chr_summary[!is.na(chr_summary$end), ]
  
  #calculate cumulative positions
  chr_summary$cum_pos <- cumsum(chr_summary$end)
  chr_summary$midpoint <- chr_summary$cum_pos - (chr_summary$end / 2)
  
  #create labels
  chr_summary$label <- as.character(chr_summary$Chromosome)
  
  #define chromosomes to be labeled
  standard_chrs <- c(as.character(1:15), as.character(17:20), as.character(seq(22, 28, by=2)))
  special_chrs <- c("1A", "4A", "37", "Z")
  labeled_chrs <- c(standard_chrs, special_chrs)
  
  #set empty labels for chromosomes not in the labeled_chrs list
  chr_summary$label[!(chr_summary$Chromosome %in% labeled_chrs)] <- ""
  
  return(list(breaks = chr_summary$midpoint, labels = chr_summary$label))
}

#apply the custom scale function
chr_scale <- custom_chr_scale(myDiffUp)

ManhattanPlot <-
  ggplot() +
  geom_point(data = myDiffUp[!significant,], aes(x = cum_start, y = -log10(pvalue), col = CHR_State), shape = 16, size = 1, alpha = .4) +
  geom_point(data = myDiffUp[significant,], aes(x = cum_start, y = -log10(pvalue)), shape = 16, size = 2, alpha = .75, color = "blue") +
  scale_x_continuous("Chromosome",
                     breaks = chr_scale$breaks,
                     labels = chr_scale$labels,
                     expand = c(0, 0)) +
  scale_y_continuous(expression("-log"^"10"*"(p-value)"),
                     limits=c(0,50),
                     breaks = c(10, 20, 30, 40), 
                     labels = c("10", "20", "30", "40"),
                     expand = c(0, 0)) +
  scale_color_manual(values = c("azure4", "darkgrey")) +
  geom_hline(aes(yintercept = significance), col="firebrick") +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1, "cm"),
        axis.line = element_line(colour = "black", size = .3),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 8,color="grey34"),
        axis.text.y = element_text(size = 14), axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) + 
  guides(col="none")


#save the plot
png("Fig2_ZF_Manhattan_treatment_FINAL.png",  res = 300, width = 4000, height = 1200)
ManhattanPlot
dev.off()

# DMS annotation and stats ------------------------------------------------
#annotation and stats

#load file containing only DMS
#load file containing annotation of all CpG sites 
zf_anot<- readRDS("ZF_Annotation_All_cpGs.rds")


#merge the two by position
anot_all<- merge(significant, zf_anot, by="Pos")

prop_zf_anot<-  as.data.frame(prop.table(table(anot_all$category)))

sum(prop_zf_anot$Freq)

prop_zf_anot$Var1<- factor(prop_zf_anot$Var1, levels=c("Promoter", "Intron", "Exon", "Intergenic"))

prop_zf_anot$ObsExp<- "Observed"
colnames(prop_zf_anot)<- c("Location", "Proportion", "ObsExp")


zf_p <- read_xlsx("DNAm_BroodSize_DMS_Database_FINAL.xlsx", 
                          sheet = "Significant_p-value")

zf_p_anot<- merge(zf_p, zf_anot, by="Pos")

prop_zf_p_anot<-  as.data.frame(prop.table(table(zf_p_anot$category)))

sum(prop_zf_p_anot$Freq)

prop_zf_p_anot$Var1<- factor(prop_zf_p_anot$Var1, levels=c("Promoter", "Intron", "Exon", "Intergenic"))

colnames(prop_zf_p_anot)<- c("Location", "Proportion")

prop_zf_p_anot$ObsExp<- "Significant"

zf_loc_all_p<- rbind(prop_zf_anot, prop_zf_p_anot)

#expected all CpGs AFTER filtering
full_df_anot<- merge(myDiff, zf_anot, by="Pos")

prop_full_df_anot<-  as.data.frame(prop.table(table(full_df_anot$category)))

colnames(prop_full_df_anot)<- c("Location", "Proportion")

prop_full_df_anot$ObsExp<- "Expected"

zf_loc_all_p_full<- rbind(zf_loc_all_p,prop_full_df_anot)

# Figure 3 ----------------------------------------------------------------


zf_loc_all_plot<- ggplot(zf_loc_all_p_full, aes(x=factor(ObsExp, levels = c("Observed", "Significant", "Expected")), y=Proportion,fill=Location))+
  geom_bar( position = "fill", stat="identity", width=0.73)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_brewer(palette="RdBu")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        axis.text.x = element_text(colour = "black", angle = 0, size = 11),
        aspect.ratio = 1.2)+
  xlab(NULL)+  ylab(NULL)

png("Fig3_ZF_Annotation_FINAL.png",  res = 300, width = 2000, height = 1200)
zf_loc_all_plot
dev.off()

# stats_observed_expected -------------------------------------------------

#define observed proportions
observed <- data.frame(
  Location = c("Exon", "Intergenic", "Intron", "Promoter"),
  Proportion = c(0.147651007, 0.288590604, 0.557046980, 0.006711409),
  Counts = c(22, 43, 83, 1)
)

#define expected proportions
expected <- data.frame(
  Location = c("Exon", "Intergenic", "Intron", "Promoter"),
  Proportion = c(0.12737149, 0.37011756, 0.47497528 ,0.02753567)
)

#calculate expected counts based on the total observed counts
total_observed_counts <- sum(observed$Counts)
expected$Counts <- expected$Proportion * total_observed_counts

#combine into a single data frame
result <- data.frame(
  Location = observed$Location,
  Observed_Proportion = observed$Proportion,
  Observed_Counts = observed$Counts,
  Expected_Proportion = expected$Proportion,
  Expected_Counts = expected$Counts
)

#find difference in proportion between observed-expected
result$dif<- result$Observed_Proportion- result$Expected_Proportion

N <- sum(result$Observed_Counts)

#run exact binomial test per annotation category
with(result, map2(Observed_Counts, Expected_Proportion, function(Observed_Counts, Expected_Proportion) {
  
  binom.test(Observed_Counts, N, p = Expected_Proportion, alternative = "two.sided")
  
}))



# DMS identification only YOUNG-------------------------------------------------------------------


file.list <- list("S1.deduplicated.bismark.cov.gz",
                  "S2.deduplicated.bismark.cov.gz",
                  "S3.deduplicated.bismark.cov.gz",
                  "S4.deduplicated.bismark.cov.gz",
                  "S5.deduplicated.bismark.cov.gz",
                  "S7.deduplicated.bismark.cov.gz",
                  "S8.deduplicated.bismark.cov.gz",
                  "S9.deduplicated.bismark.cov.gz", 
                  "S10.deduplicated.bismark.cov.gz",
                  "S11.deduplicated.bismark.cov.gz",
                  "S12.deduplicated.bismark.cov.gz", 
                  "S13.deduplicated.bismark.cov.gz",
                  "S14.deduplicated.bismark.cov.gz", 
                  "S15.deduplicated.bismark.cov.gz",
                  "S16.deduplicated.bismark.cov.gz",
                  "S17.deduplicated.bismark.cov.gz", 
                  "S18.deduplicated.bismark.cov.gz", 
                  "S19.deduplicated.bismark.cov.gz", 
                  "S20.deduplicated.bismark.cov.gz",
                  "S21.deduplicated.bismark.cov.gz",
                  "S22.deduplicated.bismark.cov.gz", 
                  "S23.deduplicated.bismark.cov.gz", 
                  "S24.deduplicated.bismark.cov.gz",
                  "S25.deduplicated.bismark.cov.gz",
                  "S27.deduplicated.bismark.cov.gz", 
                  "S28.deduplicated.bismark.cov.gz",
                  "S29.deduplicated.bismark.cov.gz",
                  "S30.deduplicated.bismark.cov.gz",
                  "S32.deduplicated.bismark.cov.gz", 
                  "S34.deduplicated.bismark.cov.gz", 
                  "S35.deduplicated.bismark.cov.gz", 
                  "S36.deduplicated.bismark.cov.gz", 
                  "S37.deduplicated.bismark.cov.gz",
                  "S39.deduplicated.bismark.cov.gz", 
                  "S40.deduplicated.bismark.cov.gz",
                  "S42.deduplicated.bismark.cov.gz", 
                  "S43.deduplicated.bismark.cov.gz", 
                  "S44.deduplicated.bismark.cov.gz",
                  "S47.deduplicated.bismark.cov.gz",
                  "S49.deduplicated.bismark.cov.gz",
                  "S57.deduplicated.bismark.cov.gz",
                  "S59.deduplicated.bismark.cov.gz",
                  "S63.deduplicated.bismark.cov.gz",
                  "S65.deduplicated.bismark.cov.gz",
                  "S68.deduplicated.bismark.cov.gz",
                  "S74.deduplicated.bismark.cov.gz",
                  "S75.deduplicated.bismark.cov.gz",
                  "S80.deduplicated.bismark.cov.gz",
                  "S81.deduplicated.bismark.cov.gz",
                  "S85.deduplicated.bismark.cov.gz")

myobj <- methRead(file.list, sample.id=list("S1","S2", "S3", "S4", "S5", "S7", "S8", "S9", "S10", "S11", "S12", "S13", 
                                            "S14", "S15", "S16", "S17", "S18", "S19", "S20","S21", "S22", "S23", "S24", 
                                            "S25", "S27", "S28", "S29", "S30", "S32", "S34", "S35", "S36", "S37", "S39", 
                                            "S40", "S42", "S43", "S44", "S47", "S49", "S57", "S59", "S63", "S65", "S68", 
                                            "S74", "S75", "S80", "S81", "S85"), 
                  pipeline = "bismarkCoverage",
                  assembly="m2021",
                  treatment=c(1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1))

# Define coverage filtering thresholds
min_coverage <- 5
max_coverage <- 99.9

# Filter by coverage
methyl_filtered <- filterByCoverage(my_obj,
                                    lo.count = min_coverage,
                                    lo.perc = NULL,
                                    hi.count = max_coverage ,
                                    hi.perc = NULL)

# Normalize coverage
methyl_normalized <- normalizeCoverage(methyl_filtered)

# Unite the samples
meth_united <- unite(methyl_normalized)

meth_percentages <- percMethylation(meth_united)

# Identify positions to keep (not 0% or 100% in all samples)
positions_to_keep <- apply(meth_percentages, 1, function(x) !all(x %in% c(0, 100)))

# Filter the methylBase object
meth_united_filtered <- meth_united[positions_to_keep,]

# Now proceed with the differential methylation analysis using meth_united_filtered
myDiff_young <- calculateDiffMeth(meth_united_filtered, covariates=covariates, overdispersion=c("none"), adjust=c("bonferroni"))

# DMS identification only OLD-------------------------------------------------------------------

file.list <- list("S6.deduplicated.bismark.cov.gz",
                  "S26.deduplicated.bismark.cov.gz",
                  "S31.deduplicated.bismark.cov.gz",
                  "S33.deduplicated.bismark.cov.gz",
                  "S38.deduplicated.bismark.cov.gz",
                  "S41.deduplicated.bismark.cov.gz",
                  "S45.deduplicated.bismark.cov.gz",
                  "S46.deduplicated.bismark.cov.gz",
                  "S48.deduplicated.bismark.cov.gz",
                  "S50.deduplicated.bismark.cov.gz",
                  "S51.deduplicated.bismark.cov.gz",
                  "S52.deduplicated.bismark.cov.gz",
                  "S53.deduplicated.bismark.cov.gz",
                  "S54.deduplicated.bismark.cov.gz",
                  "S55.deduplicated.bismark.cov.gz",
                  "S56.deduplicated.bismark.cov.gz",
                  "S58.deduplicated.bismark.cov.gz",
                  "S60.deduplicated.bismark.cov.gz",
                  "S61.deduplicated.bismark.cov.gz",
                  "S62.deduplicated.bismark.cov.gz",
                  "S64.deduplicated.bismark.cov.gz",
                  "S66.deduplicated.bismark.cov.gz",
                  "S67.deduplicated.bismark.cov.gz",
                  "S69.deduplicated.bismark.cov.gz",
                  "S70.deduplicated.bismark.cov.gz",
                  "S71.deduplicated.bismark.cov.gz",
                  "S72.deduplicated.bismark.cov.gz",
                  "S73.deduplicated.bismark.cov.gz",
                  "S76.deduplicated.bismark.cov.gz",
                  "S77.deduplicated.bismark.cov.gz",
                  "S78.deduplicated.bismark.cov.gz",
                  "S79.deduplicated.bismark.cov.gz",
                  "S82.deduplicated.bismark.cov.gz",
                  "S83.deduplicated.bismark.cov.gz",
                  "S84.deduplicated.bismark.cov.gz",
                  "S86.deduplicated.bismark.cov.gz",
                  "S87.deduplicated.bismark.cov.gz",
                  "S88.deduplicated.bismark.cov.gz",
                  "S89.deduplicated.bismark.cov.gz",
                  "S90.deduplicated.bismark.cov.gz",
                  "S91.deduplicated.bismark.cov.gz",
                  "S92.deduplicated.bismark.cov.gz",
                  "S93.deduplicated.bismark.cov.gz",
                  "S94.deduplicated.bismark.cov.gz",
                  "S95.deduplicated.bismark.cov.gz",
                  "S96.deduplicated.bismark.cov.gz",
                  "S97.deduplicated.bismark.cov.gz",
                  "S98.deduplicated.bismark.cov.gz",
                  "S99.deduplicated.bismark.cov.gz",
                  "S100.deduplicated.bismark.cov.gz")


treatment<- c(1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 1)

myobj <- methRead(file.list, 
                  sample.id = list("S6", "S26", "S31", "S33", "S38", "S41", "S45", "S46", "S48", "S50", 
                                   "S51", "S52", "S53", "S54", "S55", "S56", "S58", "S60", "S61", "S62", 
                                   "S64", "S66", "S67", "S69", "S70", "S71", "S72", "S73", "S76", "S77", 
                                   "S78", "S79", "S82", "S83", "S84", "S86", "S87", "S88", "S89", "S90", 
                                   "S91", "S92", "S93", "S94", "S95", "S96", "S97", "S98", "S99", "S100"), 
                  pipeline = "bismarkCoverage",
                  assembly = "m2021",
                  treatment = treatment)

# Define coverage filtering thresholds
min_coverage <- 10
max_coverage <- 99.9

# Filter by coverage
methyl_filtered <- filterByCoverage(myobj,
                                    lo.count = min_coverage,
                                    lo.perc = NULL,
                                    hi.count = max_coverage ,
                                    hi.perc = NULL)

# Normalize coverage
methyl_normalized <- normalizeCoverage(methyl_filtered)

# Unite the samples
meth_united <- unite(methyl_normalized)

meth_percentages <- percMethylation(meth_united)

# Identify positions to keep (not 0% or 100% in all samples)
positions_to_keep <- apply(meth_percentages, 1, function(x) !all(x %in% c(0, 100)))

# Filter the methylBase object
meth_united_filtered <- meth_united[positions_to_keep,]

# Now proceed with the differential methylation analysis using meth_united_filtered
myDiff_old <- calculateDiffMeth(meth_united_filtered, covariates=covariates, overdispersion=c("none"), adjust=c("bonferroni"))

#comparison young and old
summary(comparedf(myDiff_young,myDiff_old))

myDiff_young$Age<- "Young"
myDiff_old$Age<- "Old"

full_df<- rbind(myDiff_young, myDiff_old)
nrow(full_df)
#94053

full_df$Pos<- paste(full_df$chr, full_df$start, sep="_")

shared_pos <- full_df %>%
  group_by(Pos) %>%
  summarize(n_ages = n_distinct(Age)) %>%
  filter(n_ages > 1) %>%
  pull(Pos)

# Then, create new variables for shared positions
full_df_shared <- full_df %>%
  filter(Pos %in% shared_pos) %>%
  pivot_wider(
    id_cols = Pos,
    names_from = Age,
    values_from = meth.diff
  )

#25764 shared sites

# Figure S1-------------------------------------------------
young_old_methdiff_plot<- ggplot(full_df_shared) +
  geom_point(aes(x = Young,y = Old), alpha=0.1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_smooth(aes(x = Young,y = Old), method="lm", size=1,col="magenta")+
  xlim(-26,26)+
  ylim(-26,26)+
  xlab("% difference in DNA methylation young")+
  ylab("% difference in methylation old")+
  theme(aspect.ratio=1, panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        axis.text.x = element_text(colour = "black", angle = 0, size = 11))+
  annotate("text", label = ~ "r=0.25" , x = -14, y = 20, size = 4)+
  annotate("text", label = "p<0.001" , x = -14, y = 16, size = 4)
ggsave("C:/Users/MWP-/Desktop/PhD/ZF_DNAm_2024/Methylkit/DMS_BroodSize_Final/YoungvsOld/oldvsYoung_MethDiff.png", plot = last_plot(), width=5, height=5)

png("FigS1_ZF_Young_Old_MethDiff_FINAL.png",  res = 300, width = 2000, height = 1200)
young_old_methdiff_plot
dev.off()


cor.test(full_df_shared$Old, full_df_shared$Young)


# Figure S2 ---------------------------------------------------------------
#FIND PROPORTION OF POST-FILTERING CPGS BEING DMS

#number of  CpGs per chromosome post-filtering
prop_chrom_anot <- myDiff %>%
  group_by(chr) %>%
  summarise(count = n(), .groups = 'drop') 

#chromosome names
chr_names <-read_excel("DNAm_BroodSize_DMS_Database_FINAL.xlsx", 
                                 sheet = "Chr_names")

prop_chrom_anot<- merge(prop_chrom_anot, chr_names, by="chr")

#number of DMS per chromosome 
all_chrom<- tabyl(significant$chr)

colnames(all_chrom)<- c("chr", "n", "percent")

all_chrom<- merge(all_chrom, chr_names, by="chr")

all_chrom<- all_chrom[,2:4]

colnames(all_chrom)<-c("n_obs", "percent","Chromosome") 

all_chrom_obs_exp<- merge(all_chrom, prop_chrom_anot, by="Chromosome")

all_chrom_obs_exp$prop<- all_chrom_obs_exp$n_obs/all_chrom_obs_exp$count

all_chrom_obs_exp$Chromosome <- factor(all_chrom_obs_exp$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '9', '10', '11', '12', '13', '14', '15','17', '18', '19',  '21', '26',  '28', '29'))


prom_dms_chr_plot<- ggplot(all_chrom_obs_exp, aes(x = Chromosome, y = prop)) +
  geom_segment(aes(x = Chromosome, xend = Chromosome, y = 0, yend = prop), color = "grey") +
  geom_point(color = "royalblue", size = 4) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.0012)) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    axis.text.x = element_text(colour = "black", size = 11),
    aspect.ratio = 0.5
  ) +
  ylab("Proportion of DMS") +
  xlab("Chromosome")  

png("FigS2_ZF_Prop_DMS_Chr_FINAL.png",  res = 300, width = 2000, height = 1200)
prom_dms_chr_plot
dev.off()


# Figure S3 ---------------------------------------------------------------

#(log) coverage per chromosome 

#data frame with coverage information per CpG per sample
d<-readRDS("CovZF_Apr2024.rds")

d$Chromosome <- factor(d$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30','31', '32','33', '34', '35', '36','37', 'W', 'Z'))
d$Chromosome <- factor(d$Chromosome, levels = rev(levels(d$Chromosome)))

#plot (log)coverage per chromosome
ggplot(d, aes(y = Chromosome, x = log(Coverage), fill = Chromosome))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1)+
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.1)+
  labs(y = "Chromosome",
       x = "log Coverage")+
  theme(aspect.ratio = 1.6, panel.border = element_rect(colour = "black", fill=NA, size=0.75), legend.position = "none")+
  scale_fill_viridis_d()+
  scale_y_discrete(expand=c(0,0))

png("FigS3_ZF_log_Coverage_Chr.png",  res = 300, width = 2000, height = 1200)
prom_dms_chr_plot
dev.off()

# Figure S4 ---------------------------------------------------------------
cov_n <- read_excel("DNAm_BroodSize_DMS_Database_FINAL.xlsx", 
                    sheet = "Coverage")

cov_n$Chromosome <- factor(cov_n$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32','33', '34', '35', '37','Z'))



median_cov<- ggplot(cov_n) +
  geom_point(aes(x = num_sites, y = median_coverage, 
                 color = ifelse(Chromosome %in% c("3", "29"), "highlight", "other")), 
             size = 3, alpha = 0.8) +
  geom_text(aes(x = num_sites, y = avg_coverage, 
                label = ifelse(Chromosome %in% c("3", "29"), Chromosome, "")), 
            hjust = -0.5, vjust = 0.6, size = 4, color = "black") +
  scale_color_manual(values = c("highlight" = "red", "other" = "gray")) +
  labs(x = "Number of sites after filtering", y = "Median coverage")  +
  theme(
    aspect.ratio = 0.6,
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  ) +
  ylim(1000, 2200) +
  guides(color = "none")

png("FigS4_ZF_Median_Coverage_Chr_FINAL.png",  res = 300, width = 2000, height = 1200)
median_cov
dev.off()




# Figure S5 ---------------------------------------------------------------

#file with DNAm information for all samples for all DMS and sample info (brood size)
m_info<- read_xlsx("DNAm_BroodSize_DMS_Database_FINAL.xlsx", 
                  sheet = "DNAm_BroodSize_DMS_SampleInfo")

m_info$StandRearNestSize<- as.factor(m_info$StandRearNestSize)

m_average_brood <- m_info %>%
  group_by(StandRearNestSize) %>%
  summarise(mean=mean(MethylationPercentage),
            se= sd(MethylationPercentage)/sqrt(n()))
#meth.diff if calculated as methylation in large broods-methylation in small broods


m_all<- ggplot(m_info,aes(x=StandRearNestSize, y=MethylationPercentage, col=StandRearNestSize))+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1)+
  geom_boxplot(notch = F,  outlier.size = -1,lwd=1, alpha = 0.7,show.legend = F)+ 
  guides(col="none")+
  theme(aspect.ratio = 1.4, panel.border = element_rect(colour = "black", fill=NA, size=0.75))+
  scale_color_manual(values=c("darkolivegreen", "darkorange1"))+
  stat_summary(fun = mean, geom = "errorbar", 
               aes(ymax = after_stat(y), ymin = after_stat(y), color = StandRearNestSize), 
               linetype = "dashed", size = 0.8)+
  scale_x_discrete(labels = c("2" = "small", "6" = "large")) +
  labs(x = "Brood size",
       y = "% DNA methylation of DMS")

png("FigS5_ZF_Methylation_Brood_Size_FINAL.png",  res = 300, width = 2000, height = 1200)
m_all
dev.off()


m_info_positive <- m_info[m_info$meth.diff > 0, ]
m_info_negative <- m_info[m_info$meth.diff < 0, ]

length(unique(m_info_positive$Pos))
length(unique(m_info_negative$Pos))

#calculate mean-median per brood size for CpGs increasing and decreasing in DNAm with brood size

m_info_positive_mean <- m_info_positive %>%
  group_by(StandRearNestSize) %>%
  summarise(mean=mean(MethylationPercentage),
            se= sd(MethylationPercentage)/sqrt(n()))

m_info_negative_mean <- m_info_negative %>%
  group_by(StandRearNestSize) %>%
  summarise(mean=mean(MethylationPercentage),
            se= sd(MethylationPercentage)/sqrt(n()))

m_info_positive_median <- m_info_positive %>%
  group_by(StandRearNestSize) %>%
  summarise(median=mean(MethylationPercentage),
            se= sd(MethylationPercentage)/sqrt(n()))

m_info_negative_median <- m_info_negative %>%
  group_by(StandRearNestSize) %>%
  summarise(mean=median(MethylationPercentage),
            se= sd(MethylationPercentage)/sqrt(n()))


