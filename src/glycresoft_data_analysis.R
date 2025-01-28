library(tidyverse)
library(lawstat)

# load pre-written functions
setwd("..")
source(file="src/functions/two_group_test.R")
source(file="src/functions/glycan_monomer_extractor.R")
source(file="src/functions/impute_small.R")

# ..................................................................... ####
# read in GlycReSoft results files ####
file_list <- list.files(path = "data/", pattern = "-glycopeptides")

for (i in 1:length(file_list)) {
  print(i)
  temp_fl <- read.csv(paste0("data/", file_list[i]))
  temp_fl$Sample <- rep(str_match(file_list[i], "GLYCO-(.*)\\(")[2], nrow(temp_fl))
  if (i==1) {
    raw_data <- temp_fl
    } else {
    raw_data <- rbind(raw_data, temp_fl)
    }
  }

# read in list of glycopeptides Identified by FragPipe
fp_gps <- read.csv("data/glycopeptides_by_fragpipe.csv")

# generate / extract categorical variables ####
# add sample name and grouping columns
raw_data$group <- str_match(raw_data$Sample, "a")[,1]
raw_data$group[is.na(raw_data$group)] <- "t"

# peptide, glycan, and glycopeptide columns are also generated and added
tmp_peptide_name <- rep(0,nrow(raw_data))
for (i in 1:nrow(raw_data)) {
  print(i)
  tmp_peptide_name[i] <- strsplit(gsub("\\(N-Glycosylation\\)|\\(Oxidation\\)|\\(Carbamidomethyl\\)", "*", raw_data$glycopeptide[i]), "\\{")[[1]][1]
  }
raw_data$Peptide <- tmp_peptide_name

tmp_glycan_name <- rep(0,nrow(raw_data))
for (i in 1:nrow(raw_data)) {
  print(i)
  tmp_1 <- strsplit(gsub("\\(N-Glycosylation\\)|\\(Oxidation\\)|\\(Carbamidomethyl\\)", "*", raw_data$glycopeptide[i]), "\\{")[[1]][2]
  tmp_glycan_name[i] <- gsub("[[:punct:]]","",tmp_1)
  }
raw_data$Glycan <- tmp_glycan_name

raw_data$Glycopeptide <- paste(raw_data$Peptide, raw_data$Glycan)

# extract glycan composition based on glycan name
tmp_glycan_monomers <- matrix(rep(0,4*nrow(raw_data)), ncol=4)
for (i in 1:nrow(raw_data)) {
  print(i)
  tmp_glycan_monomers[i,] <- glycan_monomer_extractor(raw_data$glycopeptide[i])
  }
raw_data$Hex <- tmp_glycan_monomers[,1]
raw_data$HexNAc <- tmp_glycan_monomers[,2]
raw_data$Neu5Ac <- tmp_glycan_monomers[,3]
raw_data$Fuc <- tmp_glycan_monomers[,4]

# calculate glycan type based on monomers
raw_data$glycan_type <- "Misc"
raw_data$glycan_type[raw_data$HexNAc > 3 & raw_data$Hex > 2] <- "Complex"
raw_data$glycan_type[raw_data$HexNAc == 3 & raw_data$Hex > 3] <- "Hybrid"
raw_data$glycan_type[raw_data$HexNAc == 2 & raw_data$Hex > 3] <- "Oligomannose"

# filter data: throw out NAs
# and glycopeptides that FP did not validate
# FP validation: filter out GlycReSoft glycopeptide matches that FragPipe did not identify - in fp_gps YES means exact match, ??? means the backbone has been identified, but not the same glycoform (same glycan structure)
backbone_to_filter <- fp_gps$Peptide[fp_gps$ACCEPTED == "NO"]

raw_data %>%
  filter(total_signal > 1) %>%
  filter(!gsub("\\*", "", Peptide) %in% backbone_to_filter) %>%
  select(Glycopeptide, Peptide, Glycan, glycan_type, protein_name, Sample, group, Hex, HexNAc, Neu5Ac, Fuc, apex_time, peptide_start, total_signal) -> filtered_data

# ..................................................................... ####
# GPROTEIN GROUPS ####
# dealing with glycosite homology ####
# first step: determine glycosylation sites with sequence homology (duplicates)
# glycosite homology is (for this specific dataset) due to different entries in uniprot for IGG and IGHG - IGG and IGG heavy chain
glycopeptide_duplicates <- filtered_data[duplicated(paste(filtered_data$Glycopeptide,filtered_data$Sample)),]

glycosite_homologs <- unique(glycopeptide_duplicates$Peptide)

# second step: create protein groups for duplicate glycopeptides
protein_group_names <- as.data.frame(rep(NA,length(glycosite_homologs)))
names(protein_group_names) <- "protein_groups"
rownames(protein_group_names) <- glycosite_homologs
for (i in 1:length(glycosite_homologs)) {
  print(i)
  filtered_data%>%
    filter(Peptide == glycosite_homologs[i]) -> tmp_protein_groups
  protein_group_names[i,] <- paste(unique(tmp_protein_groups$protein_name),collapse=" / ")
}

# third step: remove duplicate glycopeptides
filtered_data_no_d <- filtered_data[!duplicated(paste(filtered_data$Glycopeptide,filtered_data$Sample)),]

# NORMALIZATION ####
# total area normalization of the data ####
# this is arguably the most frequently used method (and very simple) for glycopeptide data in the literature, briefly, each signal in a sample is normalized by the total sample signal
filtered_data_no_d%>%
  group_by(Sample)%>%
  summarise(sum(total_signal)) -> total_sample_intensity

for (i in 1:nrow(filtered_data_no_d)) {
  print(i)
  filtered_data_no_d$TA_normalized_signal[i] <- filtered_data_no_d$total_signal[i] / total_sample_intensity$`sum(total_signal)`[total_sample_intensity$Sample %in% filtered_data_no_d$Sample[i]]
}

# IMPUTATION ####
# the purpose of this step is ONLY to provide a basis for the statistical comparison of glycopeptides, which means that when a gp is found reliably in one group, but not at all in an other one, there is something to compare to (for a test or fc calculation), these values are NOT stored or used for anything else
filtered_data_no_d%>%
  filter(group == "t")%>%
  group_by(Glycopeptide)%>%
  summarise(n()) -> test1
filtered_data_no_d%>%
  filter(group == "a")%>%
  group_by(Glycopeptide)%>%
  summarise(n()) -> test2

unique_glycopeptides <- unique(filtered_data_no_d$Glycopeptide)

count_of_valid_values <- as.data.frame(matrix(rep(0, length(unique_glycopeptides)*2), ncol=2))
names(count_of_valid_values) <- c("t","a")
rownames(count_of_valid_values) <- unique_glycopeptides

for (i in 1:length(unique_glycopeptides)) {
  print(i)
  if (sum(test1$Glycopeptide %in% unique_glycopeptides[i]) == 0) {} else {
    count_of_valid_values[i,1] <- test1$`n()`[test1$Glycopeptide %in% unique_glycopeptides[i]]
  }
  if (sum(test2$Glycopeptide %in% unique_glycopeptides[i]) == 0) {} else {
    count_of_valid_values[i,2] <- test2$`n()`[test2$Glycopeptide %in% unique_glycopeptides[i]]
  }
}

# the rules of filtering and imputation are as follows:
# 1: a gp is considered if it is found in at least 5 samples in any group
# 2: a gp is imputed if i, rule 1 is true, and ii, it is found in less than 3 samples
# 3: imputation is group-based, and is drawn from a small normal distribution (sometime referred as to "Perseus" imputation, or QRILC)
# IMPORTANT: imputed values are not stored, they only exist temporarily in the statistical test for cycle
# these are the glycopeptides that need to be imputed
glycopeptide_to_impute <- count_of_valid_values[count_of_valid_values$t > 4 & count_of_valid_values$a < 3 | count_of_valid_values$t < 3 & count_of_valid_values$a > 4,]
gp_impute_tumor <- rownames(glycopeptide_to_impute[glycopeptide_to_impute$t < 3,])
gp_impute_adjacent <- rownames(glycopeptide_to_impute[glycopeptide_to_impute$a < 3,])

# number of samples only found in a single sample
sum(count_of_valid_values$t + count_of_valid_values$a == 1)

# ..................................................................... ####
# GLYCOPEPTIDES ####
# ... run statistical tests for glycopeptides ####
# the two_group_test is a function that tests for normality and variance equality and chooses the appropriate statistical test to run
# usage: two_group_test(data_1 (TUMOR), data_2 (NORMAL), log_data = T or F)
# log_data determines how FC is calculated if T then mean_1 / mean_2, if F then mean_1 - mean_2
# output: list with p-value, test_type, fold_change

length(unique(filtered_data_no_d$Glycopeptide))

# eligible glycopeptides are the ones with sufficient number of values / group
eligible_glycopeptides <- unique(c(test1$Glycopeptide[test1$`n()` > 4],test2$Glycopeptide[test2$`n()` > 4]))

NvsT_comparison <- as.data.frame(matrix(rep(NA,3*length(eligible_glycopeptides)),ncol=3))
names(NvsT_comparison) <- c("p-value","test_type","fold-change")
row.names(NvsT_comparison) <- eligible_glycopeptides

for (i in 1:length(eligible_glycopeptides)) {
  print(i)
  filtered_data_no_d%>%
    filter(group == "t")%>%
    filter(Glycopeptide == eligible_glycopeptides[i])%>%
    select(TA_normalized_signal) -> gp_in_tumor
  filtered_data_no_d%>%
    filter(group == "a")%>%
    filter(Glycopeptide == eligible_glycopeptides[i])%>%
    select(TA_normalized_signal) -> gp_in_adjacent
  NvsT_comparison$mean.t[i] <- mean(gp_in_tumor$TA_normalized_signal, na.rm = T)
  NvsT_comparison$sd.t[i] <- sd(gp_in_tumor$TA_normalized_signal, na.rm = T)
  NvsT_comparison$mean.n[i] <- mean(gp_in_adjacent$TA_normalized_signal, na.rm = T)
  NvsT_comparison$sd.n[i] <- sd(gp_in_adjacent$TA_normalized_signal, na.rm = T)
  if (nrow(gp_in_adjacent) < 3) {
    tmp_2test <- two_group_test(gp_in_tumor$TA_normalized_signal,c(gp_in_adjacent$TA_normalized_signal, impute_small(filtered_data_no_d$TA_normalized_signal,5)),log_data=F)
    NvsT_comparison[i,1] <- tmp_2test[[1]]
    NvsT_comparison[i,2] <- tmp_2test[[2]]
    NvsT_comparison[i,3] <- tmp_2test[[3]]
    } else if (nrow(gp_in_tumor) < 3) {
    tmp_2test <- two_group_test(c(gp_in_tumor$TA_normalized_signal,impute_small(filtered_data_no_d$TA_normalized_signal,5)),gp_in_adjacent$TA_normalized_signal,log_data=F)
    NvsT_comparison[i,1] <- tmp_2test[[1]]
    NvsT_comparison[i,2] <- tmp_2test[[2]]
    NvsT_comparison[i,3] <- tmp_2test[[3]]
    } else {
    tmp_2test <- two_group_test(gp_in_tumor$TA_normalized_signal,gp_in_adjacent$TA_normalized_signal,log_data=F)
    NvsT_comparison[i,1] <- tmp_2test[[1]]
    NvsT_comparison[i,2] <- tmp_2test[[2]]
    NvsT_comparison[i,3] <- tmp_2test[[3]]
    }
}

sum(rownames(NvsT_comparison) %in% rownames(glycopeptide_to_impute))
sum(count_of_valid_values$t > 3 & count_of_valid_values$a > 3)

sum(p.adjust(NvsT_comparison$`p-value`, method="BH") < 0.05)

filtered_data_no_d$Glycopeptide %in% rownames(NvsT_comparison)
NvsT_comparison$protein <- NA
for (i in 1:nrow(NvsT_comparison)) {
  print(i)
  NvsT_comparison$protein[i] <- filtered_data_no_d$protein_name[filtered_data_no_d$Glycopeptide %in% rownames(NvsT_comparison)[i]]
}

# ... proteomics normalization ####
# here, I am trying to normalize glycopeptide expression based on proteomics
# i need to find a way to match them
proteomics_fc <- read.csv("proteomics_NvsT.csv")

filtered_data_eligible <- filtered_data_no_d[filtered_data_no_d$Glycopeptide %in% eligible_glycopeptides,]


proteomics_pnames <- gsub("^[^|]*\\|[^|]*\\|([^;]*).*", "\\1", proteomics_fc$X)
glycoprot_pnames <- gsub("^[^|]*\\|[^|]*\\|([^ ]*).*", "\\1", filtered_data_eligible$protein_name)

NvsT_comparison$adj_pvalue <- p.adjust(NvsT_comparison$`p-value`, method="BH")

for (i in 1:length(NvsT_comparison$adj_pvalue)) {
  NvsT_comparison$Protein[i] <- gsub("^[^|]*\\|[^|]*\\|([^ ]*).*", "\\1", unique(filtered_data_eligible$protein_name[filtered_data_eligible$Glycopeptide %in% row.names(NvsT_comparison)[i]]))
  if (sum(proteomics_pnames %in% NvsT_comparison$Protein[i]) == 0) {
    NvsT_comparison$Prot.pval[i] <- NA
    NvsT_comparison$Prot.fc[i] <- NA
  } else {
  NvsT_comparison$Prot.pval[i] <- proteomics_fc$adj.P.Val[proteomics_pnames %in% NvsT_comparison$Protein[i]]
  NvsT_comparison$Prot.fc[i] <- proteomics_fc$logFC[proteomics_pnames %in% NvsT_comparison$Protein[i]]
  }}

# write.csv(NvsT_comparison, file="results/glycopeptide_comparison_data.csv")

# ............................................................. ####
# METRICS ####
# calculate metrics ####
# ... calculate sia/fuc/galylation for complex and hybrid glycans ####
# calculate sialylation
unique_peptides <- unique(filtered_data_no_d$Peptide)
samples <- unique(filtered_data_no_d$Sample)
sialylation_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)), ncol=length(samples)))
names(sialylation_df) <- samples
row.names(sialylation_df) <- unique_peptides

for (i in 1:length(unique_peptides)) {
  print(i)
  for (j in 1:length(samples)) {
    print(j)
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(Sample == samples[j]) -> tmp_s1
    tmp_s2 <- sum(tmp_s1$TA_normalized_signal*tmp_s1$Neu5Ac,na.rm=T)
    tmp_s3 <- sum(tmp_s1$TA_normalized_signal*(tmp_s1$HexNAc-2-(1-(tmp_s1$Hex-tmp_s1$HexNAc))),na.rm=T)
    sialylation_df[i,j] <- tmp_s2/tmp_s3
  }
}

# calculate fucosylation
fucosylation_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)), ncol=length(samples)))
names(fucosylation_df) <- samples
row.names(fucosylation_df) <- unique_peptides

for (i in 1:length(unique_peptides)) {
  print(i)
  for (j in 1:length(samples)) {
    print(j)
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(Sample == samples[j]) -> tmp_f1
    fucosylation_df[i,j] <- sum(tmp_f1$TA_normalized_signal*tmp_f1$Fuc,na.rm=T)/sum(tmp_f1$TA_normalized_signal,na.rm=T)
  }
}

# calculate galactosylation
galactosylation_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)), ncol=length(samples)))
names(galactosylation_df) <- samples
row.names(galactosylation_df) <- unique_peptides

for (i in 1:length(unique_peptides)) {
  print(i)
  for (j in 1:length(samples)) {
    print(j)
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(Sample == samples[j]) -> tmp_g1
    galactosylation_df[i,j] <- sum(tmp_g1$TA_normalized_signal*(tmp_g1$Hex-3),na.rm=T)/sum(tmp_g1$TA_normalized_signal*(tmp_g1$HexNAc-2),na.rm=T)
  }
}

# ... calculate glycan type distribution ####
glycan_type_complex_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)),ncol=length(samples)))
names(glycan_type_complex_df) <- samples
rownames(glycan_type_complex_df) <- unique_peptides
glycan_type_hm_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)),ncol=length(samples)))
names(glycan_type_hm_df) <- samples
rownames(glycan_type_hm_df) <- unique_peptides
glycan_type_h_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)),ncol=length(samples)))
names(glycan_type_h_df) <- samples
rownames(glycan_type_h_df) <- unique_peptides

for (i in 1:length(unique_peptides)) {
  print(i)
  for (j in 1:length(samples)) {
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(glycan_type == "Complex")%>%
      filter(Sample == samples[j])  -> temp_gtc
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(glycan_type == "Oligomannose")%>%
      filter(Sample == samples[j]) -> temp_gthm
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(glycan_type == "Hybrid")%>%
      filter(Sample == samples[j]) -> temp_gth
    if(nrow(temp_gtc)>0) {
      glycan_type_complex_df[i,j] <- sum(temp_gtc$TA_normalized_signal)/sum(temp_gtc$TA_normalized_signal,temp_gthm$TA_normalized_signal,temp_gth$TA_normalized_signal)
      } else {}
    if(nrow(temp_gthm)>0) {
      glycan_type_hm_df[i,j] <- sum(temp_gthm$TA_normalized_signal)/sum(temp_gtc$TA_normalized_signal,temp_gthm$TA_normalized_signal,temp_gth$TA_normalized_signal)
      } else {}
    if(nrow(temp_gth)>0) {
      glycan_type_h_df[i,j] <- sum(temp_gth$TA_normalized_signal)/sum(temp_gtc$TA_normalized_signal,temp_gthm$TA_normalized_signal,temp_gth$TA_normalized_signal)
    } else {}
  }
}

# run statistical tests for calculated metrics ####
# not clear how this will work, with many identical values (0 or 1)
# okay so as it turns out identical values really mess up the statistical tests, this is a problem
# this should probably work better for only eligible glycopeptides
NvsT_comparison_metrics <- as.data.frame(matrix(rep(NA,6*length(unique_peptides)),ncol=6))
names(NvsT_comparison_metrics) <- c("Sia p-value","Fuc p-value","Gal p-value","Complex p-value","Oligmonannose p-value","Híbrid p-value")
row.names(NvsT_comparison_metrics) <- unique_peptides

for (i in 1:length(unique_peptides)) {
  print(i)
  sia_a <- unlist(sialylation_df[unique_peptides[i],str_detect(names(sialylation_df),"a")])
  sia_t <- unlist(sialylation_df[unique_peptides[i],!str_detect(names(sialylation_df),"t")])
  if (sum(!is.na(sia_a)) < 3 | length(unique(sia_a)) < 3 | sum(!is.na(sia_t)) < 3 | length(unique(sia_t)) < 3) {
    print("too much missing")
  } else {
    NvsT_comparison_metrics[i,1] <- two_group_test(sia_a[!is.na(sia_a)], sia_t[!is.na(sia_t)],log_data=F)[[1]]
  }
}
for (i in 1:length(unique_peptides)) {
  print(i)
  fuc_a <- unlist(fucosylation_df[unique_peptides[i],str_detect(names(fucosylation_df),"a")])
  fuc_t <- unlist(fucosylation_df[unique_peptides[i],!str_detect(names(fucosylation_df),"t")])
  if (sum(!is.na(fuc_a)) < 3  | length(unique(fuc_a)) < 3 | sum(!is.na(fuc_t)) < 3 | length(unique(fuc_t)) < 3) {
    print("too much missing")
  } else {
    NvsT_comparison_metrics[i,2] <- two_group_test(fuc_a[!is.na(fuc_a)], fuc_t[!is.na(fuc_t)],log_data=F)[[1]]
  }
}
# I have no idea why this happens - this is a quick fix
galactosylation_df[galactosylation_df == "Inf"] <- NA
galactosylation_df[galactosylation_df == "NaN"] <- NA
for (i in 1:length(unique_peptides)) {
  print(i)
  gal_a <- unlist(galactosylation_df[unique_peptides[i],str_detect(names(galactosylation_df),"a")])
  gal_t <- unlist(galactosylation_df[unique_peptides[i],!str_detect(names(galactosylation_df),"t")])
  if (sum(!is.na(gal_a)) < 3  | length(unique(gal_a)) < 3 | sum(!is.na(gal_t)) < 3 | length(unique(gal_t)) < 3) {
    print("too much missing")
  } else {
    NvsT_comparison_metrics[i,3] <- two_group_test(gal_a[!is.na(gal_a)], gal_t[!is.na(gal_t)],log_data=F)[[1]]
  }
}
for (i in 1:length(unique_peptides)) {
  print(i)
  complex_a <- unlist(glycan_type_complex_df[unique_peptides[i],str_detect(names(glycan_type_complex_df),"a")])
  complex_t <- unlist(glycan_type_complex_df[unique_peptides[i],!str_detect(names(glycan_type_complex_df),"t")])
  if (sum(!is.na(complex_a)) < 3| length(unique(complex_a)) < 3 | sum(!is.na(complex_t)) < 3 |length(unique(complex_t)) < 3) {
    print("too much missing")
  } else {
    NvsT_comparison_metrics[i,4] <- two_group_test(complex_a[!is.na(complex_a)], complex_t[!is.na(complex_t)],log_data=F)[[1]]
  }
}
for (i in 1:length(unique_peptides)) {
  print(i)
  hm_a <- unlist(glycan_type_hm_df[unique_peptides[i],str_detect(names(glycan_type_hm_df),"a")])
  hm_t <- unlist(glycan_type_hm_df[unique_peptides[i],!str_detect(names(glycan_type_hm_df),"t")])
  if (sum(!is.na(hm_a)) < 3  | length(unique(hm_a)) < 3 | sum(!is.na(hm_t)) < 3 | length(unique(hm_t)) < 3) {
    print("too much missing")
  } else {
    NvsT_comparison_metrics[i,5] <- two_group_test(hm_a[!is.na(hm_a)], hm_t[!is.na(hm_t)],log_data=F)[[1]]
  }
}
for (i in 1:length(unique_peptides)) {
  print(i)
  h_a <- unlist(glycan_type_h_df[unique_peptides[i],str_detect(names(glycan_type_h_df),"a")])
  h_t <- unlist(glycan_type_h_df[unique_peptides[i],!str_detect(names(glycan_type_h_df),"t")])
  if (sum(!is.na(h_a)) < 3  | length(unique(h_a)) < 3 | sum(!is.na(h_t)) < 3 | length(unique(h_t)) < 3) {
    print("too much missing")
  } else {
    NvsT_comparison_metrics[i,6] <- two_group_test(h_a[!is.na(h_a)], h_t[!is.na(h_t)],log_data=F)[[1]]
  }
}

# calculate OVERALL metrics ####
# calculate OVERALL sialylation
overall_sialylation <- as.data.frame(matrix(rep(NA,length(samples)),ncol=length(samples)))
names(overall_sialylation) <- samples
for (i in 1:length(samples)) {
  print(i)
  filtered_data_no_d%>%
      filter(Sample == samples[i]) -> tmp_s11
    tmp_s22 <- sum(tmp_s11$TA_normalized_signal*tmp_s11$Neu5Ac,na.rm=T)
    tmp_s33 <- sum(tmp_s11$TA_normalized_signal*(tmp_s11$HexNAc-2-(1-(tmp_s11$Hex-tmp_s11$HexNAc))),na.rm=T)
    overall_sialylation[,i] <- tmp_s22/tmp_s33
}
# changes in sialylation are NOT significant - for TA normalized data
sia_overall_a <- unlist(overall_sialylation[,str_detect(names(overall_sialylation),"a")])
sia_overall_t <- unlist(overall_sialylation[,!str_detect(names(overall_sialylation),"t")])
two_group_test(sia_overall_a,sia_overall_t,log_data=F)

# calculate OVERALL fucosylation
overall_fucosylation <- as.data.frame(matrix(rep(NA,length(samples)),ncol=length(samples)))
names(overall_fucosylation) <- samples
for (i in 1:length(samples)) {
  print(i)
  filtered_data_no_d%>%
    filter(Sample == samples[i]) -> tmp_s11
  overall_fucosylation[,i] <- sum(tmp_s11$TA_normalized_signal*tmp_s11$Fuc,na.rm=T)/sum(tmp_s11$TA_normalized_signal,na.rm=T)
}
# changes in fucosylation are NOT significant - for TA normalized data
fuc_overall_a <- unlist(overall_fucosylation[,str_detect(names(overall_fucosylation),"a")])
fuc_overall_t <- unlist(overall_fucosylation[,!str_detect(names(overall_fucosylation),"t")])
two_group_test(fuc_overall_a,fuc_overall_t,log_data=F)

# calculate OVERALL galactosylation
overall_galactosylation <- as.data.frame(matrix(rep(NA,length(samples)),ncol=length(samples)))
names(overall_galactosylation) <- samples
for (i in 1:length(samples)) {
  print(i)
  filtered_data_no_d%>%
    filter(Sample == samples[i])%>%
    filter(glycan_type == "Complex") -> tmp_s11
  overall_galactosylation[,i] <- sum(tmp_s11$TA_normalized_signal*(tmp_s11$Hex-3),na.rm=T)/sum(tmp_s11$TA_normalized_signal*(tmp_s11$HexNAc-2),na.rm=T)
}
overall_galactosylation[overall_galactosylation > 1] <- 1
# changes in galactosylation are NOT significant - for TA normalized data
gal_overall_a <- unlist(overall_galactosylation[,str_detect(names(overall_galactosylation),"a")])
gal_overall_t <- unlist(overall_galactosylation[,!str_detect(names(overall_galactosylation),"t")])
two_group_test(gal_overall_a,gal_overall_t,log_data=F)

# calculate OVERALL metrics
glycan_type_complex_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)),ncol=length(samples)))
names(glycan_type_complex_df) <- samples
rownames(glycan_type_complex_df) <- unique_peptides
glycan_type_hm_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)),ncol=length(samples)))
names(glycan_type_hm_df) <- samples
rownames(glycan_type_hm_df) <- unique_peptides
glycan_type_h_df <- as.data.frame(matrix(rep(NA, length(unique_peptides)*length(samples)),ncol=length(samples)))
names(glycan_type_h_df) <- samples
rownames(glycan_type_h_df) <- unique_peptides

for (i in 1:length(unique_peptides)) {
  print(i)
  for (j in 1:length(samples)) {
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(glycan_type == "Complex")%>%
      filter(Sample == samples[j])  -> temp_gtc
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(glycan_type == "Oligomannose")%>%
      filter(Sample == samples[j]) -> temp_gthm
    filtered_data_no_d%>%
      filter(Peptide == unique_peptides[i])%>%
      filter(glycan_type == "Hybrid")%>%
      filter(Sample == samples[j]) -> temp_gth
    if(nrow(temp_gtc)>0) {
      glycan_type_complex_df[i,j] <- sum(temp_gtc$TA_normalized_signal)/sum(temp_gtc$TA_normalized_signal,temp_gthm$TA_normalized_signal,temp_gth$TA_normalized_signal)
    } else {}
    if(nrow(temp_gthm)>0) {
      glycan_type_hm_df[i,j] <- sum(temp_gthm$TA_normalized_signal)/sum(temp_gtc$TA_normalized_signal,temp_gthm$TA_normalized_signal,temp_gth$TA_normalized_signal)
    } else {}
    if(nrow(temp_gth)>0) {
      glycan_type_h_df[i,j] <- sum(temp_gth$TA_normalized_signal)/sum(temp_gtc$TA_normalized_signal,temp_gthm$TA_normalized_signal,temp_gth$TA_normalized_signal)
    } else {}
  }
}

# ............................................................. ####
# Create plots ####
# ... boxplot of overall fucosylation ####
fucosylation_box <- as.data.frame(cbind(t(overall_fucosylation),
!grepl("a", names(overall_fucosylation))))

names(fucosylation_box) <- c("Values", "Class")

fucosylation_box$Class <- as.factor(fucosylation_box$Class)
levels(fucosylation_box$Class) <- c("Normal","Tumor")

png(filename = paste0("results/fucosylation_box", ".png"), units = "px", width = 1000, height = 2000, type = "cairo", res = 300)
ggplot(fucosylation_box, aes(Class, Values, fill = Class)) +
  geom_boxplot()+
  labs(title = "Fucosylation", y = "TA normalized LFQ Intensity", x = "") +
  geom_jitter(width = 0, size = 1.5) +
  stat_summary(fun=mean, geom="point", shape=21, size=2.5, color="black", fill="white", stroke = 1) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.position = "none",
        plot.margin = unit(c(10,15,10,15), "points"),
        plot.title = element_text(colour = "black", size = 18, hjust = 0, family = "sans"),
        plot.tag = element_text(colour = "black", size = 20, face = "bold", family = "sans"),
        axis.title.x = element_text(color = "black", size = 16, family = "sans", vjust = -0.5),
        axis.text = element_text(color = "black", size = 16, family = "sans"),
        axis.title.y = element_text(color = "black", size = 16, family = "sans", vjust = 2.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# ... boxplot of overall sialylation ####
sialylation_box <- as.data.frame(cbind(t(overall_sialylation),
                                        !grepl("a", names(overall_sialylation))))

names(sialylation_box) <- c("Values", "Class")

sialylation_box$Class <- as.factor(sialylation_box$Class)
levels(sialylation_box$Class) <- c("Normal","Tumor")

png(filename = paste0("results/sialylation_box", ".png"), units = "px", width = 1000, height = 2000, type = "cairo", res = 300)
ggplot(sialylation_box, aes(Class, Values, fill = Class)) +
  geom_boxplot()+
  labs(title = "Sialylation", y = "TA normalized LFQ Intensity", x = "") +
  geom_jitter(width = 0, size = 1.5) +
  stat_summary(fun=mean, geom="point", shape=21, size=2.5, color="black", fill="white", stroke = 1) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.position = "none",
        plot.margin = unit(c(10,15,10,15), "points"),
        plot.title = element_text(colour = "black", size = 18, hjust = 0, family = "sans"),
        plot.tag = element_text(colour = "black", size = 20, face = "bold", family = "sans"),
        axis.title.x = element_text(color = "black", size = 16, family = "sans", vjust = -0.5),
        axis.text = element_text(color = "black", size = 16, family = "sans"),
        axis.title.y = element_text(color = "black", size = 16, family = "sans", vjust = 2.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# ... boxplot for overall galactosylation ####
galactosylation_box <- as.data.frame(cbind(t(overall_galactosylation),
                                        !grepl("a", names(overall_galactosylation))))

names(galactosylation_box) <- c("Values", "Class")

galactosylation_box$Class <- as.factor(galactosylation_box$Class)
levels(galactosylation_box$Class) <- c("Normal","Tumor")

png(filename = paste0("results/galactosylation_box", ".png"), units = "px", width = 1000, height = 2000, type = "cairo", res = 300)
ggplot(galactosylation_box, aes(Class, Values, fill = Class)) +
  geom_boxplot()+
  labs(title = "Galactosylation", y = "TA normalized LFQ Intensity", x = "") +
  geom_jitter(width = 0, size = 1.5) +
  stat_summary(fun=mean, geom="point", shape=21, size=2.5, color="black", fill="white", stroke = 1) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(legend.position = "none",
        plot.margin = unit(c(10,15,10,15), "points"),
        plot.title = element_text(colour = "black", size = 18, hjust = 0, family = "sans"),
        plot.tag = element_text(colour = "black", size = 20, face = "bold", family = "sans"),
        axis.title.x = element_text(color = "black", size = 16, family = "sans", vjust = -0.5),
        axis.text = element_text(color = "black", size = 16, family = "sans"),
        axis.title.y = element_text(color = "black", size = 16, family = "sans", vjust = 2.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# write_results
write.csv(NvsT_comparison, file = "results/NvsT_comparison_glycopeptide.csv")
