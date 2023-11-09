# Compare AMR: BLAST vs AMR++v2
# Sophia Kenney - 9 November 2023

#load packages
library(BiocManager)
library(tidyverse)
library(tidylog)
library(phyloseq)
library(microViz)
library(dada2)
library(decontam)

# ---- Prepare AMR Ouput for psobject ----
#read in tables 
app <- read.delim("tables/amr/v2meltdf.txt", sep = "\t")
blast <- read.delim("tables/amr/blastcounts.txt", sep = "\t")
met <- read.csv("tables/meta/turkey_metadata.csv")
gene <- read.csv("tables/meta/megares_full_annotations_v2.00.csv")

#transform to long format
app <- app %>%
  pivot_wider(names_from = Project_ID, values_from = Counts) %>%
  as.data.frame()

blast <- blast %>%
  select(sample, pattern, genecount) %>%
  pivot_wider(names_from = sample, values_from = genecount, values_fill = 0)%>%
  as.data.frame()

#collaps rows in blast by gene
blast <- cbind(blast[,2:ncol(blast)], as.data.frame(str_split_fixed(blast$pattern, "\\|", 4))[4]) #extract gene name from amr++ output and add to count matrix

rownames(app) <- app$Gene
app <- app[,2:ncol(app)]

rownames(blast) <- blast$V4
blast <- blast[,1:ncol(blast)-1]

#save as matrices
mat1 <- as.matrix(app)
mat2 <- as.matrix(blast)

#save as otutab
otu1 <- phyloseq::otu_table(mat1, taxa_are_rows = TRUE)
otu2 <- phyloseq::otu_table(mat2, taxa_are_rows = TRUE)

# ---- Prepare Meta and Tax for phyloseq ----
#filter for samples shotgun sequenced
met1 <- rbind(met %>%
               filter(Project_ID %in% colnames(app)), 
              #add controls sequenced when samples were resequnced
              data.frame(Project_ID = c("T780", "T781"),
                        Sample_Type = c("DNA_PC", "DNA_NC"),
                        Animal_ID = c("DNA_PC", "DNA_NC"),
                        Timepoint = c(NA, NA),
                        Treatment = c("DNA_PC", "DNA_NC")))

met2 <- met %>%
  filter(Project_ID %in% colnames(blast)) # controls for AMC were sequenced

#save as sample data
samp1 <- phyloseq::sample_data(met1)
rownames(samp1) <- samp1$Project_ID

samp2 <- phyloseq::sample_data(met2)
rownames(samp2) <- samp2$Project_ID


#Prepare taxa table
#filter from megdb for just those entries relevant to genes detected
taxDF1 <- gene %>%
  select(type, class, group) %>%
  unique() %>%
  filter(group %in% rownames(mat1))

taxDF2 <- gene %>%
  select(type, class, group) %>%
  unique() %>%
  filter(group %in% rownames(mat2))

tax1 <- as.matrix(taxDF1)
rownames(tax1) <- taxDF1$group
taxtab1 <- phyloseq::tax_table(tax1)
taxa_names(taxtab1) #check to make sure gene names

tax2 <- as.matrix(taxDF2)
rownames(tax2) <- taxDF2$group
taxtab2 <- phyloseq::tax_table(tax2)
taxa_names(taxtab2) #check to make sure gene names

# ---- Create psobj ----

ps_app <- phyloseq::phyloseq(otu1, taxtab1, samp1) # create ps for AMR++ out
ps_blN <- phyloseq::phyloseq(otu2, taxtab2, samp2) # create ps for BlastN out

#save these
saveRDS(ps_app, "rdata/ps/ps_app.rds")
saveRDS(ps_blN, "rdata/ps/ps_blN.rds")

# ---- Decontam and Filtering----

# how many controls at each step (extraction and other)
ps_app %>% samdat_tbl() %>% group_by(Treatment) %>% summarize(n()) #2 extraction controls + air/sock NCs
ps_blN %>% samdat_tbl() %>% group_by(Treatment) %>% summarize(n()) #air/sock NCs

# validate
psA <- tax_fix(ps_app)
psB <- tax_fix(ps_blN)

##filter for relative abundance##

# transform to relative abundance
psAr <- transform_sample_counts(psA, function(x) x / sum(x))
psBr <- transform_sample_counts(psB, function(x) x / sum(x))

# remove taxa with total relative abundance less than 10e-5
psAr <- filter_taxa(psAr, function(x) mean(x) > 1e-5, TRUE) # 1187 genes
psBr <- filter_taxa(psBr, function(x) mean(x) > 1e-5, TRUE) # 425

##decontam## 

# add sampling variables
psAr <- psAr %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"NC"), true = "Control", false = "Sample")
  ) %>%
  ps_mutate(is.neg = if_else(SampleBinary == "Control", true = TRUE, false = FALSE))

psBr <- psBr %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"NC"), true = "Control", false = "Sample")
  ) %>%
  ps_mutate(is.neg = if_else(SampleBinary == "Control", true = TRUE, false = FALSE))

#check
sample_data(psAr)
sample_data(psBr)

## use negative controls to identify contaminants - none identified
contamdf.prevA <- isContaminant(psAr, method="prevalence", neg="is.neg", threshold = 0.5) # more aggressive threshold is 0.5
table(contamdf.prevA$contaminant) # this identifies 6 as contaminants

contamdf.prevB <- isContaminant(psBr, method="prevalence", neg="is.neg", threshold = 0.5) # more aggressive threshold is 0.5
table(contamdf.prevB$contaminant) # this identifies 20 as contaminants


#get list of contaminants
contamsA <- rownames(contamdf.prevA[contamdf.prevA$contaminant == "TRUE",])
contamsA

contamsB <- rownames(contamdf.prevB[contamdf.prevB$contaminant == "TRUE",])
contamsB

# remove contaminant sequences 
#from relabun filtered ps
nocontamAr <- prune_taxa(!rownames(psAr@tax_table) %in% contamsA, psAr )
nocontamBr <- prune_taxa(!rownames(psBr@tax_table) %in% contamsB, psBr )

#for counts
nocontamA <- prune_taxa(rownames(psAr@tax_table), psA) #drop taxa filtered based on relabun
nocontamA <- prune_taxa(!rownames(psAr@tax_table) %in% contamsA, nocontamA ) #drop taxa deemed contaminants

nocontamB <- prune_taxa(rownames(psBr@tax_table), psB) #drop taxa filtered based on relabun
nocontamB <- prune_taxa(!rownames(psBr@tax_table) %in% contamsB, nocontamB ) #drop taxa deemed contaminants

## filter SNP genes
nocontamA <- nocontamA %>% tax_select(tax_list = "SNP", strict_matches = FALSE, deselect = TRUE)
nocontamAr <- nocontamAr %>% tax_select(tax_list = "SNP", strict_matches = FALSE, deselect = TRUE)

nocontamB <- nocontamB %>% tax_select(tax_list = "SNP", strict_matches = FALSE, deselect = TRUE)
nocontamBr <- nocontamBr %>% tax_select(tax_list = "SNP", strict_matches = FALSE, deselect = TRUE)

#save these
saveRDS(nocontamAr,"rdata/ps/psA_deconfilt_relabun.rds")
saveRDS(nocontamBr,"rdata/ps/psB_deconfilt_relabun.rds")
saveRDS(nocontamA, "rdata/ps/psA_deconfilt_counts.rds")
saveRDS(nocontamB, "rdata/ps/psB_deconfilt_counts.rds")


# ---- Exploratory Plots ----
nocontamA %>%
  comp_barplot(tax_level = "class") +
  facet_wrap(~Sample_Type, scales = "free")

nocontamB %>%
  comp_barplot(tax_level = "class") +
  facet_wrap(~Sample_Type, scales = "free")

nocontamA %>%
  ps_filter(!str_detect(Treatment, "NC")) %>%
  ps_filter(!str_detect(Treatment, "PC")) %>%
  tax_transform("clr") %>%
  ord_calc() %>%
  ord_plot(color = "Sample_Type")

nocontamB %>%
  ps_filter(!str_detect(Treatment, "NC")) %>%
  tax_transform("clr") %>%
  ord_calc() %>%
  ord_plot(color = "Sample_Type")


# save work
save.image("rdata/comparealign.rdata")

