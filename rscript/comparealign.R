# Compare AMR: BLAST vs AMR++v2
# Sophia Kenney - 9 November 2023

#load packages
library(tidyverse)
library(tidylog)
library(phyloseq)


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

#filter out snp confirmation files
blast <- blast %>%
  filter(!str_detect(V4, "SNPConfirmation"))

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

# ---- Decontam ----




# save work
save.image("rdata/comparealign.rdata")

