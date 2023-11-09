# Clean AMR++ Output : v2 and v3 

# ---- load packages
library(tidyverse)
library(tidylog)
library(reshape2)
library(ggplot2)

#read in tables
t2 <- read.csv("raw/metagenomic/resistome/app-out/AMR_analytic_matrix-v2.csv")
t3 <- read.csv("raw/metagenomic/resistome/app-out/AMR_analytic_matrix-v3.csv")
v2 <- read.csv("tables/meta/megares_full_annotations_v2.00.csv")
v3 <- read.csv("tables/meta/megares_annotations_v3.00.csv")
met <- read.csv("tables/meta/turkey_metadata.csv")

#add the controls to met
met <- rbind(met %>%
  filter(Project_ID %in% colnames(t2)), data.frame(Project_ID = c("T780", "T781"),
                                                   Sample_Type = c("DNA_PC", "DNA_NC"),
                                                   Animal_ID = c("DNA_PC", "DNA_NC"),
                                                   Timepoint = c(NA, NA),
                                                   Treatment = c("DNA_PC", "DNA_NC")))


#edit rownames
t2 <- cbind(as.data.frame(t2), as.data.frame(str_split_fixed(rownames(as.data.frame(t2)), "\\|", 5))[5])%>% #extract gene name from amr++ output and add to count matrix
  aggregate(.~V5, FUN= sum) #sum counts by gene

t3 <- cbind(as.data.frame(t3), as.data.frame(str_split_fixed(rownames(as.data.frame(t3)), "\\|", 5))[5])%>% #extract gene name from amr++ output and add to count matrix
  aggregate(.~V5, FUN= sum) #sum counts by gene

rownames(t2) <- t2$V5
rownames(t3) <- t3$V5

t2 <- t2[,2:ncol(t2)]
t3 <- t3[,2:ncol(t3)]

#melt 
t2m <- reshape2::melt((as.matrix(t2)))
colnames(t2m) <- c("Gene","Project_ID", "Counts")

write.table(t2m, "tables/amr/v2meltdf.txt", sep = "\t")

t3m <- reshape2::melt((as.matrix(t3)))
colnames(t3m) <- c("Gene","Project_ID", "Counts")

write.table(t3m, "tables/amr/v3meltdf.txt", sep = "\t")

v3genes <- t3m %>%
  filter(Counts!=0)%>%
  select(Gene)%>%
  unique()

v2genes <- t2m %>%
  filter(Counts!=0)%>%
  select(Gene)%>%
  unique()

diff <- v3genes %>%
  filter(!Gene %in% v2genes$Gene )

write.table(diff, "tables/amr/v3-v2compare.txt", sep = "\t")

save.image("rdata/compareversions.rdata")


