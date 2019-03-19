---
title: "HVC vs Shelf all animals"
author: "H. Na Choe"
date: "March 7, 2019"
output: html_document
---

Goal
   To find differentially expressed genes in the song nuclei of post hatch day 30 (PHD30) females treated with vehicle.
   Input all HVC and Shelf data generated from Kallisto pseudo-alignments
   Output three items:
       1) Differential expression of all genes significantly up or downregulated genes
         -> Genes with the same name are merged, due to unknown isoform status of genes in finch
       2) Shared differentially expressed genes across the different treatment groups/sex
       3) Pull all un-named entries whose differential expression is significant into a separate file to annotate back to the finch genome manually.




```{r setup, include=FALSE}
#Calling packages from Bioconductor and tidyverse
library(mixtools)
library(seqinr)
library(ape)
library(edgeR)
library(DESeq)
library(lattice)
library(MASS)
library(gplots)
library(calibrate)
library(plotly) 
library(tidyverse) 
library(ggrepel)

```

```{r}



```

```{r}

#Files have been sorted by sex+treatment. File names have their Library ID, the animal ID, and the region

#HVC first
Fem1_Veh_HVC <- read.table("C://Users/bhana/Desktop/RNA_seq/Female Veh/B7_HVC_abundance.tsv", header=T)
Fem2_Veh_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female Veh/D4_HVC_abundance.tsv", header=T)
Fem3_Veh_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female Veh/H5_HVC_abundance.tsv", header=T)

Fem4_E2_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_E2/B1_HVC_abundance.tsv", header=T)
Fem5_E2_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_E2/F3_HVC_abundance.tsv", header=T)
Fem6_E2_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_E2/F4_HVC_abundance.tsv", header=T)

Fem7_Exe_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_Exem/C5_HVC_abundance.tsv", header=T)
Fem8_Exe_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_Exem/F0_HVC_abundance.tsv", header=T)
Fem9_Exe_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_Exem/I2_HVC_abundance.tsv", header=T)

```
Mal1_Veh_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Veh/G1_HVC_abundance.tsv", header=T)
Mal2_Veh_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Veh/J1_HVC_abundance.tsv", header=T)
Mal3_Veh_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Veh/O4_HVC_abundance.tsv", header=T)

Mal4_E2_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_E2/B3_HVC_abundance.tsv", header=T)
Mal5_E2_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_E2/F1_HVC_abundance.tsv", header=T)
Mal6_E2_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_E2/F2_HVC_abundance.tsv", header=T)

Mal7_Exe_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Exem/C0_HVC_abundance.tsv", header=T)
Mal8_Exe_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Exem/H6_HVC_abundance.tsv", header=T)
Mal9_Exe_HVC <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Exem/I1_HVC_abundance.tsv", header=T)

```{r}
#Shelf second
Fem1_Veh_shelf <- read.table("C://Users/bhana/Desktop/RNA_seq/Female Veh/B7_shelf_abundance.tsv", header=T)
Fem2_Veh_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female Veh/D4_shelf_abundance.tsv", header=T)
Fem3_Veh_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female Veh/H5_shelf_abundance.tsv", header=T)

Fem4_E2_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_E2/B1_Shelf_abundance.tsv", header=T)
Fem5_E2_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_E2/F3_shelf_abundance.tsv", header=T)
Fem6_E2_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_E2/F4_shelf_abundance.tsv", header=T)

Fem7_Exe_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_Exem/C5_shelf_abundance.tsv", header=T)
Fem8_Exe_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_Exem/F0_shelf_abundance.tsv", header=T)
Fem9_Exe_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Female_Exem/I2_shelf_abundance.tsv", header=T)

```
Mal1_Veh_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Veh/G1_shelf_abundance.tsv", header=T)
Mal2_Veh_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Veh/J1_shelf_abundance.tsv", header=T)
Mal3_Veh_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Veh/O4_shelf_abundance.tsv", header=T)

Mal4_E2_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_E2/B3_shelf_abundance.tsv", header=T)
Mal5_E2_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_E2/F1_shelf_abundance.tsv", header=T)
Mal6_E2_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_E2/F2_shelf_abundance.tsv", header=T)

Mal7_Exe_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Exem/C0_shelf_abundance.tsv", header=T)
Mal8_Exe_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Exem/H6_shelf_abundance.tsv", header=T)
Mal9_Exe_shelf <- read.table("C://Users/bhana/Desktop/RNA_Seq/Male_Exem/I1_shelf_abundance.tsv", header=T)
```
```

```{r}
#Any Genscan IDs without refseq names are "NA" (first pass)
#"NA" are then given a unique numbered identifiers (example:unown1, unown2, unown3, ...)
#Contains isoforms/repeats of a genename since genscans are all unique

GenScanNames <- read_tsv("genscan_names.txt") 
RefAllNames <- read_tsv("Genscan_to_Refseq_all_indexing.tsv") 
JoinNames <- left_join(GenScanNames,RefAllNames, by = "target")
JoinNamesUniqueNA <- JoinNames
currentnum <- 1;
for (i in 1:length(JoinNamesUniqueNA$target)) {
  if (is.na(JoinNamesUniqueNA$gene[i])) { 
    JoinNamesUniqueNA$gene[i] <- paste('unown', currentnum, sep=""); 
    currentnum = currentnum + 1
  }
}

```

```{r}
#Merge isoforms from Genscan IDs to gene names from Human/Rat/Mouse/Finch/Chicken/Turkey refseq
#All counts associated with Genscan IDs are merged according to their refseq name or unown# name

MergeCounts <- function(abundanceFile){
  Abundance_to_Names <- left_join(JoinNamesUniqueNA, abundanceFile, by = c("target" = "target_id"));
  UniqueNAmes <- unique(Abundance_to_Names$gene);
  Result <- NULL;
  for (i in 1:length(UniqueNAmes)) {
    AllShared <- Abundance_to_Names[which(Abundance_to_Names$gene == UniqueNAmes [i]),]
    Sum <- round(apply(as.matrix(AllShared$est_counts), 2, sum))
    Row_to_add <- c(AllShared$gene[1], Sum)
    Result <- rbind(Result, Row_to_add)
  }
  colnames(Result) <- c("gene","est_counts")
  return(data.frame(Result))
}

```

```{r}

#Compile merged names/merged counts into a table for differential expression analysis by sex+treatment to compare between song nucleus to surround (Total 6 tables per nucleus)

#Change the names and input variable as needed for each group 
Combined_kallisto_output_FemVehHVC <- list("Fem1_Veh_HVC"=Fem1_Veh_HVC, "Fem2_Veh_HVC"=Fem2_Veh_HVC, "Fem3_Veh_HVC"=Fem3_Veh_HVC, "Fem1_Veh_shelf"=Fem1_Veh_shelf, "Fem2_Veh_shelf"=Fem2_Veh_shelf, "Fem3_Veh_shelf"=Fem3_Veh_shelf)

aggregate_counts <- function(ListAll_abundance_for_region) {
 countsTable <- MergeCounts(data.frame(ListAll_abundance_for_region[[1]]))
   for (i in 2:length(ListAll_abundance_for_region)) {
     abundanceFile = ListAll_abundance_for_region[[i]]
    countsTable <- left_join(countsTable, MergeCounts(data.frame(abundanceFile)), by = "gene")
   }
 colnames(countsTable) <- c("gene",names( ListAll_abundance_for_region) [1:length( ListAll_abundance_for_region)])
 #This is the only way I know how to reliably drop the "gene" column without losing information
 countsTable_tbl <- column_to_rownames(countsTable, var = "gene")
 countsTable_names <- rownames(countsTable_tbl)
 #This drops the rownames which has to be rescued with "countsTable_names"
 countsTable_tbl <- sapply(countsTable_tbl, function(x) as.numeric(as.character(x)))
 rownames(countsTable_tbl) <- countsTable_names
 return(countsTable_tbl)
}

#For testing in case things break
#tmp7 <- aggregate_counts(Combined_kallisto_output)
#tmp6test <- sapply(tmp6, function(x) as.numeric(as.character(x)))
#rownames(tmp6test) <- rownames(tmp6[0])


```

```{r}

#Make the aggregate count for Sex/treatmeant (batch 1/6)
Female_Veh_HVCshelf_counts <- aggregate_counts(Combined_kallisto_output_FemVehHVC)

```

```{r}

#make the labels, needs to match countsTable
diffExpression_labels <- data.frame(
  FileName = c("Fem1_Veh_HVC","Fem2_Veh_HVC","Fem3_Veh_HVC",
               "Fem1_shelf_HVC","Fem2_shelf_HVC","Fem3_shelf_HVC"), 
  Subject = c("A","B","C","A","B","C"), 
  Area = c(0,0,0,1,1,1))

```

```{r}

#Differential Expression function borrowed from Kevin/Yue/Ken
DiffExp <- function (Labels, Output_of_aggregation) {
  Area <- factor(Labels$Area);
  Subject <- factor(Labels$Subject);
  design <- model.matrix(~Subject+Area)
  cds <- newCountDataSet(Output_of_aggregation,Area);
  cds <- estimateSizeFactors(cds);
  cds <- estimateDispersions(cds);
  d <- nbinomTest(cds,"0","1")
  e.litter <- DGEList(counts=Output_of_aggregation)
  e.litter <- estimateGLMCommonDisp(e.litter,design)
  e.litter <- estimateGLMTrendedDisp(e.litter,design)
  e.litter <- estimateGLMTagwiseDisp(e.litter,design)
  fit <- glmFit(e.litter, design);
  lrt <- glmLRT(fit);
  diff <- topTags(lrt,n=dim(lrt)[1])$table
  result <- merge(merge(diff,Output_of_aggregation,by=0,sort=F),d, by.x="Row.names", by.y="id",sort=F)
  colnames(result)[1] <- "id"
  return(result)
}


Diff_FemVeh_HVC <- DiffExp(diffExpression_labels, Female_Veh_HVCshelf_counts)


```




## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
