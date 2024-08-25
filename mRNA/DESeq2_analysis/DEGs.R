# author: Badreldin Afifi
# 12.03.2023
# root-knot nematode (Meloidogyne incognita)-infected tomato (Solanum lycopersicum) roots 

library(dplyr) 
library(limma)
library(edgeR )
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(apeglm)

# single read counts
R1 <- read.table("HISAT2_pair_end.tabular", header = T)
rownames(R1) <- R1$Geneid
R1 <- R1[-1] # rename the rownames

# Pair read counts
R2 <- read.table("HISAT2_single_end.tabular", header = T)
rownames(R2) <- R2$Geneid
R2 <- R2[-1]

# combind
Allmatrix <- cbind(R1,R2)
 
#order the columns
ctrlnamess <- c("SRR4423503","SRR4423504","SRR4423507","SRR4423508",
                "SRR4423511","SRR4423512", "SRR4423515","SRR4423516",
                "SRR4423519","SRR4423520")

Infnamess <- c("SRR4423505","SRR4423506","SRR4423509","SRR4423510",
               "SRR4423513","SRR4423514","SRR4423517","SRR4423518",
               "SRR4423521","SRR4423522")

cotrl_data <- Allmatrix %>% select(all_of(ctrlnamess))
Infec_data <- Allmatrix %>% select(all_of(Infnamess))

# combine all data
samples <- cbind(Infec_data, cotrl_data)

# remove all 0 data

samples <- samples[ rowSums(samples) > 50,]
head(samples)


# names aggregation
{
id <- sub("\\..*", "",row.names(samples))
samples <- cbind(samples, id)

# check duplciation of of gene symbols?  
x=duplicated(samples$id)  
sum(x)

#solutions : aggregation
agg.samples <- aggregate(samples, list(samples$id),FUN=mean)

rownames(agg.samples)=agg.samples$Group.1
agg.samples=agg.samples[- 1]
agg.samples=agg.samples[- 21]
}

write.table(samples, file = "HISAT2_matrix.txt",row.names = T,col.names = T,quote = F)

#====================[Explore The Data]=========================================
samples <- read.table("HISAT2_matrix.txt")
#***For your own data, you will have to decide which columns will be useful in the analysis
samplnames <- colnames(samples)

## source_name_ch1 seems to contain factors we might need for the analysis. Let's pick just those columns
condition <- c(rep("Infected", 10),rep("Control", 10))
stage <- c(rep("stage1", 2),rep("stage2", 2),rep("stage3", 2),rep("stage4", 2),rep("stage5", 2),
           rep("stage1", 2),rep("stage2", 2),rep("stage3", 2),rep("stage4", 2),rep("stage5", 2))

sampleInfo <- as.data.frame(bind_cols(samplnames,condition,stage)) # meta data
names(sampleInfo)[1]="sample_names" 
names(sampleInfo)[2]="sample_source"
names(sampleInfo)[3]="stage"

#***Unsupervised analysis is a good way to get an understanding of the sources of variation in the data.
#It can also identify potential outlier samples.

corMatrix <- cor(samples,use="c") ## argument use="c" stops an error if there are any missing data points
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,annotation_col=sampleInfo)    

#**** A complementary approach is to use Principal Components Analysis (PCA).
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
pca <- prcomp(t(samples),)
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=sample_source,label= sample_names)) + 
  geom_point() +
  geom_text(hjust = 0, nudge_x = 0.05)

#==============================[DESeq2]=============================================
#Design specifies how the counts from each gene depend on our variables in "sampleInfo"
#For this dataset the factor we care about is our infection status (sample_source)
#tidy=TRUE argument, which tells DESeq2 to output the results table with rownames as a first #column called 'row.

head(sampleInfo)

# Stage 1
{
st1_samp <- samples[,c(1:2, 11:12)]
st1_Info <- sampleInfo[c(1:2, 11:12),]

#Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData= round(st1_samp), 
                              colData= st1_Info, 
                              design= ~  sample_source) # tidy = F
dds <- DESeq(dds)

#             *****Results*******
resultsNames(dds)
res <- results(dds, cooksCutoff= T) #coef=1
res <- res[order(res$padj),] #Sort summary list by p-value
summary(res) #summary of results


# changing p-adj and LFC values 
#sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1 , na.rm=TRUE)
#results(dds, alpha=0.05, lfcThreshold = 1)

res05 <- subset(res, res$padj < 0.05 & abs(res$log2FoldChange) > 1)
res05 <- res05[order(res05$padj),] #Sort summary list by p-value
summary(res05)

write.table(res, file = "stage1_Deseq2.txt",row.names = T,col.names = T,quote = F) 
write.table(res05, file = "stage1_DEGs.txt",row.names = T,col.names = T,quote = F) 
write.table(rownames(res05), file = "stage1_DEGs[names].txt",row.names = T,col.names = T,quote = F) 


#               **** Volcano plot *******

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Stage 1", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj< 0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, padj< 0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

}

# Stage 2
{
  st_samp <- samples[,c(3:4, 13:14)]
  st_Info <- sampleInfo[c(3:4, 13:14),]
  
  #Construct DESEQDataSet Object
  dds <- DESeqDataSetFromMatrix(countData= round(st_samp), 
                                colData= st_Info, 
                                design= ~  sample_source) # tidy = F
  dds <- DESeq(dds)
  
  #             *****Results*******
  resultsNames(dds)
  res <- results(dds, cooksCutoff= T) #coef=1
  res <- res[order(res$padj),] #Sort summary list by p-value
  summary(res) #summary of results
  
  
  # changing p-adj and LFC values 
  #sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1 , na.rm=TRUE)
  #results(dds, alpha=0.05, lfcThreshold = 1)
  
  res05 <- subset(res, res$padj < 0.05 & abs(res$log2FoldChange) > 1)
  res05 <- res05[order(res05$padj),] #Sort summary list by p-value
  res05up <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
  res05down <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
  ummary(res05)
  
  write.table(res, file = "stage2_Deseq2.txt",row.names = T,col.names = T,quote = F) 
  write.table(res05, file = "stage2_DEGs.txt",row.names = T,col.names = T,quote = F) 
  write.table(rownames(res05), file = "stage2_DEGs[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05up), file = "stage2_DEGsup[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05down), file = "stage2_DEGsdown[names].txt",row.names = F,col.names = F,quote = F) 
  
  
  #               **** Volcano plot *******
  
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Stage 2", xlim=c(-6,6)))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj< 0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, padj< 0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  
}

# Stage 3
{
  st_samp <- samples[,c(5:6, 15:16)]
  st_Info <- sampleInfo[c(5:6, 15:16),]
  
  #Construct DESEQDataSet Object
  dds <- DESeqDataSetFromMatrix(countData= round(st_samp), 
                                colData= st_Info, 
                                design= ~  sample_source) # tidy = F
  dds <- DESeq(dds)
  
  #             *****Results*******
  resultsNames(dds)
  res <- results(dds, cooksCutoff= T) #coef=1
  res <- res[order(res$padj),] #Sort summary list by p-value
  summary(res) #summary of results
  
  
  # changing p-adj and LFC values 
  #sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1 , na.rm=TRUE)
  #results(dds, alpha=0.05, lfcThreshold = 1)
  
  res05 <- subset(res, res$padj < 0.05 & abs(res$log2FoldChange) > 1)
  res05 <- res05[order(res05$padj),] #Sort summary list by p-value
  res05up <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
  res05down <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
  summary(res05)
  
  write.table(res, file = "stage3_Deseq2.txt",row.names = T,col.names = T,quote = F) 
  write.table(res05, file = "stage3_DEGs.txt",row.names = T,col.names = T,quote = F) 
  write.table(rownames(res05), file = "stage3_DEGs[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05up), file = "stage3_DEGsup[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05down), file = "stage3_DEGsdown[names].txt",row.names = F,col.names = F,quote = F) 
  
  
  #               **** Volcano plot *******
  
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Stage 3", xlim=c(-6,6)))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj< 0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, padj< 0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  
}

# Stage 4
{
  st_samp <- samples[,c(7:8, 17:18)]
  st_Info <- sampleInfo[c(7:8, 17:18),]
  
  #Construct DESEQDataSet Object
  dds <- DESeqDataSetFromMatrix(countData= round(st_samp), 
                                colData= st_Info, 
                                design= ~  sample_source) # tidy = F
  dds <- DESeq(dds)
  
  #             *****Results*******
  resultsNames(dds)
  res <- results(dds, cooksCutoff= T) #coef=1
  res <- res[order(res$padj),] #Sort summary list by p-value
  summary(res) #summary of results
  
  
  # changing p-adj and LFC values 
  #sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1 , na.rm=TRUE)
  #results(dds, alpha=0.05, lfcThreshold = 1)
  
  res05 <- subset(res, res$padj < 0.05 & abs(res$log2FoldChange) > 1)
  res05 <- res05[order(res05$padj),] #Sort summary list by p-value
  res05up <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
  res05down <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
  summary(res05)
  
  write.table(res, file = "stage4_Deseq2.txt",row.names = T,col.names = T,quote = F) 
  write.table(res05, file = "stage4_DEGs.txt",row.names = T,col.names = T,quote = F) 
  write.table(rownames(res05), file = "stage4_DEGs[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05up), file = "stage4_DEGsup[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05down), file = "stage4_DEGsdown[names].txt",row.names = F,col.names = F,quote = F) 
  
  
  #               **** Volcano plot *******
  
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Stage 4", xlim=c(-6,6)))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj< 0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, padj< 0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  
}

# Stage 5
{
  st_samp <- samples[,c(9:10, 19:20)]
  st_Info <- sampleInfo[c(9:10, 19:20),]
  
  #Construct DESEQDataSet Object
  dds <- DESeqDataSetFromMatrix(countData= round(st_samp), 
                                colData= st_Info, 
                                design= ~  sample_source) # tidy = F
  dds <- DESeq(dds)
  
  #             *****Results*******
  resultsNames(dds)
  res <- results(dds, cooksCutoff= T) #coef=1
  res <- res[order(res$padj),] #Sort summary list by p-value
  summary(res) #summary of results
  
  
  # changing p-adj and LFC values 
  #sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1 , na.rm=TRUE)
  #results(dds, alpha=0.05, lfcThreshold = 1)
  
  res05 <- subset(res, res$padj < 0.05 & abs(res$log2FoldChange) > 1)
  res05 <- res05[order(res05$padj),] #Sort summary list by p-value
  res05up <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
  res05down <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
  summary(res05)
  
  write.table(res, file = "stage5_Deseq2.txt",row.names = T,col.names = T,quote = F) 
  write.table(res05, file = "stage5_DEGs.txt",row.names = T,col.names = T,quote = F) 
  write.table(rownames(res05), file = "stage5_DEGs[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05up), file = "stage5_DEGsup[names].txt",row.names = F,col.names = F,quote = F) 
  write.table(rownames(res05down), file = "stage5_DEGsdown[names].txt",row.names = F,col.names = F,quote = F) 
  
  
  #               **** Volcano plot *******
  
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Stage 5", xlim=c(-6,6)))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj< 0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  with(subset(res, padj< 0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
  
}
