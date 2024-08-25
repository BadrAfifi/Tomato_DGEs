library(GEOquery)
library(hgu133plus2.db)
library(hgu133acdf)
library(limma)
library(hgu133plus2cdf)
library(GSEABase)
library(GOstats)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(curl)
library(RCurl)
library(affy)
library(readr)
library(hgu133a.db)
library(genefilter)
library(multtest)
library(affyPLM)
library(pheatmap)
library(pacman)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(Biobase)
library(stats)
library(DESeq2)
library(edgeR)
# Load count files
control_counts <- read.table("unifs1.txt")
treated_counts <- read.table("infs1.txt")

All <- as.data.frame(c(control_counts,treated_counts))
All <- All[-3]

rownames(All) <- All$V1
All <- All[-1]

names(All)[1] <- "V1"
names(All)[2] <- "V2"

All
All1 <- All

All2 <- as.data.frame(c(All,All1))
rownames(All2) <- row.names(All)

x <- as.matrix(All2)
#==============================================================


# 2. Construct DGEList and add group information at the same time
sampleinfo <- data.frame(
  sample = colnames(x),
  condition = c("ctrl", "treat", "ctrl", "treat"),
  batch = c("I", "I", "II", "II")
)#Set batch, remove batch effect

coldata <- data.frame(row.names = sampleinfo$sample, 
                      condition = factor(sampleinfo$condition, levels = c("ctrl","treat")), 
                      batch = factor(sampleinfo$batch))


y <- DGEList(counts=x, group=coldata$condition) #Note that the difference comparison direction is the latter of condition
#- the former, that is, treat-ctrl, which is opposite to DESeq2

# 3. Filter low expression genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# 4. TMM standardization
y <- calcNormFactors(y)

# 5.Estimating dispersions 
# Batch effects can be removed here 
#design <- model.matrix(~ coldata$batch + coldata$condition)

design <- model.matrix(~ coldata$condition)
y <- estimateDisp(y,design)

# 6. Calculate the difference
## 6.1 exact test
#et <- exactTest(y)
#res <- topTags(et,n=30000)$table
#

## 6.2 quasi-likelihood (QL) F-test, this method is less strict than the likelihood ratio test
#fit <- glmQLFit(y,design)
#qlf <- glmQLFTest(fit,coef=2)
#res <- topTags(qlf,n=30000)$table

# 6.3 likelihood ratio test
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
res <- topTags(lrt,n=30000)$table

res$gene_id <- rownames(res)


# 7. Screen differentially expressed genes

p <- log10(0.05) * (-1)
fc <- log2(2)
degs <- res %>%
  select(gene_id,logFC,logCPM,LR,PValue,FDR) %>%  #
  mutate(logp = log10(PValue) * (-1)) %>%
  mutate(type = if_else(logp > p, if_else(logFC > fc, "Up", if_else(logFC < (-fc), "Down", "N.s")), "N.s")) %>%
  arrange(type)

table(degs$type)
write.table(degs, "degs_edgeR.txt", row.names = F, sep='\t')

#[2] Visualization
# [A] ----------- Volcano (UP/Down)-----------  

degs %>%
  filter(logFC < 10) %>%
  filter(logp < 50) %>%
  ggplot(aes(logFC, logp)) +
  geom_point(aes(colour = type), alpha = 1, size = 2) +
  scale_colour_manual(values = c("#00FF00", "grey", "#FF0000")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
   geom_hline(yintercept = (0.05), linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = seq(-4, 4, 2), labels = seq(-4, 4, 2), limits = c(-5, 5)) +
  labs(x = quote(log[2] ~ FoldChange), y = quote(-log[10] ~ pvalue), colour = "") +
  theme_light() # +
  #theme( axis.ticks = element_line(colour = "black"),
    #axis.text = element_text(colour = "black", size = 14),
    #axis.title = element_text(size = 14),
    #legend.key.height = unit(1, "cm"), legend.key.width = unit(1, "cm"))




