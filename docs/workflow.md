# RNA-Seq Differential Expression Pipeline

```bash
## Quality control of FASTQ files

fastqc *.fastq.gz -o ../results/qc/

## Align reads using Hisat2

hisat2 -x genome_index -U sample.fastq.gz -S output.sam

## Convert SAM to BAM

samtools view -bS output.sam > output.bam

## Sort BAM

samtools sort output.bam -o output_sorted.bam

## Index BAM

samtools index output_sorted.bam
```

```R
library(DESeq2)
library(pheatmap)
library(qqman)

## Load sample metadata

sampleInfo <- read.table("assignmentsampledata.txt", header=TRUE)
row.names(sampleInfo) <- sampleInfo$sampleName

## Create DESeq dataset

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleInfo,
  directory = "data/",
  design = ~ condition
)

## Run DE analysis

dds <- DESeq(dds)

## Normalised counts

norm_counts <- counts(dds, normalized=TRUE)

## PCA

rld <- rlog(dds)
pca <- prcomp(t(assay(rld)))

plot(pca$x[,1], pca$x[,2],
     col = as.integer(colData(rld)$condition),
     main="PCA Plot")

## Differential expression

res <- results(dds)

## Volcano plot

log10.pval <- -log10(res$padj)
log2.fc <- res$log2FoldChange

plot(log2.fc, log10.pval,
     xlab="log2 Fold Change",
     ylab="-log10 p-value",
     main="Volcano Plot")

## Heatmap

topGenes <- res$padj < 0.05
pheatmap(assay(rld)[which(topGenes),],
         scale="row",
         show_rownames=FALSE)

set.seed(42)

## Load required libraries (install separately if needed)

library(DESeq2)
library(topGO)
library(org.Hs.eg.db)

## Define contrast

res <- results(dds, contrast = c(factor, treatment, reference))

## Order results by adjusted p-value

res <- res[order(res$padj), ]

## Define significance threshold

padj_threshold <- 0.05
topGenes <- res$padj < padj_threshold

## Create output directory

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

## Save all results

write.table(
res,
file = paste0("results/tables/diff_genes_", treatment, "_vs_", reference, ".txt"),
sep = "\t",
quote = FALSE,
row.names = TRUE
)

## Save significant genes only

write.table(
res[topGenes, ],
file = paste0("results/tables/sig_diff_genes_", treatment, "_vs_", reference, ".txt"),
sep = "\t",
quote = FALSE,
row.names = TRUE
)

## Save normalised counts

write.table(
countsnorm,
file = "results/tables/normalised_counts.txt",
sep = "\t",
quote = FALSE,
row.names = TRUE
)

## Remove NA values

res_clean <- res[!is.na(res$padj), ]
topGenes_clean <- topGenes[!is.na(res$padj)]

## Create gene list (1 = significant, 0 = not)

geneList <- factor(as.integer(topGenes_clean))
names(geneList) <- rownames(res_clean)

## Map GO terms

allGO2genes <- annFUN.org(
whichOnto = "BP",
mapping = "org.Hs.eg.db",
ID = "symbol"
)

## Create topGO dataset

GOdata <- new(
"topGOdata",
ontology = "BP",
allGenes = geneList,
annot = annFUN.GO2genes,
GO2genes = allGO2genes
)

## Run Fisher test

go_results <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

## Generate GO results table

go_table <- GenTable(
GOdata,
fisher = go_results,
orderBy = "fisher",
topNodes = length(score(go_results))
)

# Save GO results

write.csv(
go_table,
file = paste0("results/tables/GO_enrichment_", treatment, "_vs_", reference, ".csv"),
row.names = FALSE
)
```
