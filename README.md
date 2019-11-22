# Helicobacter-human dual-species transcriptomics in vitro and in vivo

# Table of Contents

1. [Pre-analysis setup](#preanalysis)
    1. [Set bash software paths](#preanalysis_bashsoftware)
2. [Confirm cag KO in H. pylori](#confirmcagko)
    1. [Download 4 H. pylori cagKO fastq files from the SRA](#confirmcagko_sra)
    2. [Align H. pylori cagKO fastq files to H. pylori reference genome](#confirmcagko_align)
    3. [Find reads that mapped to H. pylori](#confirmcagko_findhpylorireads)
    4. [Subset fastqs to only contain H. pylori mapping reads](#confirmcagko_assemble)
3. [In vivo human RNA-Seq analysis](#invivohuman)
    1. [Quantify human genes from in vivo RNA-Seq samples](#invivohuman_quant)
    2. [Create human counts and TPM dataframes](#invivohuman_countstpm)
    3. [Choose top candidate genes for IPA](#invivohuman_filtergenes)
    

# Pre-analysis <a name="preanalysis"></a>
## Set project-wide bash software and directory paths <a name="preanalysis_bashsoftware"></a>

```{bash, eval = F}
THREADS=4

BOWTIE2_BIN_DIR=/usr/local/packages/bowtie2-2.3.4.3
KALLISTO_BIN_DIR=/usr/local/packages/kallisto-0.45.0
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin

REFERENCES_DIR=/local/aberdeen2rw/julie/Matt_dir/EHPYL/references/
```

## Download references and create combined human-H. pylori reference <a name="preanalysis_references"></a>

```{bash, eval = F}
wget -O "$REFERENCES_DIR"/hpylori26995.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/525/GCF_000008525.1_ASM852v1/GCF_000008525.1_ASM852v1_genomic.fna.gz
wget -O "$REFERENCES_DIR"/hsapiensGRCh38.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
gunzip "$REFERENCES_DIR"/hpylori26995.fna.gz
gunzip "$REFERENCES_DIR"/hsapiensGRCh38.fna.gz
cat "$REFERENCES_DIR"/hpylori26995.fna "$REFERENCES_DIR"/hsapiensGRCh38.fna > "$REFERENCES_DIR"/combined_hsapiensGRCh38_hpylori26695.fna


```


# Confirm cag KO in H. pylori <a name="confirmcagko"></a>
## Set bash input and directory paths <a name="confirmcagko_setpaths"></a>

```{bash, eval = F}
SRR_ID_LIST=/local/projects-t3/EBMAL/mchung_dir/EHPYL/genomic_srr.list
READS_DIR=/local/scratch/mchung
REFERENCES_DIR=/local/aberdeen2rw/julie/Matt_dir/EHPYL/references/
OUTPUT_DIR=
```


### Download 4 H. pylori cagKO fastq files from the SRA (2019-06-25) <a name="confirmcagko_sra"></a>

#### Commands:
```{bash, eval = F}
while read SRR_ID
do
  qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" --gzip -O "$READS_DIR"
done < "$SRR_ID_LIST"
```

### Align H. pylori cagKO fastq files to H. pylori reference genome (2019/06/25) <a name="confirmcagko_align"></a>

```{bash, eval = F}
"$BOWTIE2_BIN_DIR"/bowtie2-build --large-index "$REF_FNA" "$REF_FNA"
while read SRR
do
echo ""$BOWTIE2_BIN_DIR"/bowtie2 --threads "$THREADS" -x "$REF_FNA" -1 "$READS_DIR"/"$SRR"_1.fastq.gz -2 "$READS_DIR"/"$SRR"_2.fastq.gz | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=5G -N bowtie2 -wd "$OUTPUT_DIR"
done < "$SRR_LIST"
```

```{bash, eval = F}
SRR_LIST=/local/projects-t3/EBMAL/mchung_dir/EHPYL/genomic_srr.list
THREADS=16
REF_FNA=/local/aberdeen2rw/julie/Matt_dir/EHPYL/references/NC_000915.1.fa
FASTQ_DIR=/local/scratch/mchung/
OUTPUT_DIR=/local/aberdeen2rw/julie/Matt_dir/EHPYL/genomic_bam
```

### Find reads that mapped to H. pylori (2019/06/25) <a name="confirmcagko_findhpylorireads"></a>

```{bash, eval = F}
for BAM in $(find "$BAM_DIR" -name "*[.]bam")
do
  "$SAMTOOLS_BIN_DIR"/samtools view "$BAM" | grep "$query" | cut -f1 | sort -n | uniq > "$BAM".subset.reads &
done
```

```{bash, eval = F}
QUERY=NC_000915.1
BAM_DIR=/local/aberdeen2rw/julie/Matt_dir/EHPYL/genomic_bam
```

## Subset fastqs to only contain H. pylori mapping reads (2019/06/26) <a name="confirmcagko_makehpylorimapfastq"></a>

```{bash, eval = F}
while read SRR
do
echo -e "/usr/local/packages/seqtk-1.2/bin/seqtk subseq "$FASTQ_DIR"/"$SRR"_1.fastq.gz "$SUBSET_READS_LIST_DIR"/"$SRR".bam.subset.reads > "$FASTQ_DIR"/"$SRR"_1.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=5G -N seqtk -wd "$FASTQ_DIR"
echo -e "/usr/local/packages/seqtk-1.2/bin/seqtk subseq "$FASTQ_DIR"/"$SRR"_2.fastq.gz "$SUBSET_READS_LIST_DIR"/"$SRR".bam.subset.reads > "$FASTQ_DIR"/"$SRR"_2.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=5G -N seqtk -wd "$FASTQ_DIR"
done < "$SRR_LIST"
```

```{bash, eval = F}
FASTQ_DIR=/local/scratch/mchung
SUBSET_READS_LIST_DIR=/local/aberdeen2rw/julie/Matt_dir/EHPYL/genomic_bam
SRR_LIST=/local/projects-t3/EBMAL/mchung_dir/EHPYL/genomic_srr.list
```

## Assemble H. pylori genomes using H. pylori-mapped fastqs (2019/06/26) <a name="confirmcagko_assemble"></a>

```{bash, eval = F}
while read SRR
do
mkdir "$OUTPUT_DIR"/"$SRR"
qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N spades -wd "$OUTPUT_DIR" -b y /usr/local/packages/spades-3.13.0/bin/spades.py -1 "$FASTQ_DIR"/"$SRR"_1.subset.fastq -2 "$FASTQ_DIR"/"$SRR"_2.subset.fastq -o "$OUTPUT_DIR"/"$SRR"
done < "$SRR_LIST"
```

```{bash, eval = F}
SRR_LIST=/local/projects-t3/EBMAL/mchung_dir/EHPYL/genomic_srr.list
THREADS=16
FASTQ_DIR=/local/scratch/mchung/
OUTPUT_DIR=/local/aberdeen2rw/julie/Matt_dir/EHPYL/genomic_spades
```

## In vivo human RNA-Seq analysis <a name="invivohuman"></a>
### Quantify human genes from in vivo RNA-Seq samples (2019/10/14) <a name="invivohuman_quant"></a>
```{bash, eval = F}
for FASTQ in $(find "$FASTQ_DIR"/*R/ILLUMINA_DATA/ -name *"fastq.gz" | sed "s/_R[12].fastq.gz//g" | sort -n | uniq)
do

  SAMPLE="$(echo "$(dirname "$FASTQ")" | sed "s/\\/ILLUMINA_DATA.*//g" | sed "s/.*\\///g")"
  mkdir "$OUTPUT_DIR"/"$SAMPLE"
  
  echo -e "qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N kallisto -wd "$OUTPUT_DIR"/"$SAMPLE" -b y "$KALLISTO_BIN_DIR"/kallisto quant -t "$THREADS" -i "$NUC_TRANSCRIPT_FNA".kallisto.index -o "$OUTPUT_DIR"/"$SAMPLE"/"$(basename "$FASTQ")" "$FASTQ"_R1.fastq.gz "$FASTQ"_R2.fastq.gz"
  
done
```

```{bash, eval = F}
NUC_TRANSCRIPT_FNA=/local/aberdeen2rw/julie/Matt_dir/EHPYL/references/combined_hsapiensGRCh38_hpylori26695.cds.fna
FASTQ_DIR=/local/aberdeen2ro/ESTAD
OUTPUT_DIR=/local/projects-t3/EBMAL/mchung_dir/ESTAD/kallisto
THREADS=4
```

### Create counts and TPM dataframes (2019/10/16) <a name="invivohuman_countstpm"></a>
```{r}
counts.dir <- "Z:/EBMAL/mchung_dir/ESTAD/kallisto"
output.dir <- "Z:/EBMAL/mchung_dir/ESTAD/analysis"

counts.files <- list.files(counts.dir, recursive = T, pattern = "abundance.tsv")

samples <- unique(gsub("\\/.*","",counts.files))

counts <- as.data.frame(matrix(0,
                               nrow=nrow(read.delim(paste0(counts.dir,"/",counts.files[1]))),
                               ncol=length(samples)))
rownames(counts) <- read.delim(paste0(counts.dir,"/",counts.files[1]))[,1]
colnames(counts) <- samples

genelength <- read.delim(paste0(counts.dir,"/",counts.files[1]))[,3]

for(i in 1:length(counts.files)){
  counts.file <- read.delim(paste0(counts.dir,"/",counts.files[i]))
  sample <- gsub("\\/.*","",counts.files[i])
  
  counts[,which(colnames(counts) == sample)] <- counts[,which(colnames(counts) == sample)] + counts.file[match(rownames(counts),counts.file[,1]),4]
}

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(counts,
            paste0(output.dir,"/human_invivo_counts.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
write.table(tpm,
            paste0(output.dir,"/human_invivo_tpm.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

### Choose top candidate genes for IPA (2019/10/16) <a name="invivohuman_filtergenes"></a>

IPA allows only 3000 genes for core analysis. The 190,198 human transcripts in the analysis were pared down based on two criteria: (1) the transcript must have a CPM value in all samples greater than 2000 reads in the lowest sequenced sample and (2) the gene must have a log2TPM ratio value >2 or <-2 in at least one sample. After both filters are applied, the log2TPM ratio of the 2,697 genes that met both criteria were used as an input to IGV.

```{r}
library(edgeR)

cpm.cutoff <- 2000/min(colSums(counts)) * 1000000

y <- DGEList(counts = counts, group = colnames(counts))
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= 1

tpm <- tpm[keep,]
ratiotpm <- tpm[rowSums(tpm) != 0,]
ratiotpm <- as.data.frame(t(apply(ratiotpm,1,function(x){as.numeric(as.character(x))/mean(as.numeric(as.character(x)))})))
colnames(ratiotpm) <- samples

log2ratiotpm <- log2(ratiotpm)
log2ratiotpm <- log2ratiotpm[rowSums(log2ratiotpm < -2 | log2ratiotpm > 2) > 0,]

write.table(log2ratiotpm,
            paste0(output.dir,"/human_invivo_log2ratiotpm.tsv"),,
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

---
title: "Untitled"
author: "Matthew Chung"
date: "October 17, 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
output.dir <- "Z:/EBMAL/mchung_dir/ESTAD/analysis"
```

## R Markdown


```{r,fig.height=2,fig.width=7}
require(ggplot2)
require(gridExtra)
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(20)

hmcolkey.plot <- ggplot() + 
  geom_raster(aes(x=seq(-5,5,0.5), y=seq(-5,5,0.5), fill = seq(-5,5,0.5)))+
  scale_fill_gradientn(name = "activation z-score",
                       colours=hmcol,
                       breaks=c(-5,0,5))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolkey.plot)
pdf(paste0(output.dir,"/hmcolorkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()
```

```{r,fig.height=7,fig.width=4}
require(gplots)
require(varhandle)

canonicalpathway.raw.path <- "Z:/EBMAL/mchung_dir/ESTAD/analysis/comparison_canonicalpathway.txt"
canonicalpathway.raw <- read.delim(canonicalpathway.raw.path,
                                   comment.char = "Â",
                                   header = F,
                                   col.names = seq(1,5))
canonicalpathway <- canonicalpathway.raw
rownames(canonicalpathway) <- canonicalpathway[,1]
colnames(canonicalpathway) <- unfactor(canonicalpathway[1,])
canonicalpathway <- canonicalpathway[-1,-1]
canonicalpathway <- apply(canonicalpathway,c(1,2),function(x){as.numeric(as.character(x))})

canonicalpathway <- canonicalpathway[1:30,]

pdf(paste0(output.dir,"/hm_canonicalpathway.pdf"),
    height=7,
    width=4)
heatmap.2(as.matrix(canonicalpathway),
          col=hmcol,
          trace="none",
          Rowv = F,
          Colv = F,
          lwid=c(1,4),
          lhei = c(1,5),
          breaks = seq(-5,5,by=.5),
          dendrogram = "none")
dev.off()
```

```{r,fig.height=11,fig.width=4}
require(gplots)
diseasesfunction.raw.path <- "Z:/EBMAL/mchung_dir/ESTAD/analysis/comparison_diseasesfunctions.txt"
diseasesfunction.raw <- read.delim(diseasesfunction.raw.path,
                                   comment.char = "Â",
                                   header = F,
                                   col.names = seq(1,5))
diseasesfunction <- diseasesfunction.raw
rownames(diseasesfunction) <- diseasesfunction[,1]
colnames(diseasesfunction) <- unfactor(diseasesfunction[1,])
diseasesfunction <- diseasesfunction[-1,-1]
diseasesfunction <- apply(diseasesfunction,c(1,2),function(x){as.numeric(as.character(x))})

diseasesfunction <- diseasesfunction[1:30,]

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(20)

pdf(paste0(output.dir,"/hm_diseasesfunction.pdf"),
    height=7,
    width=4)
heatmap.2(as.matrix(diseasesfunction),
          col=hmcol,
          trace="none",
          Rowv = F,
          Colv = F,
          lwid=c(1,4),
          lhei = c(1,5),
          breaks = seq(-5,5,by=.5),
          dendrogram = "none")
dev.off()
```

```{r,fig.height=7,fig.width=4}
require(gplots)
upstreamanalysis.raw.path <- "Z:/EBMAL/mchung_dir/ESTAD/analysis/comparison_upstreamanalysis.txt"
upstreamanalysis.raw <- read.delim(upstreamanalysis.raw.path,
                                   comment.char = "Â",
                                   header = F,
                                   col.names = seq(1,5))
upstreamanalysis <- upstreamanalysis.raw
rownames(upstreamanalysis) <- upstreamanalysis[,1]
colnames(upstreamanalysis) <- unfactor(upstreamanalysis[1,])
upstreamanalysis <- upstreamanalysis[-1,-1]
upstreamanalysis <- apply(upstreamanalysis,c(1,2),function(x){as.numeric(as.character(x))})

upstreamanalysis <- upstreamanalysis[1:30,]

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(20)

pdf(paste0(output.dir,"/hm_upstreamanalysis.pdf"),
    height=7,
    width=4)
heatmap.2(as.matrix(upstreamanalysis),
          col=hmcol,
          trace="none",
          Rowv = F,
          Colv = F,
          lwid=c(1,4),
          lhei = c(1,5),
          breaks = seq(-5,5,by=.5),
          dendrogram = "none")
dev.off()
```

