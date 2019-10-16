# Helicobacter-human dual-species transcriptomics in vitro and in vivo

## Table of Contents
1. [In vivo human RNA-Seq analysis](#invivohuman)
    1. [Quantify human genes from in vivo RNA-Seq samples](#invivohuman_quant)
    2. [Create human counts and TPM dataframes](#invivohuman_countstpm)
    3. [Choose top candidate genes for IPA](#invivohuman_filtergenes)

## In vivo human RNA-Seq analysis <a name="invivohuman"></a>
### Quantify human genes from in vivo RNA-Seq samples (2019/10/14) <a name="invivohuman_quant"></a>
```{bash}
for FASTQ in $(find "$FASTQ_DIR"/*R/ILLUMINA_DATA/ -name *"fastq.gz" | sed "s/_R[12].fastq.gz//g" | sort -n | uniq)
do

  SAMPLE="$(echo "$(dirname "$FASTQ")" | sed "s/\\/ILLUMINA_DATA.*//g" | sed "s/.*\\///g")"
  mkdir "$OUTPUT_DIR"/"$SAMPLE"
  
  echo -e "qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N kallisto -wd "$OUTPUT_DIR"/"$SAMPLE" -b y /usr/local/packages/kallisto-0.45.0/kallisto quant -t "$THREADS" -i "$NUC_TRANSCRIPT_FNA".kallisto.index -o "$OUTPUT_DIR"/"$SAMPLE"/"$(basename "$FASTQ")" "$FASTQ"_R1.fastq.gz "$FASTQ"_R2.fastq.gz"
  
done
```

```{bash}
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
            paste0(output.dir,"/human_invivo_counts.tsv"),,
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
write.table(tpm,
            paste0(output.dir,"/human_invivo_tpm.tsv"),,
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

