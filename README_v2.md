# Helicobacter-human dual-species transcriptomics in vitro and in vivo

# Table of Contents
<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
	- [Software](#software)
	- [Directories](#directories)
	- [Create directories](#create-directories)
	- [Set up reference files](#set-up-reference-files)
		- [Download H. pylori and human reference files](#download-h-pylori-and-human-reference-files)
		- [Create combined H. pylori and human reference](#create-combined-h-pylori-and-human-reference)
		- [Create mapping file from H. pylori GFF](#create-mapping-file-from-h-pylori-gff)
- [Confirm cag KO in H. pylori](#confirm-cag-ko-in-h-pylori)
	- [Create list of SRAs to download](#create-list-of-sras-to-download)
	- [Download H. pylori genomic FASTQs from the SRA](#download-h-pylori-genomic-fastqs-from-the-sra)
	- [Align H. pylori genomic FASTQs to H. pylori reference](#align-h-pylori-genomic-fastqs-to-h-pylori-reference)
	- [Assess depth over H. pylori genome](#assess-depth-over-h-pylori-genome)
	- [Visualize depth over H. pylori genome](#visualize-depth-over-h-pylori-genome)
		- [Set R inputs](#set-r-inputs)
		- [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo)
		- [Construct data frame of the depth at H. pylori genome from positions 500,000-600,000](#construct-data-frame-of-the-depth-at-h-pylori-genome-from-positions-500000-600000)
		- [Calculate average depth across the H. pylori knockout region at positions 547,165-584,551](#calculate-average-depth-across-the-h-pylori-knockout-region-at-positions-547165-584551)
		- [Create plot data frame](#create-plot-data-frame)
		- [Plot depth in H. pylori genome from positions 500,000-600,000](#plot-depth-in-h-pylori-genome-from-positions-500000-600000)
- [Identify SNP differences between the sequenced WT and cagKO H. pylori strains](#identify-snp-differences-between-the-sequenced-wt-and-cagko-h-pylori-strains)
	- [Generate MPILEUP files for each of the strains](#generate-mpileup-files-for-each-of-the-strains)
	- [Identify SNPs between the genomic reference files](#identify-snps-between-the-genomic-reference-files)
	- [Filter VCF for SNP/INDEL positions that are WT in two samples and homozygous in two samples](#filter-vcf-for-snpindel-positions-that-are-wt-in-two-samples-and-homozygous-in-two-samples)
	- [Identify whether SNP/INDELs occur in genic positions](#identify-whether-snpindels-occur-in-genic-positions)
		- [Prepare VCF for ANNOVAR](#prepare-vcf-for-annovar)
			- [Convert filtered VCF to ANNOVARs input format](#convert-filtered-vcf-to-annovars-input-format)
			- [Fix contig names in AVINPUT file to match contig names in GTF](#fix-contig-names-in-avinput-file-to-match-contig-names-in-gtf)
		- [Prepare transcript FASTA for ANNOVAR](#prepare-transcript-fasta-for-annovar)
			- [Format GTF for ANNOVAR by fixing all identical transcript_id values](#format-gtf-for-annovar-by-fixing-all-identical-transcript_id-values)
			- [Convert GTF to genePred format](#convert-gtf-to-genepred-format)
			- [Format genomic FASTA for ANNOVAR](#format-genomic-fasta-for-annovar)
			- [Construct transcript FASTA](#construct-transcript-fasta)
		- [Create an ANNOVAR database](#create-an-annovar-database)
		- [Find functional impact of each SNP/INDEL](#find-functional-impact-of-each-snpindel)
	- [Create a summary table for all SNPs/INDELs with functional consequences](#create-a-summary-table-for-all-snpsindels-with-functional-consequences)
		- [Set R inputs](#set-r-inputs-1)
		- [View sessionInfo](#view-sessioninfo)
		- [Change start positions in VCF to match the ANNOVAR output format](#change-start-positions-in-vcf-to-match-the-annovar-output-format)
		- [Parse VCF to show percentages of alternative basecalls for each SNPs/INDELs in each sample](#parse-vcf-to-show-percentages-of-alternative-basecalls-for-each-snpsindels-in-each-sample)
		- [Count the number of functional consequences caused by SNPs/INDELs](#count-the-number-of-functional-consequences-caused-by-snpsindels)
		- [Create a summary table for all SNPs/INDELs that cause functional changes](#create-a-summary-table-for-all-snpsindels-that-cause-functional-changes)
- [Conduct in vitro H. pylori and human RNA-Seq analysis](#conduct-in-vitro-h-pylori-and-human-rna-seq-analysis)
	- [Create SRA ID and sample map file](#create-sra-id-and-sample-map-file)
	- [Download FASTQs from SRA](#download-fastqs-from-sra)
	- [Assign InterPro descriptions and GO terms to H. pylori genes](#assign-interpro-descriptions-and-go-terms-to-h-pylori-genes)
		- [Create CDS FASTA for each H. pylori gene](#create-cds-fasta-for-each-h-pylori-gene)
		- [Assign InterPro descriptions and GO terms using InterProScan](#assign-interpro-descriptions-and-go-terms-using-interproscan)
		- [Convert the InterProScan output to a geneinfo mapping file](#convert-the-interproscan-output-to-a-geneinfo-mapping-file)
	- [Quantify transcipts](#quantify-transcipts)
		- [H. pylori](#h-pylori)
			- [Align FASTQs with splicing to combined human and H. pylori refedsInputs](#align-fastqs-with-splicing-to-combined-human-and-h-pylori-refedsinputs)
			- [Count number of reads mapping to human and H. pylori in each BAM](#count-number-of-reads-mapping-to-human-and-h-pylori-in-each-bam)
			- [Find reads that mapped to H. pylori](#find-reads-that-mapped-to-h-pylori)
			- [Subset FASTQs to only contain H. pylori-mapping reads](#subset-fastqs-to-only-contain-h-pylori-mapping-reads)
			- [Align subset FASTQs with no splicing to combined human and H. pylori reference](#align-subset-fastqs-with-no-splicing-to-combined-human-and-h-pylori-reference)
			- [Sort BAM files](#sort-bam-files)
			- [Index BAM files](#index-bam-files)
			- [Quantify H. pylori genes from BAM files](#quantify-h-pylori-genes-from-bam-files)
		- [Human](#human)
			- [Quantify human genes directly from FASTQ files](#quantify-human-genes-directly-from-fastq-files)
	- [Create group files](#create-group-files)
		- [H. pylori](#h-pylori-1)
		- [Human](#human-1)
	- [Identify differentially expressed genes](#identify-differentially-expressed-genes)
		- [H. pylori](#h-pylori-2)
			- [Set R inputs](#set-r-inputs-2)
			- [Load R functions](#load-r-functions)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo)
			- [Create counts data frame for H. pylori genes](#create-counts-data-frame-for-h-pylori-genes)
			- [Exclude genes with 0 unique positions](#exclude-genes-with-0-unique-positions)
			- [Create TPM data frame for H. pylori genes](#create-tpm-data-frame-for-h-pylori-genes)
			- [Set group levels](#set-group-levels)
			- [Conduct saturation analysis for H. pylori genes](#conduct-saturation-analysis-for-h-pylori-genes)
			- [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally)
			- [Identify differentially expressed genes with pairwise comparisons](#identify-differentially-expressed-genes-with-pairwise-comparisons)
			- [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff)
			- [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules)
			- [Plot a heatmap of the module containing the cag pathogenicity genes](#plot-a-heatmap-of-the-module-containing-the-cag-pathogenicity-genes)
			- [Plot a heatmap of the module containing the genes that differ because of SNPs/INDELs](#plot-a-heatmap-of-the-module-containing-the-genes-that-differ-because-of-snpsindels)
		- [Human](#human-2)
			- [Set R inputs](#set-r-inputs-3)
			- [Load R functions](#load-r-functions-1)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-1)
			- [Create counts data frame for human genes](#create-counts-data-frame-for-human-genes)
			- [Create TPM data frame for human genes](#create-tpm-data-frame-for-human-genes)
			- [Set group levels](#set-group-levels-1)
			- [Conduct saturation analysis for human genes](#conduct-saturation-analysis-for-human-genes)
			- [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally-1)
			- [Identify differentially expressed genes with pairwise comparisons](#identify-differentially-expressed-genes-with-pairwise-comparisons-1)
			- [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff-1)
			- [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules-1)
			- [Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs](#create-a-list-of-genes-and-their-logfc-and-fdr-values-for-each-wgcna-module-as-ipa-inputs)
	- [Combine cyan F + darkred T and cyan T + darkred F module gene lists for WGCNA analyses due to their similar expression patterns](#combine-cyan-f--darkred-t-and-cyan-t--darkred-f-module-gene-lists-for-wgcna-analyses-due-to-their-similar-expression-patterns)
	- [Assess whether excluding 24 h samples for humans allows other clusters to be resolved](#assess-whether-excluding-24-h-samples-for-humans-allows-other-clusters-to-be-resolved)
		- [Set R inputs](#set-r-inputs-4)
		- [Load R functions](#load-r-functions-2)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-2)
		- [Create counts data frame for human genes](#create-counts-data-frame-for-human-genes-1)
		- [Remove 24 h samples from analysis](#remove-24-h-samples-from-analysis)
		- [Create TPM data frame for human genes](#create-tpm-data-frame-for-human-genes-1)
		- [Set group levels](#set-group-levels-2)
		- [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally-2)
		- [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules-2)
			- [Find soft power value for WGCNA](#find-soft-power-value-for-wgcna)
			- [Identify expression modules](#identify-expression-modules)
			- [Merge similar expression modules](#merge-similar-expression-modules)
			- [Plot WGCNA expression modules as a heatmap](#plot-wgcna-expression-modules-as-a-heatmap)
		- [Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs](#create-a-list-of-genes-and-their-logfc-and-fdr-values-for-each-wgcna-module-as-ipa-inputs-1)
- [Conduct in vivo human RNA-Seq analysis](#conduct-in-vivo-human-rna-seq-analysis)
	- [Quantify human genes directly from FASTQ files](#quantify-human-genes-directly-from-fastq-files-1)
	- [Find genes in tumor and metaplasia samples that are most up- and down-regulated relative to the tumor adjacent sample](#find-genes-in-tumor-and-metaplasia-samples-that-are-most-up--and-down-regulated-relative-to-the-tumor-adjacent-sample)
		- [Set R inputs](#set-r-inputs-5)
		- [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo-1)
		- [Create counts data frame](#create-counts-data-frame)
		- [Check how many read counts were quantified from each organism](#check-how-many-read-counts-were-quantified-from-each-organism)
			- [H. pylori](#h-pylori-3)
			- [Human](#human-3)
		- [Create TPM data frame](#create-tpm-data-frame)
		- [Calculate the log2TPM ratios between the two tumor samples and the metaplasia sample versus the tumor adjacent sample](#calculate-the-log2tpm-ratios-between-the-two-tumor-samples-and-the-metaplasia-sample-versus-the-tumor-adjacent-sample)
		- [For each comparison, find the top 3000 up- and down-regulated genes relative to the tumor adjacent sample](#for-each-comparison-find-the-top-3000-up--and-down-regulated-genes-relative-to-the-tumor-adjacent-sample)
	- [Assess which in vivo human samples are most similar to the in vitro human samples](#assess-which-in-vivo-human-samples-are-most-similar-to-the-in-vitro-human-samples)
		- [Set R inputs](#set-r-inputs-6)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-3)
		- [Combine counts data frame for in vivo and in vitro human data sets](#combine-counts-data-frame-for-in-vivo-and-in-vitro-human-data-sets)
		- [Create TPM data frame for in vivo and in vitro human data sets](#create-tpm-data-frame-for-in-vivo-and-in-vitro-human-data-sets)
		- [Conduct a hierarchical cluster analysis on the TPM values across all in vivo and in vitro human samples](#conduct-a-hierarchical-cluster-analysis-on-the-tpm-values-across-all-in-vivo-and-in-vitro-human-samples)
		- [Conduct a PCA on the TPM values across all in vivo and in vitro human samples](#conduct-a-pca-on-the-tpm-values-across-all-in-vivo-and-in-vitro-human-samples)

<!-- /MarkdownTOC -->

# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software
```{bash, eval = F}
JULIA_DEPOT_PATH=/home/mattchung/.julia
PYTHON_LIB_PATH=/usr/local/packages/python-3.5/lib

JAVA_BIN_DIR=/usr/bin
JULIA_BIN_DIR=/usr/local/bin
R_BIN_DIR=/usr/local/packages/r-3.6.0/bin

ANNOVAR_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/annovar
BOWTIE2_BIN_DIR=/usr/local/packages/bowtie2-2.3.4.3
FADU_BIN_DIR=/home/mattchung/scripts/FADU
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
INTERPROSCAN_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/interproscan-5.34-73.0
KALLISTO_BIN_DIR=/usr/local/packages/kallisto-0.45.0
KENTUTILS_BIN_DIR=/usr/local/packages/kentutils/bin
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
VARSCAN_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/varscan
```

## Directories
```{bash, eval = F}
SCRIPTS_DIR=/home/mattchung/scripts/
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/EHPYL
```

## Create directories
```{bash, eval = F}
mkdir -p "$WORKING_DIR"/annovar
mkdir -p "$WORKING_DIR"/bam
mkdir -p "$WORKING_DIR"/fadu
mkdir -p "$WORKING_DIR"/kallisto
mkdir -p "$WORKING_DIR"/plots
mkdir -p "$WORKING_DIR"/reads
mkdir -p "$WORKING_DIR"/references
mkdir -p "$WORKING_DIR"/varscan
```

## Set up reference files

### Download H. pylori and human reference files

```{bash, eval = F}
wget -O "$WORKING_DIR"/references/hpylori26995.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/525/GCF_000008525.1_ASM852v1/GCF_000008525.1_ASM852v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/hpylori26995.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/525/GCF_000008525.1_ASM852v1/GCF_000008525.1_ASM852v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/hpylori26995.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/525/GCF_000008525.1_ASM852v1/GCF_000008525.1_ASM852v1_genomic.gtf.gz
gunzip "$WORKING_DIR"/references/hpylori26995.fna.gz
gunzip "$WORKING_DIR"/references/hpylori26995.gff.gz
gunzip "$WORKING_DIR"/references/hpylori26995.gtf.gz

wget -O "$WORKING_DIR"/references/hsapiensGRCh38.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz
gunzip "$WORKING_DIR"/references/hsapiensGRCh38.fna.gz
```

### Create combined H. pylori and human reference

```{bash, eval = F}
cat "$WORKING_DIR"/references/hpylori26995.fna "$WORKING_DIR"/references/hsapiensGRCh38.fna > "$WORKING_DIR"/references/combined_hsapiensGRCh38_hpylori26695.fna
```

### Create mapping file from H. pylori GFF

##### Inputs
```{bash, eval = F}
GFF3="$WORKING_DIR"/references/hpylori26995.gff
```

##### Commands
```{bash, eval = F}
"$R_BIN_DIR"/Rscript "$SCRIPTS_DIR"/gff3_to_map.R "$GFF3"
```

# Confirm cag KO in H. pylori

## Create list of SRAs to download

##### Commands
```{bash, eval = F}
vim "$WORKING_DIR"/genomic_srr.list
```

```{bash, eval = F}
SRR5410345
SRR5410346
SRR7191641
SRR7191642
```

## Download H. pylori genomic FASTQs from the SRA

##### Inputs
```{bash, eval = F}
SRR_ID_LIST="$WORKING_DIR"/genomic_srr.list
```

##### Commands
```{bash, eval = F}
while read SRR_ID
do
  qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$WORKING_DIR"/reads -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" --gzip -O "$WORKING_DIR"/reads
done < "$SRR_ID_LIST"
```

## Align H. pylori genomic FASTQs to H. pylori reference 

##### Inputs
```{bash, eval = F}
SRR_ID_LIST="$WORKING_DIR"/genomic_srr.list
THREADS=16
REF_FNA="$WORKING_DIR"/references/hpylori26995.fna
OUTPUT_DIR="$WORKING_DIR"/bam
```

##### Commands
```{bash, eval = F}
"$BOWTIE2_BIN_DIR"/bowtie2-build --large-index "$REF_FNA" "$REF_FNA"
while read SRR
do
	echo ""$BOWTIE2_BIN_DIR"/bowtie2 --threads "$THREADS" -x "$REF_FNA" -1 "$WORKING_DIR"/reads/"$SRR"_1.fastq.gz -2 "$WORKING_DIR"/reads/"$SRR"_2.fastq.gz | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=5G -N bowtie2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

## Assess depth over H. pylori genome

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
```

##### Commands
```{bash, eval = F}
for BAM in $(find "$BAM_DIR" -name "*[.]bam")
do
  "$SAMTOOLS_BIN_DIR"/samtools depth -@ "$THREADS" -aa -d 1000000 "$BAM" > "$BAM".depth
done
```

## Visualize depth over H. pylori genome

### Set R inputs
```{R}
BAM.DIR="Z:/EBMAL/mchung_dir/EHPYL/bam"
OUTPUT.DIR="Z:/EBMAL/mchung_dir/EHPYL/plots"
```

### Load R packages and view sessionInfo
```{R}
library(ggplot2)
library(reshape)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C				   LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.2.0 reshape_0.8.8

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2	 rstudioapi_0.10  magrittr_1.5     tidyselect_0.2.5 munsell_0.5.0    colorspace_1.4-1 R6_2.4.0	   rlang_0.4.0     
 [9] stringr_1.4.0    plyr_1.8.4	 dplyr_0.8.3	tools_3.5.0	grid_3.5.0	 gtable_0.3.0     withr_2.1.2	yaml_2.2.0	
[17] lazyeval_0.2.2   assertthat_0.2.1 tibble_2.1.3     crayon_1.3.4     purrr_0.3.2	reshape2_1.4.3   glue_1.3.1	 labeling_0.3    
[25] stringi_1.4.3    compiler_3.5.0   pillar_1.4.2     scales_1.0.0     pkgconfig_2.0.2
```

### Construct data frame of the depth at H. pylori genome from positions 500,000-600,000
```{R}
kovalidation_depth.files <- list.files(BAM,DIR, full.names = T,pattern="*.depth")
kovalidation_depth.df <- as.data.frame(matrix(nrow = 100001,
							    ncol = length(kovalidation_depth.files)))
rownames(kovalidation_depth.df) <- 500000:600000
colnames(kovalidation_depth.df) <- gsub("[.].*","",basename(kovalidation_depth.files))
for(i in 1:length(kovalidation_depth.files)){
  kovalidation_depth.df[,i] <- read.delim(kovalidation_depth.files[i],header = F)[500000:600000,3]
}
```

### Calculate average depth across the H. pylori knockout region at positions 547,165-584,551
```{R}
colMeans(kovalidation_depth.df[rownames(kovalidation_depth.df) %in% 547165:584551,])
```

```{R, eval = F}
SRR5410345 SRR5410346 SRR7191641 SRR7191642 
  135.0958   200.4621     0.0000     0.000
```

### Create plot data frame
```{R}
kovalidation_depth.df$position <- rownames(kovalidation_depth.df)
kovalidation_depth.melt.df <- melt(kovalidation_depth.df)
```

### Plot depth in H. pylori genome from positions 500,000-600,000 
```{R}
kovalidation.plot <- ggplot(kovalidation_depth.melt.df,aes(x=as.numeric(as.character(position)),
									     y=value,
									     group=variable,
									     color = variable,
									     fill = variable))+
  # geom_line()+
  geom_area()+
  facet_grid(rows = vars(variable))+
  guides(group = F, color = F, fill = F)+
  labs(x="position (bp)", y="sequencing depth")+
  scale_x_continuous(expand=c(0,0))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 4))

pdf(paste0(OUTPUT.DIR,"/SupFig1.pdf"),
    height=7,
    width=7)
print(kovalidation.plot)
dev.off()

png(paste0(OUTPUT.DIR,"/SupFig1.png"),
    height=7,
    width=7,
	units = "in",res=300)
print(kovalidation.plot)
dev.off()

print(kovalidation.plot)
```

![image](/images/SupFig1.png)

# Identify SNP differences between the sequenced WT and cagKO H. pylori strains

## Generate MPILEUP files for each of the strains

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
REF_FNA="$WORKING_DIR"/references/hpylori26995.fna
```

##### Commands
```{bash, eval = F}
BAM_FILES="$(find "$BAM_DIR" -name "*sortedbyposition.bam" | sort | tr '\n' ' ')"

"$SAMTOOLS_BIN_DIR"/samtools mpileup -d 1000000 -f "$REF_FNA" -o "$WORKING_DIR"/varscan/genomic.mpileup $BAM_FILES
```

## Identify SNPs between the genomic reference files

##### Inputs
```{bash, eval = F}
MPILEUP="$WORKING_DIR"/varscan/genomic.mpileup
SRR_ID_LIST="$WORKING_DIR"/genomic_srr.list
```

##### Commands
```{bash, eval = F}
"$JAVA_BIN_DIR"/java -jar "$VARSCAN_BIN_DIR"/VarScan.v2.4.4.jar mpileup2cns "$MPILEUP" --p-value 0.01 --output-vcf 1 --vcf-sample-list "$SRR_ID_LIST" > "$(echo "$MPILEUP" | sed "s/[.]mpileup$/.vcf/g")" 
```

```{bash, eval = F}
Got the following sample list:
SRR5410345	SRR5410346	SRR7191641	SRR7191642
Min coverage:   8
Min reads2:     2
Min var freq:   0.2
Min avg qual:   15
P-value thresh: 0.01
Reading input from /local/projects-t3/EBMAL/mchung_dir/EHPYL/varscan/genomic.mpileup
1667375 bases in pileup file
256 variant positions (144 SNP, 114 indel)
44 were failed by the strand-filter
256 variant positions reported (144 SNP, 114 indel)
```
## Filter VCF for SNP/INDEL positions that are WT in two samples and homozygous in two samples

##### Inputs
```{bash, eval = F}
VCF="$WORKING_DIR"/varscan/genomic.vcf
```

##### Commands
```{bash, eval = F}
awk '$5 != "." && $7 == "PASS" {print $0}' $VCF | grep WT=2 | grep HOM=2 > "$(echo "$VCF" | sed "s/[.]vcf$/.filtered.vcf/g")"
```

## Identify whether SNP/INDELs occur in genic positions

### Prepare VCF for ANNOVAR

#### Convert filtered VCF to ANNOVARs input format

##### Inputs
```{bash, eval = F}
FILTERED_VCF="$WORKING_DIR"/varscan/genomic.filtered.vcf
```

##### Commands
```{bash, eval = F}
"$ANNOVAR_BIN_DIR"/convert2annovar.pl -format vcf4 "$FILTERED_VCF" > "$(echo "$VCF" | sed "s/[.]vcf$/.filtered.avinput/g")"
```

```{bash, eval = F}
NOTICE: Read 40 lines and wrote 37 different variants at 40 genomic positions (33 SNPs and 7 indels)
NOTICE: Among 40 different variants at 40 positions, 0 are heterozygotes, 37 are homozygotes
NOTICE: Among 33 SNPs, 19 are transitions, 14 are transversions (ratio=1.36)
```

#### Fix contig names in AVINPUT file to match contig names in GTF

##### Inputs
```{bash, eval = F}
AVINPUT="$WORKING_DIR"/varscan/genomic.filtered.avinput
```

##### Commands
```{bash, eval = F}
sed -i "s/gi|15644634|ref|NC_000915.1|/NC_000915.1/g" "$AVINPUT"
```

### Prepare transcript FASTA for ANNOVAR

#### Format GTF for ANNOVAR by fixing all identical transcript_id values

Most transcript ID values in the H. pylori gtf were "unknown_transcript_1"

##### Inputs
```{bash, eval = F}
GTF="$WORKING_DIR"/references/hpylori26995.gtf
```

##### Commands
```{bash, eval = F}
rm "$(echo "$GTF" | sed "s/[.]gtf$/.annovar.gtf/g")"

IFS=""

cat "$GTF" | while read LINE; do
	if [[ $LINE == *"unknown_transcript_1"* ]]; then
		GENE_ID="$(echo "$LINE" | cut -f9 | sed "s/.;.*//g" | sed "s/gene_id .//g")"
		echo "$(echo "$LINE" | sed "s/unknown_transcript_1/"$GENE_ID"/g")" >> "$(echo "$GTF" | sed "s/[.]gtf$/.annovar.gtf/g")"
	else
		echo "$LINE" >> "$(echo "$GTF" | sed "s/[.]gtf$/.annovar.gtf/g")"
	fi
done
```

#### Convert GTF to genePred format

##### Inputs
```{bash, eval = F}
ANNOVAR_GTF="$WORKING_DIR"/references/hpylori26995.annovar.gtf
```

##### Commands
```{bash, eval = F}
"$KENTUTILS_BIN_DIR"/gtfToGenePred -genePredExt "$ANNOVAR_GTF" "$(echo "$ANNOVAR_GTF" | sed "s/[.]gtf$/.refGene.txt/g")"
```

#### Format genomic FASTA for ANNOVAR

##### Inputs
```{bash, eval = F}
REF_FNA="$WORKING_DIR"/references/hpylori26995.fna
```

##### Commands
```{bash, eval = F}
sed  "s/>.*/>NC_000915.1/g" "$REF_FNA" > "$(echo "$REF_FNA" | sed "s/[.]fna/.annovar.fna/g")"
```

#### Construct transcript FASTA 

##### Inputs
```{bash, eval = F}
ANNOVAR_REF_FNA="$WORKING_DIR"/references/hpylori26995.annovar.fna
REFGENE="$WORKING_DIR"/references/hpylori26995.annovar.refGene.txt
```

##### Commands
```{bash, eval = F}
"$ANNOVAR_BIN_DIR"/retrieve_seq_from_fasta.pl --format refGene --seqfile "$ANNOVAR_REF_FNA" "$REFGENE" --out "$(echo "$REFGENE" | sed "s/[.]refGene.txt$/.mRNA.fna/g")"
```

### Create an ANNOVAR database

##### Inputs
```{bash, eval = F}
DB_NAME=hpylori26995
REFGENE="$WORKING_DIR"/references/hpylori26995.annovar.refGene.txt
TRANSCRIPT_FNA="$WORKING_DIR"/references/hpylori26995.annovar.mRNA.fna
```

##### Commands
```{bash, eval = F}
mkdir "$ANNOVAR_BIN_DIR"/"$DB_NAME"/
cp "$REFGENE" "$ANNOVAR_BIN_DIR"/"$DB_NAME"/"$DB_NAME"_refGene.txt
cp "$TRANSCRIPT_FNA" "$ANNOVAR_BIN_DIR"/"$DB_NAME"/"$DB_NAME"_refGeneMrna.fa
```

### Find functional impact of each SNP/INDEL

##### Inputs
```{bash, eval = F}
DB_NAME=hpylori26995
AV_INPUT="$WORKING_DIR"/varscan/genomic.filtered.avinput
```

##### Commands
```{bash, eval = F}
"$ANNOVAR_BIN_DIR"/annotate_variation.pl --geneanno --buildver "$DB_NAME" "$AV_INPUT" "$ANNOVAR_BIN_DIR"/"$DB_NAME" --outfile "$WORKING_DIR"/annovar/annovar
```

## Create a summary table for all SNPs/INDELs with functional consequences


### Set R inputs
```{R}
ANNOVAR_OUTPUT.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/annovar/annovar.exonic_variant_function"
GFF_MAP.PATH <-"Z:/EBMAL/mchung_dir/EHPYL/references/hpylori26995.gff.map"
VCF.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/varscan/genomic.vcf"
```

### View sessionInfo
```{R}
sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C				   LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_3.5.0 tools_3.5.0    yaml_2.2.
```

### Change start positions in VCF to match the ANNOVAR output format

```{R}
vcf <- read.delim(VCF.PATH, comment.char = "#", header = F)
vcf <- vcf[grep("WT=2;HET=0;HOM=2",vcf[,8]),]

vcf[,2] <- apply(vcf, 1, function(x){
  if(nchar(x[4]) > nchar(x[5])){
    return(as.numeric(as.character(x[2])) + 1)
  }else{
    return(as.numeric(as.character(x[2])))
  }
})
```

### Parse VCF to show percentages of alternative basecalls for each SNPs/INDELs in each sample
```{R}
for(i in 10:13){
  vcf[,i] <- unlist(lapply(vcf[,i], function(x){unlist(strsplit(as.character(x), split = ":"))[7]}))
}
```

### Count the number of functional consequences caused by SNPs/INDELs
```{R}
annovar_output <- read.delim(ANNOVAR_OUTPUT.PATH, header = F)
table(annovar_output[,2])
```

```{R, eval = F}
 frameshift deletion frameshift insertion    nonsynonymous SNV		 stopgain	 synonymous SNV		  unknown 
			 3			  1			  0			  1			  2			  0
```

### Create a summary table for all SNPs/INDELs that cause functional changes

```{R}
annovar_output <- annovar_output[annovar_output[,2] != "nonsynonymous SNV" & annovar_output[,2] != "unknown",]

gff_map <- read.delim(GFF_MAP.PATH)
gff_map <- gff_map[match(gsub(":.*","",annovar_output[,3]),gsub("gene-","",gsub("[|].*","",gff_map$ID))),]

vcf <- vcf[match(annovar_output[,5],vcf[,2]),]

annovar_functionalinfo <- as.data.frame(cbind(gsub(":.*","",annovar_output[,3]),
							    annovar_output[,2],
							    annovar_output[,5:8],
							    gsub("%2C.*","",gff_map$old_locus_tag),
							    gff_map$product,
							    vcf[,10:13]))
rownames(annovar_functionalinfo) <- NULL
colnames(annovar_functionalinfo) <- c("gene_id","mutation_type","start","stop","ref","alt","old_locus_tag","product",
						  "SRR5410345","SRR5410346","SRR7191641","SRR7191642")

print(annovar_functionalinfo)
```

```{R, eval = F}
     gene_id	  mutation_type   start    stop ref alt old_locus_tag						 product SRR5410345
1 HP_RS02780		 stopgain  597715  597715   G   A	  HP0565				 YkgB family protein	 100%
2 HP_RS03230	 synonymous SNV  703308  703308   C   A	  HP0655 outer membrane protein assembly factor BamA	   0%
3 HP_RS03300 frameshift insertion  721378  721378   -   A	  HP0671			    outer membrane protein	   0%
4 HP_RS04245  frameshift deletion  922126  922127  TA   -	  HP0870		     flagellar hook protein FlgE     96.61%
5 HP_RS04335  frameshift deletion  943757  943757   G   -	  HP0889		   iron ABC transporter permease	   0%
6 HP_RS07010  frameshift deletion 1487681 1487682  GA   -	    <NA>		 phosphoethanolamine transferase	75.3%
7 HP_RS07365	 synonymous SNV 1559723 1559723   G   A	  HP1487			  ABC transporter permease	   0%
  SRR5410346 SRR7191641 SRR7191642
1	 100%	   0%	   0%
2	   0%	 100%	 100%
3	   0%     97.44%     97.54%
4     97.24%	   0%	   0%
5	   0%     98.48%     97.95%
6     77.23%     10.37%	4.93%
7	   0%     99.54%	 100%
```

# Conduct in vitro H. pylori and human RNA-Seq analysis

## Create SRA ID and sample map file

##### Commands
```{bash, eval = F}
vim "$WORKING_DIR"/srr_sample.map
```

```{bash, eval = F}
SRR5423603	N87_a	polyA
SRR5423604	N87_a	polyA
SRR5422947	N87_b	polyA
SRR5422948	N87_b	polyA
SRR5423599	N87_c	polyA
SRR5422749	N87_Hp_cag_2h_a	agss
SRR5422857	N87_Hp_cag_2h_a	polyA
SRR5422858	N87_Hp_cag_2h_a	polyA
SRR5422746	N87_Hp_cag_2h_b	agss
SRR5423069	N87_Hp_cag_2h_b	polyA
SRR5423070	N87_Hp_cag_2h_b	polyA
SRR5423636	N87_Hp_cag_2h_c	polyA
SRR5422747	N87_Hp_cag_4h_a	agss
SRR5422856	N87_Hp_cag_4h_a	polyA
SRR5422748	N87_Hp_cag_4h_b	agss
SRR5423058	N87_Hp_cag_4h_b	polyA
SRR5422946	N87_Hp_cag_4h_c	polyA
SRR5422750	N87_Hp_cag_24h_a	agss
SRR5422751	N87_Hp_cag_24h_a	agss
SRR5422752	N87_Hp_cag_24h_a	agss
SRR5423659	N87_Hp_cag_24h_a	polyA
SRR5422753	N87_Hp_cag_24h_b	agss
SRR5422754	N87_Hp_cag_24h_b	agss
SRR5422755	N87_Hp_cag_24h_b	agss
SRR5423632	N87_Hp_cag_24h_b	polyA
SRR5423658	N87_Hp_cag_24h_c	polyA
SRR5422756	N87_Hp_ko_2h_a	agss
SRR5423638	N87_Hp_ko_2h_a	polyA
SRR5422833	N87_Hp_ko_2h_b	agss
SRR5423059	N87_Hp_ko_2h_b	polyA
SRR5422928	N87_Hp_ko_2h_c	polyA
SRR5422831	N87_Hp_ko_4h_a	agss
SRR5423655	N87_Hp_ko_4h_a	polyA
SRR5422832	N87_Hp_ko_4h_b	agss
SRR5422927	N87_Hp_ko_4h_b	polyA
SRR5422945	N87_Hp_ko_4h_c	polyA
SRR5422834	N87_Hp_ko_24h_b	agss
SRR5422835	N87_Hp_ko_24h_b	agss
SRR5422836	N87_Hp_ko_24h_b	agss
SRR5424793	N87_Hp_ko_24h_b	polyA
SRR5422837	N87_Hp_ko_24h_c	agss
SRR5422838	N87_Hp_ko_24h_c	agss
SRR5422839	N87_Hp_ko_24h_c	agss
SRR5424794	N87_Hp_ko_24h_c	polyA
SRR5422043	Hp_cag_0h_a	ribozero
SRR5422045	Hp_cag_0h_b	ribozero
SRR5422044	Hp_cag_0h_c	ribozero
SRR5422048	Hp_cag_2h_a	ribozero
SRR5422046	Hp_cag_2h_b	ribozero
SRR5422047	Hp_cag_2h_c	ribozero
SRR5422049	Hp_cag_4h_a	ribozero
SRR5422050	Hp_cag_4h_b	ribozero
SRR5422741	Hp_cag_4h_c	ribozero
SRR5422742	Hp_cag_24h_a	ribozero
SRR5422743	Hp_cag_24h_b	ribozero
SRR5422744	Hp_cag_24h_c	ribozero
SRR5422840	Hp_ko_0h_a	ribozero
SRR5422841	Hp_ko_0h_b	ribozero
SRR5422842	Hp_ko_0h_c	ribozero
SRR5422843	Hp_ko_0h_c	ribozero
SRR5422844	Hp_ko_0h_c	ribozero
SRR5422845	Hp_ko_2h_a	ribozero
SRR5422846	Hp_ko_2h_b	ribozero
SRR5422847	Hp_ko_2h_c	ribozero
SRR5422848	Hp_ko_4h_a	ribozero
SRR5422849	Hp_ko_4h_b	ribozero
SRR5422850	Hp_ko_4h_c	ribozero
SRR5422851	Hp_ko_24h_a	ribozero
SRR5422852	Hp_ko_24h_b	ribozero
SRR5422853	Hp_ko_24h_c	ribozero
SRR5422854	Hp_ko_24h_c	ribozero
SRR5422855	Hp_ko_24h_c	ribozero
```

## Download FASTQs from SRA

##### Inputs
```{bash, eval = F}
READS_DIR="$WORKING_DIR"/reads
SRR2SAMPLE_MAP="$WORKING_DIR"/srr_sample.map
```

##### Commands
```{bash, eval = F}
cut -f1 "$SRR2SAMPLE_MAP" | while read SRR_ID
do
	qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --gzip --split-files "$SRR_ID" -O "$READS_DIR"
done
```

## Assign InterPro descriptions and GO terms to H. pylori genes

### Create CDS FASTA for each H. pylori gene

##### Inputs
```{bash, eval = F}
FNA="$WORKING_DIR"/references/hpylori26995.fna
GFF3="$WORKING_DIR"/references/hpylori26995
FEAT_TYPE=CDS
ID_ATTR=ID
```

##### Commands
```{bash, eval = F}
"$SCRIPTS_DIR"/createnuccdsfasta.sh -n "$FNA" -g "$GFF3" -f "$FEAT_TYPE" -i "$ID_ATTR" > "$(echo "$FNA" | sed "s/[.]fna/cds.fna/")"
```

### Assign InterPro descriptions and GO terms using InterProScan

##### Inputs
```{bash, eval = F}
CDS_FNA="$WORKING_DIR"/references/hpylori26995.cds.fna
THREADS=4
```

##### Commands
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_PATH":"$LD_LIBRARY_PATH"\n"$INTERPROSCAN_BIN_DIR"/interproscan.sh -i "$CDS_FNA" -f tsv -o "$CDS_FNA".interproscan.tsv --seqtype n --goterms --iprlookup" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=20G -N interproscan -wd "$(dirname "$CDS_FNA")"
```

### Convert the InterProScan output to a geneinfo mapping file

##### Inputs
```{bash, eval = F}
IPRSCAN_OUTPUT="$WORKING_DIR"/references/hpylori26995.cds.fna.interproscan.tsv
```

##### Commands
```{bash, eval = F}
"$R_BIN_DIR"/Rscript "$SCRIPTS_DIR"/interproscan2geneinfo.R "$IPRSCAN_OUTPUT"
```

## Quantify transcipts

### H. pylori

#### Align FASTQs with splicing to combined human and H. pylori refedsInputs
```{bash, eval = F}
SRR_ID_LIST="$WORKING_DIR"/PRJNA378649.srr.list
REF_FNA="$WORKING_DIR"/references/combined_hsapiensGRCh38_hpylori26695.fna
OUTPUT_DIR="$WORKING_DIR"/bam
READS_DIR="$WORKING_DIR"/reads
THREADS=4
```

##### Commands
```{bash, eval = F}
"$HISAT2_BIN_DIR"/hisat2-build "$REF_FNA" "$REF_FNA"
while read SRR
do
	echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" -1 "$READS_DIR"/"$SRR"_1.fastq.gz -2 "$READS_DIR"/"$SRR"_2.fastq.gz | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
```

#### Count number of reads mapping to human and H. pylori in each BAM

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
```

##### Commnads
```{bash, eval = F}
rm "$WORKING_DIR"/mapped_stats.tsv
for BAM in $(find "$BAM_DIR" -name "*[.]bam")
do

	echo -e ""$(basename "$BAM" | sed "s/[.]bam$//g")"\t"$("$SAMTOOLS_BIN_DIR"/samtools view -F 4 -F 256 "$BAM" | grep NC_000915.1 | wc -l)"\t"$("$SAMTOOLS_BIN_DIR"/samtools view -F 4 -F 256 "$BAM" | grep -v NC_000915.1 | wc -l)"" >> "$WORKING_DIR"/mapped_stats.tsv &
done
```

#### Find reads that mapped to H. pylori

##### Inputs
```{bash, eval = F}
QUERY=NC_000915.1
BAM_DIR="$WORKING_DIR"/bam
```

##### Commands
```{bash, eval = F}
for BAM in $(find "$BAM_DIR" -name "*[.]bam")
do
  "$SAMTOOLS_BIN_DIR"/samtools view "$BAM" | grep "$QUERY" | awk '{print $1}' | sort -n | uniq > "$BAM"."$QUERY".reads
done
```

#### Subset FASTQs to only contain H. pylori-mapping reads

##### Inputs
```{bash, eval = F}
SRR_ID_LIST="$WORKING_DIR"/PRJNA378649.srr.list
READS_DIR="$WORKING_DIR"/reads
BAM_DIR="$WORKING_DIR"/bam
```

##### Commands
```{bash, eval = F}
while read SRR
do
	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$READS_DIR"/"$SRR"_1.fastq "$BAM_DIR"/"$SRR".bam."$QUERY".reads > "$READS_DIR"/"$SRR"_1.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$READS_DIR"
	echo -e ""$SEQTK_BIN_DIR"/seqtk subseq "$READS_DIR"/"$SRR"_2.fastq "$BAM_DIR"/"$SRR".bam."$QUERY".reads > "$READS_DIR"/"$SRR"_2.subset.fastq" | qsub -P jdhotopp-lab -l mem_free=2G -N seqtk -wd "$READS_DIR"
done < "$SRR_ID_LIST"
```

#### Align subset FASTQs with no splicing to combined human and H. pylori reference

##### Inputs
```{bash, eval = F}
SRR_ID_LIST="$WORKING_DIR"/PRJNA378649.srr.list
REF_FNA="$WORKING_DIR"/references/combined_hsapiensGRCh38_hpylori26695.fna
BAM_DIR="$WORKING_DIR"/bam
READS_DIR="$WORKING_DIR"/reads
THREADS=4
```

##### Commands
```{bash, eval = F}
while read SRR
do
	echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" -1 "$READS_DIR"/"$SRR"_1.subset.fastq -2 "$READS_DIR"/"$SRR"_2.subset.fastq --no-spliced-alignment -X 1000 | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$BAM_DIR"/"$SRR".subset.bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$BAM_DIR"
done < "$SRR_ID_LIST"
```

#### Sort BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
THREADS=4
```

##### Commands
```{bash, eval = F}
for BAM in $(find $BAM_DIR -name "*[.]subset.bam")
do
	qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=2G -N sort -wd "$BAM_DIR" -b y "$SAMTOOLS_BIN_DIR"/samtools sort "$BAM" -@ "$THREADS" -o "$(echo $BAM | sed "s/[.]bam$/.sortedbyposition.bam"/g)"
done
```

#### Index BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
THREADS=4
```

##### Commands
```{bash, eval = F}
for BAM in $(find $BAM_DIR -name "*[.]subset.sortedbyposition.bam")
do
	qsub -P jdhotopp-lab -l mem_free=2G -N index -wd "$BAM_DIR" -b y $SAMTOOLS_BIN_DIR"/samtools index "$BAM"
done
```

#### Quantify H. pylori genes from BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
GFF="$WORKING_DIR"/references/hpylori26995.gff
FEAT_TYPE="CDS"
STRANDEDNESS="no"
ATTR_ID="ID"
OUTPUT_DIR="$WORKING_DIR"/fadu
```

##### Commands
```{bash, eval = F}
for BAM in $(find "$BAM_DIR" -name "*[.]sortedbyposition.bam")
do
	echo -e "export JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -g "$GFF" -b "$BAM" -o "$OUTPUT_DIR" -s "$STRANDEDNESS" -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -l mem_free=5G -N fadu -wd "$OUTPUT_DIR"
done
```

### Human

#### Quantify human genes directly from FASTQ files

##### Inputs
```{bash, eval = F}
SRR_ID_LIST="$WORKING_DIR"/PRJNA378649.srr.list
REF_TRANSCRIPT_FNA="$WORKING_DIR"/references/combined_hsapiensGRCh38_hpylori26695.cds.fna
READS_DIR="$WORKING_DIR"/reads
THREADS=4
```

##### Commands
```{bash, eval = F}
"$KALLISTO_BIN_DIR"/kallisto index -i "$REF_TRANSCRIPT_FNA".kallisto.index --make-unique "$REF_TRANSCRIPT_FNA"
while read SRR
do
	mkdir "$OUTPUT_DIR"/"$SRR"
	qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N kallisto -wd "$OUTPUT_DIR" -b y "$KALLISTO_BIN_DIR"/kallisto quant -t "$THREADS" -i "$REF_TRANSCRIPT_FNA".kallisto.index -o "$OUTPUT_DIR"/"$SRR" "$READS_DIR"/"$SRR"_1.fastq "$READS_DIR"/"$SRR"_2.fastq
done < "$SRR_ID_LIST"
```

## Create group files
### H. pylori
```{R}
vim "$WORKING_DIR"/hpylori_groups.tsv
```

```{R}
N87_Hp_cag_2h_a	N87_Hp_cag_2h	#92dae1	16	N87 + HP WT 2 h
N87_Hp_cag_2h_b	N87_Hp_cag_2h	#92dae1	16	N87 + HP WT 2 h
N87_Hp_cag_4h_a	N87_Hp_cag_4h	#5bc8d2	17	N87 + HP WT 4 h
N87_Hp_cag_4h_b	N87_Hp_cag_4h	#5bc8d2	17	N87 + HP WT 4 h
N87_Hp_cag_24h_a	N87_Hp_cag_24h	#24b6c3	15	N87 + HP WT 24 h
N87_Hp_cag_24h_b	N87_Hp_cag_24h	#24b6c3	15	N87 + HP WT 24 h
N87_Hp_ko_2h_a	N87_Hp_ko_2h	#e28e88	16	N87 + HP cag- 2 h
N87_Hp_ko_2h_b	N87_Hp_ko_2h	#e28e88	16	N87 + HP cag- 2 h
N87_Hp_ko_4h_a	N87_Hp_ko_4h	#d4564d	17	N87 + HP cag- 4 h
N87_Hp_ko_4h_b	N87_Hp_ko_4h	#d4564d	17	N87 + HP cag- 4 h
N87_Hp_ko_24h_b	N87_Hp_ko_24h	#C51D12	15	N87 + HP cag- 24 h
N87_Hp_ko_24h_c	N87_Hp_ko_24h	#C51D12	15	N87 + HP cag- 24 h
Hp_cag_0h_a	Hp_cag_0h	#BB99FF	18	HP WT 0 h
Hp_cag_0h_b	Hp_cag_0h	#BB99FF	18	HP WT 0 h
Hp_cag_0h_c	Hp_cag_0h	#BB99FF	18	HP WT 0 h
Hp_cag_2h_a	Hp_cag_2h	#9966FF	16	HP WT 2 h
Hp_cag_2h_b	Hp_cag_2h	#9966FF	16	HP WT 2 h
Hp_cag_2h_c	Hp_cag_2h	#9966FF	16	HP WT 2 h
Hp_cag_4h_a	Hp_cag_4h	#7733FF	17	HP WT 4 h
Hp_cag_4h_b	Hp_cag_4h	#7733FF	17	HP WT 4 h
Hp_cag_4h_c	Hp_cag_4h	#7733FF	17	HP WT 4 h
Hp_cag_24h_a	Hp_cag_24h	#5500FF	15	HP WT 24 h
Hp_cag_24h_b	Hp_cag_24h	#5500FF	15	HP WT 24 h
Hp_cag_24h_c	Hp_cag_24h	#5500FF	15	HP WT 24 h
Hp_ko_0h_a	Hp_ko_0h	#FFCC99	18	HP cag- 0 h
Hp_ko_0h_b	Hp_ko_0h	#FFCC99	18	HP cag- 0 h
Hp_ko_0h_c	Hp_ko_0h	#FFCC99	18	HP cag- 0 h
Hp_ko_2h_a	Hp_ko_2h	#FFB266	16	HP cag- 2 h
Hp_ko_2h_b	Hp_ko_2h	#FFB266	16	HP cag- 2 h
Hp_ko_2h_c	Hp_ko_2h	#FFB266	16	HP cag- 2 h
Hp_ko_4h_a	Hp_ko_4h	#FF9933	17	HP cag- 4 h
Hp_ko_4h_b	Hp_ko_4h	#FF9933	17	HP cag- 4 h
Hp_ko_4h_c	Hp_ko_4h	#FF9933	17	HP cag- 4 h
Hp_ko_24h_a	Hp_ko_24h	#FF8000	15	HP cag- 24 h
Hp_ko_24h_b	Hp_ko_24h	#FF8000	15	HP cag- 24 h
Hp_ko_24h_c	Hp_ko_24h	#FF8000	15	HP cag- 24 h
```

### Human
```{R}
vim "$WORKING_DIR"/human_groups.tsv
```

```{R}
N87_a	N87	#c1e31d	18	N87
N87_b	N87	#c1e31d	18	N87
N87_c	N87	#c1e31d	18	N87
N87_Hp_cag_2h_a	N87_Hp_cag_2h	#92dae1	16	N87 + HP WT 2 h
N87_Hp_cag_2h_b	N87_Hp_cag_2h	#92dae1	16	N87 + HP WT 2 h	
N87_Hp_cag_2h_c	N87_Hp_cag_2h	#92dae1	16	N87 + HP WT 2 h
N87_Hp_cag_4h_a	N87_Hp_cag_4h	#5bc8d2	17	N87 + HP WT 4 h
N87_Hp_cag_4h_b	N87_Hp_cag_4h	#5bc8d2	17	N87 + HP WT 4 h
N87_Hp_cag_4h_c	N87_Hp_cag_4h	#5bc8d2	17	N87 + HP WT 4 h
N87_Hp_cag_24h_a	N87_Hp_cag_24h	#24b6c3	15	N87 + HP WT 24 h
N87_Hp_cag_24h_b	N87_Hp_cag_24h	#24b6c3	15	N87 + HP WT 24 h
N87_Hp_cag_24h_c	N87_Hp_cag_24h	#24b6c3	15	N87 + HP WT 24 h
N87_Hp_ko_2h_a	N87_Hp_ko_2h	#e28e88	16	N87 + HP cag- 2 h
N87_Hp_ko_2h_b	N87_Hp_ko_2h	#e28e88	16	N87 + HP cag- 2 h
N87_Hp_ko_2h_c	N87_Hp_ko_2h	#e28e88	16	N87 + HP cag- 2 h
N87_Hp_ko_4h_a	N87_Hp_ko_4h	#d4564d	17	N87 + HP cag- 4 h
N87_Hp_ko_4h_b	N87_Hp_ko_4h	#d4564d	17	N87 + HP cag- 4 h
N87_Hp_ko_4h_c	N87_Hp_ko_4h	#d4564d	17	N87 + HP cag- 4 h
N87_Hp_ko_24h_b	N87_Hp_ko_24h	#C51D12	15	N87 + HP cag- 24 h
N87_Hp_ko_24h_c	N87_Hp_ko_24h	#C51D12	15	N87 + HP cag- 24 h
```

## Identify differentially expressed genes

### H. pylori

#### Set R inputs
```{R}
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/hpylori_groups.tsv"
SRR2SAMPLE_MAP.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/srr_sample.map"
GFF3_MAP.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/references/hpylori26995.gff.map"
GENEINFO.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/references/hpylori26995.cds.fna.interproscan.geneinfo.tsv"
WORKING.DIR="Z:/EBMAL/mchung_dir/EHPYL/"
```

#### Load R functions

```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}

get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}

get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
							    ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}

find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(dy) < 0.05)[1]
  if(softpower == 1){
    softpower <- which(abs(dy) < 0.05)[2]
  }
  return(softpower)
}

eigengene_invert_id <- function(tpm.de, mergedColors, mergedMEs){
  tpm.de.wgcna <- tpm.de
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de[i,]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
	tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}

wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
								 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
								  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
								  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}

get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}

functionaltermenrichment <- function(genes, geneinfo){
  for(i in 1:ncol(geneinfo)){geneinfo[,i] <- as.character(geneinfo[,i])}
  geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
  geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
  geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
  geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo$interpro_description, collapse = "|"),  split = "[|]")))),
					 gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
					 gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
					 gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  geneinfo.subset <- geneinfo[geneinfo$gene %in% genes,]
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()

  functionalterms.list.subset <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$interpro_description, collapse = "|"),  split = "[|]")))),
						  gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
						  gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
						  gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(i in 1:length(functionalterms.list)){
    for(j in 1:nrow(functionalterms.list[[i]])){
	freq.all <- functionalterms.list[[i]][j,2]
	freq.subset <- ifelse(functionalterms.list[[i]][j,1] %in% functionalterms.list.subset[[i]][,1],
				    functionalterms.list.subset[[i]][functionalterms.list.subset[[i]][,1] == as.character(functionalterms.list[[i]][j,1]),2],
				    0)
	genes.all <- nrow(geneinfo)
	genes.subset <- nrow(geneinfo.subset)

	fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
						 genes.subset - freq.subset, genes.all - genes.subset - freq.all + freq.subset),
					     nrow = 2,
					     ncol = 2)
	fisher.test <- fisher.test(fisherexact.matrix)
	
	term[length(term) + 1] <- as.character(functionalterms.list[[i]][j,1])
	clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
	genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
	pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
	correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "fdr", n = nrow(functionalterms.list[[i]]))
	oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
					    clusteroccurences,
					    genomeoccurences,
					    pvalue,
					    correctedpvalue,
					    oddsratio))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}
```

#### Load packages and view sessionInfo
```{R}
library(dendextend)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(pvclust)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C				  
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vegan_2.5-5		     lattice_0.20-35		 permute_0.9-5		  
 [4] pvclust_2.0-0		   gridExtra_2.3		   gplots_3.0.1.1		 
 [7] ggplot2_3.2.0		   ggdendro_0.1-20		 FactoMineR_1.42		
[10] edgeR_3.24.3		    limma_3.38.3		    DESeq2_1.22.2		  
[13] SummarizedExperiment_1.12.0 DelayedArray_0.8.0	    BiocParallel_1.16.6	  
[16] matrixStats_0.54.0	    Biobase_2.42.0		  GenomicRanges_1.34.0	 
[19] GenomeInfoDb_1.18.2	   IRanges_2.16.0		  S4Vectors_0.20.1	     
[22] BiocGenerics_0.28.0	   dendextend_1.12.0	    

loaded via a namespace (and not attached):
 [1] nlme_3.1-137	     bitops_1.0-6	     bit64_0.9-7		RColorBrewer_1.1-2    
 [5] tools_3.5.0		backports_1.1.4	  R6_2.4.0		   KernSmooth_2.23-15    
 [9] rpart_4.1-13	     mgcv_1.8-23		Hmisc_4.2-0		DBI_1.0.0		 
[13] lazyeval_0.2.2	   colorspace_1.4-1	 nnet_7.3-12		withr_2.1.2	     
[17] tidyselect_0.2.5	 bit_1.1-14		 compiler_3.5.0	   htmlTable_1.13.1	
[21] flashClust_1.01-2	caTools_1.17.1.2	 scales_1.0.0	     checkmate_1.9.4	 
[25] genefilter_1.64.0	stringr_1.4.0	    digest_0.6.20	    foreign_0.8-70	  
[29] XVector_0.22.0	   base64enc_0.1-3	  pkgconfig_2.0.2	  htmltools_0.3.6	 
[33] htmlwidgets_1.3	  rlang_0.4.0		rstudioapi_0.10	  RSQLite_2.1.2	   
[37] gtools_3.8.1	     acepack_1.4.1	    dplyr_0.8.3		RCurl_1.95-4.12	 
[41] magrittr_1.5	     GenomeInfoDbData_1.2.0 Formula_1.2-3	    leaps_3.0		 
[45] Matrix_1.2-14	    Rcpp_1.0.2		 munsell_0.5.0	    viridis_0.5.1	   
[49] scatterplot3d_0.3-41   stringi_1.4.3	    yaml_2.2.0		 MASS_7.3-51.4	   
[53] zlibbioc_1.28.0	  grid_3.5.0		 blob_1.2.0		 gdata_2.18.0	    
[57] crayon_1.3.4	     splines_3.5.0	    annotate_1.60.1	  locfit_1.5-9.1	  
[61] zeallot_0.1.0	    knitr_1.23		 pillar_1.4.2	     geneplotter_1.60.0    
[65] XML_3.98-1.20	    glue_1.3.1		 latticeExtra_0.6-28    data.table_1.12.2     
[69] vctrs_0.2.0		gtable_0.3.0	     purrr_0.3.2		assertthat_0.2.1	
[73] xfun_0.8		   xtable_1.8-4	     survival_2.41-3	  viridisLite_0.3.0     
[77] tibble_2.1.3	     AnnotationDbi_1.44.0   memoise_1.1.0	    cluster_2.0.7-1   
```

#### Create counts data frame for H. pylori genes

```{R}
groups <- read.delim(GROUPS.PATH, header = F)
srr2sample_map <- read.delim(SRR2SAMPLE_MAP.PATH, header = F)
fadu_output <- list.files(paste0(WORKING.DIR,"/fadu/"), recursive = T, full.names = T, pattern = "counts.txt")

srr2sample_map <- srr2sample_map[srr2sample_map[,3] != "polyA",]

counts <- as.data.frame(matrix(0,
					  nrow = nrow(read.delim(fadu_output[1],header = T)),
					  ncol = length(unique(srr2sample_map[,2]))))
rownames(counts) <- read.delim(fadu_output[2],header = T)[,1]
colnames(counts) <- unique(srr2sample_map[,2])

for(i in 1:ncol(counts)){
  srr2sample_map.subset <- srr2sample_map[srr2sample_map[,2] == colnames(counts)[i],]
  for(j in 1:nrow(srr2sample_map.subset)){
    counts.subset <- read.delim(grep(srr2sample_map.subset[j,1],fadu_output,value=T),header = T)
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),4]
  }
}
counts <- counts[grep("^HP",rownames(counts)),match(groups[,1],colnames(counts))]

write.table(counts,
		paste0(WORKING.DIR,"/hpylori_counts.tsv"),
		quote = F,
		col.names = T,
		row.names = T,
		sep = "\t")

sort(colSums(counts))
```

```{R, eval = F}
    Hp_cag_24h_c      Hp_cag_4h_a     Hp_cag_24h_b       Hp_ko_2h_c      Hp_ko_24h_b     Hp_cag_24h_a      Hp_ko_24h_a       Hp_ko_0h_b 
         5139048          6026150          7117276          7516840          7732775          8322814          8382444          9189260 
     Hp_cag_2h_b       Hp_ko_4h_b       Hp_ko_4h_a      Hp_cag_4h_c       Hp_ko_0h_a       Hp_ko_4h_c      Hp_cag_0h_c      Hp_cag_0h_b 
         9398915          9554572          9637180         11127451         11153378         11689215         11795834         12161553 
      Hp_ko_2h_a      Hp_cag_2h_c      Hp_ko_24h_c      Hp_cag_0h_a      Hp_cag_2h_a       Hp_ko_2h_b       Hp_ko_0h_c      Hp_cag_4h_b 
        12337557         12589539         12856730         13300949         13487127         13681967         14458855         15245200 
  N87_Hp_ko_4h_b N87_Hp_cag_24h_b   N87_Hp_ko_4h_a  N87_Hp_cag_4h_b  N87_Hp_cag_2h_b  N87_Hp_cag_2h_a  N87_Hp_cag_4h_a N87_Hp_cag_24h_a 
        16423134         18108534         19372683         22006944         26935535         30723695         31325897         34305328 
 N87_Hp_ko_24h_b   N87_Hp_ko_2h_a   N87_Hp_ko_2h_b  N87_Hp_ko_24h_c 
        37932073         43678354         44564830         59250547 
```

#### Exclude genes with 0 unique positions

```{R}
genelength <- read.delim(fadu_output[1],header = T)[match(rownames(counts),read.delim(fadu_output[1],header = T)[,1]),2]

counts <- counts[genelength != 0,]
genelength <- genelength[genelength != 0]
```

#### Create TPM data frame for H. pylori genes
```{R}
tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(tpm,
		paste0(WORKING.DIR,"/hpylori_tpm.tsv"),
		quote = F,
		col.names = T,
		row.names = T,
		sep = "\t")

dim(tpm)
```

```{R, eval = F}
[1] 1445   36
```

#### Set group levels

```{R}
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
groups[,5] <- factor(groups[,5], levels = unique(groups[,5]))
```

#### Conduct saturation analysis for H. pylori genes

```{R, fig.height=5, fig.width=6}
rarefy.counts <- round(counts,0)
raremax <- round(min(rowSums(t(rarefy.counts))),0)
srare <- rarefy(t(rarefy.counts),raremax) 

rarefy.raw.df <- rarecurve(t(rarefy.counts), step = round(raremax/10,0), sample = raremax)

rarefy.df <- as.data.frame(matrix(nrow = 0,
					    ncol = 5))
rarefy.points.df <- rarefy.df
for(i in 1:length(rarefy.raw.df)){
  steps <- as.numeric(gsub("N","",names(rarefy.raw.df[[i]])))
  detected_genes <- as.numeric(rarefy.raw.df[[i]])
  rarefy.df <- as.data.frame(rbind(rarefy.df,
					     cbind(as.numeric(steps),as.numeric(detected_genes),as.character(groups[i,1]),as.character(groups[i,2]),groups[i,3])))
  rarefy.points.df <- as.data.frame(rbind(rarefy.points.df,
							cbind(as.numeric(max(steps)),as.numeric(max(detected_genes)),as.character(groups[i,1]),as.character(groups[i,2],groups[i,3]))))
  
}
rarefy.plot <- ggplot()+
  geom_line(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  #geom_point(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  geom_point(mapping=aes(x=as.numeric(as.character(rarefy.points.df[,1])), y=as.numeric(as.character(rarefy.points.df[,2])),group=rarefy.points.df[,3],color=rarefy.points.df[,4]),size = 3)+
  guides(colour = F,shape = F)+
  scale_color_manual(values = levels(groups[,3]))+
  labs(x="reads mapping to genes", y="genes detected", color = "Sample")+
  #coord_cartesian(xlim=c(0,100000))+
  theme_bw()


pdf(paste0(WORKING.DIR,"/plots/hpylori_rarefication_plot.pdf"),
    height=5,
    width=6)
print(rarefy.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_rarefication_plot.png"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

![image](/images/hpylori_rarefication_plot.png)

#### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

##### DESeq2

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

deseq.groups <- groups[,1:2]
deseq.groups <- deseq.groups[,2,drop = F]

rownames(deseq.groups) <- groups[,1]
colnames(deseq.groups) <- "condition"

dds <- DESeqDataSetFromMatrix(countData = round(counts),
					colData = deseq.groups,
					design = ~ condition)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
print(paste0(length(which(keep == F))," genes excluded with CPM cutoff"))

dds <- DESeq(dds, test="LRT", reduced=~1)
#res <- results(dds)
res <- results(dds, cooksCutoff=T)
res <- res[!is.na(res$padj),]
print(paste0(nrow(res[res$padj < FDRcutoff,])," DE genes identified using DESeq2 longitudinal"))

write.table(res[res$padj < FDRcutoff,],
		paste0(WORKING.DIR,"/hpylori_counts_deseq2_longitudinal.tsv"),
		row.names = T,
		col.names = NA,
		quote = F,
		sep = "\t")
```

```{R, eval = F}
[1] "0 genes excluded with CPM cutoff"
[1] "1444 DE genes identified using DESeq2 longitudinal"
```

##### edgeR

```{R}
y <- DGEList(counts = counts, group = groups[,2])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,2]))
keep.df <- as.data.frame(table(keep))
print(paste0(length(keep.df[keep.df[,1] == F,2])," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,2])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))
qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
edgeR.longitudinal.degenes <- qlf$table[qlf$table$padj < FDRcutoff,]
print(paste0(nrow(edgeR.longitudinal.degenes)," DE genes identified using edgeR longitudinal"))

write.table(edgeR.longitudinal.degenes,
		paste0(WORKING.DIR,"/human_counts_edgeR_longitudinal.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")

counts.keep <- counts[keep,]
tpm.keep <- tpm[keep,]
```

```{R, eval = F}
[1] "0 genes excluded with CPM cutoff"
[1] "1445 DE genes identified using edgeR longitudinal
```

#### Identify differentially expressed genes with pairwise comparisons

##### DESeq2

```{R}
edgeR.pairwise.degenes <- list(N87_Hp_cag_2h_vs_N87_Hp_ko_2h=c(), 
					 N87_Hp_cag_4h_vs_N87_Hp_ko_4h=c(),
					 N87_Hp_cag_24h_vs_N87_Hp_ko_24h=c(),
					 Hp_cag_0h_vs_Hp_ko_0h=c(), 
					 Hp_cag_2h_vs_Hp_ko_2h=c(), 
					 Hp_cag_4h_vs_Hp_ko_4h=c(),
					 Hp_cag_24h_vs_Hp_ko_24h=c())

deseq2.pairwise.degenes <- edgeR.pairwise.degenes 

for(i in 1:length(deseq2.pairwise.degenes)){
  print(names(deseq2.pairwise.degenes)[i])
  group1 <- substr(names(deseq2.pairwise.degenes)[i],1,regexpr("_vs_",names(deseq2.pairwise.degenes)[i])-1)
  group2 <- substr(names(deseq2.pairwise.degenes)[i],regexpr("_vs_",names(deseq2.pairwise.degenes)[i])+4,nchar(names(deseq2.pairwise.degenes)[i]))
  
  deseq.pairwise.groups <- deseq.groups[c(grep(paste0("^",group1),deseq.groups[,1]),grep(paste0("^",group2),deseq.groups[,1])),,drop = F]
  counts.pairwise <- counts[,colnames(counts) %in% rownames(deseq.pairwise.groups)]

  dds <- DESeqDataSetFromMatrix(countData = round(counts.pairwise),
					colData = deseq.pairwise.groups,
					design = ~ condition)
  dds <- DESeq(dds, test="LRT", reduced=~1)
  keep <- rowSums(counts(dds)) >= 5
  dds <- dds[keep,]
  print(paste0(length(which(keep == F))," genes excluded with CPM cutoff"))

  #res <- results(dds)
  res <- results(dds, cooksCutoff=T)
  res <- res[!is.na(res$padj),]
  deseq.pairwise.degenes <- res[res$padj < FDRcutoff,]
  print(paste0(nrow(deseq.pairwise.degenes)," DE genes identified using DESeq2"))

  write.table(res[res$padj < FDRcutoff,],
		  paste0(WORKING.DIR,"/hpylori_counts_deseq2_",names(deseq.pairwise.degenes)[i],".tsv"),
		  row.names = T,
		  col.names = NA,
		  quote = F,
		  sep = "\t")
}
```

```{R, eval = F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "0 genes excluded with CPM cutoff"
[1] "820 DE genes identified using DESeq2"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "2 genes excluded with CPM cutoff"
[1] "1102 DE genes identified using DESeq2"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "0 genes excluded with CPM cutoff"
[1] "995 DE genes identified using DESeq2"

[1] "Hp_cag_0h_vs_Hp_ko_0h"
[1] "0 genes excluded with CPM cutoff"
[1] "935 DE genes identified using DESeq2"

[1] "Hp_cag_2h_vs_Hp_ko_2h"
[1] "0 genes excluded with CPM cutoff"
[1] "996 DE genes identified using DESeq2"

[1] "Hp_cag_4h_vs_Hp_ko_4h"
[1] "0 genes excluded with CPM cutoff"
[1] "881 DE genes identified using DESeq2"

[1] "Hp_cag_24h_vs_Hp_ko_24h"
[1] "0 genes excluded with CPM cutoff"
[1] "1057 DE genes identified using DESeq2"
```

##### edgeR

```{R}
for(i in 1:length(edgeR.pairwise.degenes)){
  print(names(edgeR.pairwise.degenes)[i])
  group1 <- substr(names(edgeR.pairwise.degenes)[i],1,regexpr("_vs_",names(edgeR.pairwise.degenes)[i])-1)
  group2 <- substr(names(edgeR.pairwise.degenes)[i],regexpr("_vs_",names(edgeR.pairwise.degenes)[i])+4,nchar(names(edgeR.pairwise.degenes)[i]))
  
  groups.pairwise <- groups[c(grep(paste0("^",group1),groups[,2]),grep(paste0("^",group2),groups[,2])),]
  groups.pairwise[,2] <- factor(groups.pairwise[,2],levels=unique(groups.pairwise[,2]))
  
  counts.pairwise <- counts[,colnames(counts) %in% groups.pairwise[,1]]
  
  cpm.cutoff <- 5/min(colSums(counts.pairwise)) * 1000000
  
  y <- DGEList(counts = counts.pairwise, group = groups.pairwise[,2])
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups.pairwise[,2]))
  keep.df <- as.data.frame(table(keep))
  print(paste0(length(keep.df[keep.df[,1] == F,2])," genes excluded with CPM cutoff"))
  
  y <- y[keep, , keep.lib.sizes = F]
  design <- model.matrix(~groups.pairwise[,2])
  y <- estimateDisp(y , design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
  edgeR.pairwise.degenes[[i]] <- rownames(qlf$table)[qlf$table$padj < FDRcutoff]
  print(paste0(length(edgeR.pairwise.degenes[[i]])," DE genes identified using edgeR"))


  write.table(qlf$table[qlf$table$padj < FDRcutoff,],
		  paste0(WORKING.DIR,"/hpylori_counts_edgeR_",names(edgeR.pairwise.degenes)[i],".tsv"),
		  row.names = T,
		  col.names = NA,
		  quote = F,
		  sep = "\t")
}
```

```{R, eval = F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "1 genes excluded with CPM cutoff"
[1] "831 DE genes identified using edgeR"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "1 genes excluded with CPM cutoff"
[1] "1137 DE genes identified using edgeR"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "1 genes excluded with CPM cutoff"
[1] "976 DE genes identified using edgeR"

[1] "Hp_cag_0h_vs_Hp_ko_0h"
[1] "0 genes excluded with CPM cutoff"
[1] "943 DE genes identified using edgeR"

[1] "Hp_cag_2h_vs_Hp_ko_2h"
[1] "0 genes excluded with CPM cutoff"
[1] "1009 DE genes identified using edgeR"

[1] "Hp_cag_4h_vs_Hp_ko_4h"
[1] "1 genes excluded with CPM cutoff"
[1] "838 DE genes identified using edgeR"

[1] "Hp_cag_24h_vs_Hp_ko_24h"
[1] "1 genes excluded with CPM cutoff"
[1] "1062 DE genes identified using edgeR"
```
#### Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff

##### Create sample legend for PCA and hierarchical clustering plots

```{R, fig.height = 2, fig.width = 10}
legend.plot <- ggplot(mapping=aes(x=groups[,5], y=seq(1,length(groups[,2]),1), group = groups[,1]))+
    geom_point(aes(color = groups[,5],shape=groups[,5]), size = 4)+
    scale_shape_manual(values = as.numeric(as.character(groups[!duplicated(groups[,3]),4])))+
    scale_color_manual(values = as.character(unique(groups[,3])))+
    guides(shape = guide_legend(title = "Samples", title.position = "top",nrow=2),
	     colour = guide_legend(title = "Samples", title.position = "top",nrow=2))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/hpylori_pca_hc_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_pca_hc_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/hpylori_pca_hc_samplekey.png)

##### Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff

```{R, fig.height=5, fig.width=10}
dendrogram <- as.data.frame(t(scale(t(log2(tpm.keep+1)))))

result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100)

structure <- get_dendro_structure(result)
dendro.data <- get_dendro_data(result)
bootstrap.positions <- get_dendro_bootstraps(dendro.data)
  
points.df <- as.data.frame(cbind(seq(1,length(structure),1),
					   structure))
dendrogroups <- groups[,2][result$hclust$order]
dendrocol <- groups[,3][result$hclust$order]
dendroshape <- groups[,4][result$hclust$order]
# dendrosize <- colSums(counts.keep)[result$hclust$order]
dendrosize <- 1

dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
  geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize, shape = dendroshape))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  scale_color_manual(values = levels(groups[,3]))+
  labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
  guides(colour = guide_legend(ncol = 2), size = F, shape = F)+
  #scale_x_discrete(limits = as.character(dendrolabels))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
	  axis.text.y = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank())

for(i in 1:length(result$edges$bp)){
  text <- round(result$edges$bp[i] * 100,0)
  dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[i,1] + 0.4, y=bootstrap.positions[i,2] + 0.04, size = 2)
}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_kept_dendrogram.pdf"),
    height=5,
    width=10)
print(dendrogram.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_kept_dendrogram.png"),
    height=5,
    width=10,
	units = "in",res=300)
print(dendrogram.plot)
dev.off()

print(dendrogram.plot)
```

![image](/images/hpylori_tpm_kept_dendrogram.png)
##### Conduct a PCA on the TPM values of all genes that passed CPM cutoff

```{R,fig.height=5,fig.width=5}
pca.df <- t(scale(t(log2(tpm.keep + 1))))
pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.keep) - 1)

pca.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = groups[,2],size = 1, shape = dendroshape))+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
	 x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
	 y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  guides(color = F,size = F, shape = F)+
  # guides(colour = guide_legend(ncol = 2))+
  scale_color_manual(values = levels(groups[,3]))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_kept_pca.pdf"),
    height=5,
    width=5)
print(pca.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_kept_pca.png"),
    height=5,
    width=5,
	units = "in",res=300)
print(pca.plot)
dev.off()

print(pca.plot)
```

![image](/images/hpylori_tpm_kept_pca.png)


#### Divide differentially expressed genes into expression modules

##### Find soft power value for WGCNA

```{R}
tpm.de <- tpm[rownames(tpm) %in% rownames(edgeR.longitudinal.degenes),]

wgcna <- as.data.frame(t(tpm.de))
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(wgcna, powerVector = powers, verbose = 5)
softpower <- find_soft_power(sft)

text.color <- rep("black",length(sft$fitIndices[,1]))
text.color[which(sft$fitIndices[,1] == softpower)] <- "red"

scale_independence.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Scale Independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
  theme_bw()

mean_connectivity.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Mean Connectivity", x = "soft threshold (power)", y = "mean connectivity")+
  theme_bw()

wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
lay <- rbind(c(1,2))

pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_soft_power_plot.pdf"),
    width = 10, 
    height = 5)
grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_soft_power_plot.png"),
    width = 10, 
    height = 5,
	units = "in",res=300)
grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)
dev.off()

grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)

```

![image](/images/hpylori_tpm_de_wgcna_soft_power_plot.png)

##### Identify expression modules 

```{R, fig.height=5, fig.width=12}
adjacency <- adjacency(wgcna, power = softpower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");

minModuleSize <- 1
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
				     deepSplit = 2, pamRespectsDendro = FALSE,
				     minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(wgcna, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
METree = hclust(as.dist(MEDiss), method = "average")

MEDissThres = 0.25

pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_merge_eigengenes_plot.pdf"),
    width = 12, 
    height = 5)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_merge_eigengenes_plot.png"),
    width = 12, 
    height = 5,
    units = "in",res=300)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![image](/images/hpylori_tpm_de_wgcna_merge_eigengenes_plot.png)

##### Merge similar expression modules

```{R}
merge = mergeCloseModules(wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

tpm.de.wgcna <- eigengene_invert_id(tpm.de, mergedColors, mergedMEs)

write.table(tpm.de.wgcna,
		paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_modules.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")

pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_eigengene_dendrogram.pdf"),
    width = 8, 
    height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_eigengene_dendrogram.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

![image](/images/hpylori_tpm_de_wgcna_eigengene_dendrogram.png)

##### Plot WGCNA expression modules as a heatmap

```{R}
tpm.de.wgcna <- wgcna_heatmap_reorder(tpm.de.wgcna)

log2tpm.de <- log2(tpm.de.wgcna[,1:(ncol(tpm.de.wgcna) - 2)] + 1)
zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
rowcol1 <- tpm.de.wgcna$module
rowcol2 <- unlist(lapply(tpm.de.wgcna$invert,function(x){if(x == F){return("grey")}else{return("black")}}))
colcol <- as.character(groups[,3])

rowsep <- get_heatmap_separators(rowcol1)
colsep <- get_heatmap_separators(colcol)

rev(sort(table(paste0(tpm.de.wgcna$module))))
```

```{R, eval = F}
    firebrick4          brown          black  darkseagreen4    darkorange2    floralwhite     darkorange           cyan         purple 
           416            204            185            150             77             72             61             43             41 
lightsteelblue  darkslateblue          ivory  antiquewhite4       darkgrey     indianred4       thistle1         coral1     lightcoral 
            39             37             28             27             17             10              8              7              5 
         plum3     darkviolet        thistle        salmon2 palevioletred2           grey 
             4              4              3              3              3              1 
```


```{R}
rev(sort(table(paste0(tpm.de.wgcna$module,tpm.de.wgcna$invert))))
```

```{R, eval = F}
    firebrick4FALSE          brownFALSE  darkseagreen4FALSE          blackFALSE    darkorange2FALSE           blackTRUE      firebrick4TRUE 
                363                 180                 134                 130                  73                  55                  53 
   floralwhiteFALSE     darkorangeFALSE           cyanFALSE lightsteelblueFALSE  darkslateblueFALSE         purpleFALSE          ivoryFALSE 
                 50                  41                  38                  37                  31                  29                  24 
          brownTRUE  antiquewhite4FALSE     floralwhiteTRUE      darkorangeTRUE       darkgreyFALSE   darkseagreen4TRUE          purpleTRUE 
                 24                  24                  22                  20                  17                  16                  12 
    indianred4FALSE       thistle1FALSE         coral1FALSE   darkslateblueTRUE     lightcoralFALSE            cyanTRUE           ivoryTRUE 
                  9                   8                   7                   6                   5                   5                   4 
    darkorange2TRUE        thistleFALSE        salmon2FALSE          plum3FALSE palevioletred2FALSE     darkvioletFALSE   antiquewhite4TRUE 
                  4                   3                   3                   3                   3                   3                   3 
 lightsteelblueTRUE           plum3TRUE      indianred4TRUE           greyFALSE      darkvioletTRUE 
                  2                   1                   1                   1                   1 
```

###### Create sample legend for WGCNA heatmap

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,5], y=seq(1,length(groups[,2]),1), group = groups[,1]))+
    geom_line(aes(color = groups[,5]), size = 4)+
    scale_color_manual(values = as.character(unique(groups[,3])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=5))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/hpylori_hm_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_hm_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/hpylori_hm_samplekey.png)

###### Create z-score log2TPM legend for WGCNA heatmap

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
			     colours=hmcol,
			     breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/hpylori_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/hpylori_hm_zscorelog2tpmkey.png)

###### Use module assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_module_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/hpylori_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png)

###### Use inverse assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/hpylori_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png)

#### Plot a heatmap of the module containing the cag pathogenicity genes

###### Create heatmap row labels containing gene names and encoded products 

```{R}
module <- "darkslateblue"

gff3_map <- read.delim(GFF3_MAP.PATH, comment.char = "#", header = T)
gff3_map$old_locus_tag <- gsub("%2C.*","",gff3_map$old_locus_tag)

zscore.log2tpm.cag <- zscore.log2tpm.de[intersect(which(rowcol1 == module),which(rowcol2 == "grey")),]
labrow1 <- paste0(rownames(zscore.log2tpm.cag),": ", gff3_map$product[which(gff3_map$old_locus_tag %in% rownames(zscore.log2tpm.cag))])

zscore.log2tpm.cag <- zscore.log2tpm.de[intersect(which(rowcol1 == module),which(rowcol2 == "black")),]
labrow2 <- paste0(rownames(zscore.log2tpm.cag),": ", gff3_map$product[which(gff3_map$old_locus_tag %in% rownames(zscore.log2tpm.cag))])

zscore.log2tpm.cag <- zscore.log2tpm.de[rowcol1 == module,]
labrow <- c(labrow1,labrow2)
```

###### Use module assigments as the row color bar

```{R,fig.height=5.5, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_cagmodule_heatmap.pdf"),
    width = 5, 
    height = 5.5)
heatmap.2(as.matrix(zscore.log2tpm.cag),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[which(rowcol1 == module)],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.cag)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_cagmodule_heatmap.png"),
    width = 5, 
    height = 5.5,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.cag),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[which(rowcol1 == module)],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.cag)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.cag),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[which(rowcol1 == module)],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.cag)),
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/hpylori_tpm_de_wgcna_zscorelog2tpm_cagmodule_heatmap.png)

###### Use inverse assigments as the row color bar

```{R,fig.height=5.5, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_caginvert_heatmap.pdf"),
    width = 5, 
    height = 5.5)
heatmap.2(as.matrix(zscore.log2tpm.cag),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2[which(rowcol1 == module)],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.cag)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_caginvert_heatmap.png"),
    width = 5, 
    height = 5.5,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.cag),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2[which(rowcol1 == module)],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.cag)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.cag),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2[which(rowcol1 == module)],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.cag)),
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/hpylori_tpm_de_wgcna_zscorelog2tpm_caginvert_heatmap.png)

#### Plot a heatmap of the module containing the genes that differ because of SNPs/INDELs

###### Create heatmap row labels containing gene names and encoded products 

```{R}
snpindel.genes <- c("HP0565","HP0655","HP0671","HP0870","HP0889","HP1487")
zscore.log2tpm.snpindel <- zscore.log2tpm.de[rownames(zscore.log2tpm.de) %in% snpindel.genes,]

labrow <- paste0(rownames(zscore.log2tpm.snpindel),": ", gff3_map$product[gff3_map$old_locus_tag %in% rownames(zscore.log2tpm.snpindel)])
```

###### Use module assigments as the row color bar

```{R,fig.height=5.5, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_snpindelmodule_heatmap.pdf"),
    width = 5, 
    height = 5.5)
heatmap.2(as.matrix(zscore.log2tpm.snpindel),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[rownames(zscore.log2tpm.de) %in% snpindel.genes],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.snpindel)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_snpindelmodule_heatmap.png"),
    width = 5, 
    height = 5.5,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.snpindel),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[rownames(zscore.log2tpm.de) %in% snpindel.genes],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.snpindel)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.snpindel),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[rownames(zscore.log2tpm.de) %in% snpindel.genes],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.snpindel)),
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/hpylori_tpm_de_wgcna_zscorelog2tpm_snpindelmodule_heatmap.png)

###### Use inverse assigments as the row color bar

```{R,fig.height=5.5, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_snpindelinvert_heatmap.pdf"),
    width = 5, 
    height = 5.5)
heatmap.2(as.matrix(zscore.log2tpm.snpindel),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[rownames(zscore.log2tpm.de) %in% snpindel.genes],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.snpindel)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

png(paste0(WORKING.DIR,"/plots/hpylori_tpm_de_wgcna_zscorelog2tpm_snpindelinvert_heatmap.png"),
    width = 5, 
    height = 5.5,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.snpindel),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[rownames(zscore.log2tpm.de) %in% snpindel.genes],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.snpindel)),
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.snpindel),
		  col=hmcol,
		  trace="none",
		  labRow=labrow,
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1[rownames(zscore.log2tpm.de) %in% snpindel.genes],
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = seq(1:nrow(zscore.log2tpm.snpindel)),
		  colsep = colsep,
		  dendrogram = "none")
```

### Human

#### Set R inputs
```{R}
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/human_groups.tsv"
SRR2SAMPLE_MAP.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/srr_sample.map"
WORKING.DIR="Z:/EBMAL/mchung_dir/EHPYL"
```

#### Load R functions

```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}

get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}

get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
							    ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}
find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(dy) < 0.05)[1]
  if(softpower == 1){
    softpower <- which(abs(dy) < 0.05)[2]
  }
  return(softpower)
}

eigengene_invert_id <- function(tpm.de, mergedColors, mergedMEs){
  tpm.de.wgcna <- tpm.de
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de[i,]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
	tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}

wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
								 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
								  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
								  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}

get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}
```

#### Load packages and view sessionInfo
```{R}
library(dendextend)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(pvclust)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C				  
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WGCNA_1.68			fastcluster_1.1.25	    dynamicTreeCut_1.63-1	
 [4] vegan_2.5-5		     lattice_0.20-35		 permute_0.9-5		  
 [7] pvclust_2.0-0		   gridExtra_2.3		   gplots_3.0.1.1		 
[10] ggplot2_3.2.0		   ggdendro_0.1-20		 FactoMineR_1.42		
[13] edgeR_3.24.3		    limma_3.38.3		    DESeq2_1.22.2		  
[16] SummarizedExperiment_1.12.0 DelayedArray_0.8.0	    BiocParallel_1.16.6	  
[19] matrixStats_0.54.0	    Biobase_2.42.0		  GenomicRanges_1.34.0	 
[22] GenomeInfoDb_1.18.2	   IRanges_2.16.0		  S4Vectors_0.20.1	     
[25] BiocGenerics_0.28.0	   dendextend_1.12.0	    

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1	 htmlTable_1.13.1	 XVector_0.22.0	   base64enc_0.1-3	 
 [5] rstudioapi_0.10	  bit64_0.9-7		mvtnorm_1.0-11	   AnnotationDbi_1.44.0  
 [9] codetools_0.2-15	 splines_3.5.0	    leaps_3.0		  doParallel_1.0.15     
[13] impute_1.56.0	    robustbase_0.93-5	geneplotter_1.60.0     knitr_1.23		
[17] zeallot_0.1.0	    Formula_1.2-3	    annotate_1.60.1	  cluster_2.0.7-1	 
[21] GO.db_3.7.0		rrcov_1.4-7		compiler_3.5.0	   backports_1.1.4	 
[25] assertthat_0.2.1	 Matrix_1.2-14	    lazyeval_0.2.2	   acepack_1.4.1	   
[29] htmltools_0.3.6	  tools_3.5.0		gtable_0.3.0	     glue_1.3.1		
[33] GenomeInfoDbData_1.2.0 dplyr_0.8.3		Rcpp_1.0.2		 vctrs_0.2.0	     
[37] gdata_2.18.0	     preprocessCore_1.44.0  nlme_3.1-137	     iterators_1.0.12	
[41] xfun_0.8		   stringr_1.4.0	    gtools_3.8.1	     XML_3.98-1.20	   
[45] DEoptimR_1.0-8	   zlibbioc_1.28.0	  MASS_7.3-51.4	    scales_1.0.0	    
[49] RColorBrewer_1.1-2     yaml_2.2.0		 memoise_1.1.0	    rpart_4.1-13	    
[53] latticeExtra_0.6-28    stringi_1.4.3	    RSQLite_2.1.2	    genefilter_1.64.0     
[57] pcaPP_1.9-73	     foreach_1.4.7	    checkmate_1.9.4	  caTools_1.17.1.2	
[61] rlang_0.4.0		pkgconfig_2.0.2	  bitops_1.0-6	     purrr_0.3.2	     
[65] htmlwidgets_1.3	  labeling_0.3	     bit_1.1-14		 tidyselect_0.2.5	
[69] robust_0.4-18.1	  magrittr_1.5	     R6_2.4.0		   fit.models_0.5-14     
[73] Hmisc_4.2-0		DBI_1.0.0		  pillar_1.4.2	     foreign_0.8-70	  
[77] withr_2.1.2		mgcv_1.8-23		survival_2.41-3	  scatterplot3d_0.3-41  
[81] RCurl_1.95-4.12	  nnet_7.3-12		tibble_2.1.3	     crayon_1.3.4	    
[85] KernSmooth_2.23-15     viridis_0.5.1	    locfit_1.5-9.1	   grid_3.5.0		
[89] data.table_1.12.2	blob_1.2.0		 digest_0.6.20	    flashClust_1.01-2     
[93] xtable_1.8-4	     munsell_0.5.0	    viridisLite_0.3.0
```

#### Create counts data frame for human genes

```{R}
groups <- read.delim(GROUPS.PATH, header = F)
srr2sample_map <- read.delim(SRR2SAMPLE_MAP.PATH, header = F)
kallisto_output <- list.files(paste0(WORKING.DIR,"/kallisto/"), recursive = T, full.names = T, pattern = "abundance.tsv")

srr2sample_map <- srr2sample_map[srr2sample_map[,3] == "polyA",]

counts <- as.data.frame(matrix(0,
					  nrow = nrow(read.delim(kallisto_output[1],header = T)),
					  ncol = length(unique(srr2sample_map[,2]))))
rownames(counts) <- read.delim(kallisto_output[1],header = T)[,1]
colnames(counts) <- unique(srr2sample_map[,2])

for(i in 1:ncol(counts)){
  srr2sample_map.subset <- srr2sample_map[srr2sample_map[,2] == colnames(counts)[i],]
  for(j in 1:nrow(srr2sample_map.subset)){
    counts.subset <- read.delim(grep(srr2sample_map.subset[j,1],kallisto_output,value=T),header = T)
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),4]
  }
}
counts <- counts[grep("ENST",rownames(counts)),match(groups[,1],colnames(counts))]

write.table(counts,
		paste0(WORKING.DIR,"/human_counts.tsv"),
		quote = F,
		col.names = T,
		row.names = T,
		sep = "\t")
```

#### Create TPM data frame for human genes

```{R}
genelength <- read.delim(kallisto_output[1],header = T)[match(rownames(counts),read.delim(kallisto_output[1],header = T)[,1]),3]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(tpm,
		paste0(WORKING.DIR,"/human_tpm.tsv"),
		quote = F,
		col.names = T,
		row.names = T,
		sep = "\t")

dim(tpm)
```

```{R, eval = F}
[1] 188753     20
```

#### Set group levels

```{R}
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
groups[,5] <- factor(groups[,5], levels = unique(groups[,5]))
```

#### Conduct saturation analysis for human genes

```{R, fig.height=5, fig.width=6}
rarefy.counts <- round(counts,0)
raremax <- round(min(rowSums(t(rarefy.counts))),0)
srare <- rarefy(t(rarefy.counts),raremax) 

rarefy.raw.df <- rarecurve(t(rarefy.counts), step = round(raremax/10,0), sample = raremax)

rarefy.df <- as.data.frame(matrix(nrow = 0,
					    ncol = 5))
rarefy.points.df <- rarefy.df
for(i in 1:length(rarefy.raw.df)){
  steps <- as.numeric(gsub("N","",names(rarefy.raw.df[[i]])))
  detected_genes <- as.numeric(rarefy.raw.df[[i]])
  rarefy.df <- as.data.frame(rbind(rarefy.df,
					     cbind(as.numeric(steps),as.numeric(detected_genes),as.character(groups[i,1]),as.character(groups[i,2]),groups[i,3])))
  rarefy.points.df <- as.data.frame(rbind(rarefy.points.df,
							cbind(as.numeric(max(steps)),as.numeric(max(detected_genes)),as.character(groups[i,1]),as.character(groups[i,2],groups[i,3]))))
  
}
rarefy.plot <- ggplot()+
  geom_line(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  #geom_point(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  geom_point(mapping=aes(x=as.numeric(as.character(rarefy.points.df[,1])), y=as.numeric(as.character(rarefy.points.df[,2])),group=rarefy.points.df[,3],color=rarefy.points.df[,4]),size = 3)+
  guides(colour = F,shape = F)+
  scale_color_manual(values = levels(groups[,3]))+
  labs(x="reads mapping to genes", y="genes detected", color = "Sample")+
  #coord_cartesian(xlim=c(0,100000))+
  theme_bw()


pdf(paste0(WORKING.DIR,"/plots/human_rarefication_plot.pdf"),
    height=5,
    width=6)
print(rarefy.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_rarefication_plot.png"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

![image](/images/human_rarefication_plot.png)

#### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

##### DESeq2

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

deseq.groups <- groups[,1:2]
deseq.groups <- deseq.groups[,2,drop = F]

rownames(deseq.groups) <- groups[,1]
colnames(deseq.groups) <- "condition"

dds <- DESeqDataSetFromMatrix(countData = round(counts),
					colData = deseq.groups,
					design = ~ condition)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
print(paste0(length(which(keep == F))," genes excluded with CPM cutoff"))

dds <- DESeq(dds, test="LRT", reduced=~1)
#res <- results(dds)
res <- results(dds, cooksCutoff=T)
res <- res[!is.na(res$padj),]
print(paste0(nrow(res[res$padj < FDRcutoff,])," DE genes identified using DESeq2 longitudinal"))

write.table(res[res$padj < FDRcutoff,],
		paste0(WORKING.DIR,"/human_counts_deseq2_longitudinal.tsv"),
		row.names = T,
		col.names = NA,
		quote = F,
		sep = "\t")
```

```{R, eval = F}
[1] "70975 genes excluded with CPM cutoff"
[1] "13 DE genes identified using DESeq2 longitudinal"
```

##### edgeR

```{R}
y <- DGEList(counts = counts, group = groups[,2])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,2]))
keep.df <- as.data.frame(table(keep))
print(paste0(length(keep.df[keep.df[,1] == F,2])," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,2])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))
qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
edgeR.longitudinal.degenes <- qlf$table[qlf$table$padj < FDRcutoff,]
print(paste0(nrow(edgeR.longitudinal.degenes)," DE genes identified using edgeR longitudinal"))

write.table(edgeR.longitudinal.degenes,
		paste0(WORKING.DIR,"/human_counts_edgeR_longitudinal.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")

counts.keep <- counts[keep,]
tpm.keep <- tpm[keep,]
```

```{R, eval = F}
[1] "99691 genes excluded with CPM cutoff"
[1] "952 DE genes identified using edgeR longitudinal
```

#### Identify differentially expressed genes with pairwise comparisons

##### DESeq2

```{R}
edgeR.pairwise.degenes <- list(N87_Hp_cag_2h_vs_N87_Hp_ko_2h=c(), 
					 N87_Hp_cag_4h_vs_N87_Hp_ko_4h=c(),
					 N87_Hp_cag_24h_vs_N87_Hp_ko_24h=c())

deseq2.pairwise.degenes <- edgeR.pairwise.degenes 

for(i in 1:length(deseq2.pairwise.degenes)){
  print(names(deseq2.pairwise.degenes)[i])
  group1 <- substr(names(deseq2.pairwise.degenes)[i],1,regexpr("_vs_",names(deseq2.pairwise.degenes)[i])-1)
  group2 <- substr(names(deseq2.pairwise.degenes)[i],regexpr("_vs_",names(deseq2.pairwise.degenes)[i])+4,nchar(names(deseq2.pairwise.degenes)[i]))
  
  deseq.pairwise.groups <- deseq.groups[c(grep(group1,deseq.groups[,1]),grep(group2,deseq.groups[,1])),,drop = F]
  counts.pairwise <- counts[,colnames(counts) %in% rownames(deseq.pairwise.groups)]

  dds <- DESeqDataSetFromMatrix(countData = round(counts.pairwise),
					colData = deseq.pairwise.groups,
					design = ~ condition)
  dds <- DESeq(dds, test="LRT", reduced=~1)
  keep <- rowSums(counts(dds)) >= 5
  dds <- dds[keep,]
  print(paste0(length(which(keep == F))," genes excluded with CPM cutoff"))

  #res <- results(dds)
  res <- results(dds, cooksCutoff=T)
  res <- res[!is.na(res$padj),]
  deseq.pairwise.degenes <- res[res$padj < FDRcutoff,]
  print(paste0(length(deseq.pairwise.degenes[[i]])," DE genes identified using DESeq2"))

  write.table(res[res$padj < FDRcutoff,],
		  paste0(WORKING.DIR,"/human_counts_deseq2_",names(deseq.pairwise.degenes)[i],".tsv"),
		  row.names = T,
		  col.names = NA,
		  quote = F,
		  sep = "\t")
}
```

```{R, eval = F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "101766 genes excluded with CPM cutoff"
[1] "0 DE genes identified using DESeq2"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "109084 genes excluded with CPM cutoff"
[1] "0 DE genes identified using DESeq2"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "112638 genes excluded with CPM cutoff"
[1] "0 DE genes identified using DESeq2
```

##### edgeR

```{R}
for(i in 1:length(edgeR.pairwise.degenes)){
  print(names(edgeR.pairwise.degenes)[i])
  group1 <- substr(names(edgeR.pairwise.degenes)[i],1,regexpr("_vs_",names(edgeR.pairwise.degenes)[i])-1)
  group2 <- substr(names(edgeR.pairwise.degenes)[i],regexpr("_vs_",names(edgeR.pairwise.degenes)[i])+4,nchar(names(edgeR.pairwise.degenes)[i]))
  
  groups.pairwise <- groups[c(grep(group1,groups[,2]),grep(group2,groups[,2])),]
  groups.pairwise[,2] <- factor(groups.pairwise[,2],levels=unique(groups.pairwise[,2]))
  
  counts.pairwise <- counts[,colnames(counts) %in% groups.pairwise[,1]]
  
  cpm.cutoff <- 5/min(colSums(counts.pairwise)) * 1000000
  
  y <- DGEList(counts = counts.pairwise, group = groups.pairwise[,2])
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups.pairwise[,2]))
  keep.df <- as.data.frame(table(keep))
  print(paste0(length(keep.df[keep.df[,1] == F,2])," genes excluded with CPM cutoff"))
  
  y <- y[keep, , keep.lib.sizes = F]
  design <- model.matrix(~groups.pairwise[,2])
  y <- estimateDisp(y , design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
  edgeR.pairwise.degenes[[i]] <- rownames(qlf$table)[qlf$table$padj < FDRcutoff]
  print(paste0(length(edgeR.pairwise.degenes[[i]])," DE genes identified using edgeR"))


  write.table(qlf$table[qlf$table$padj < FDRcutoff,],
		  paste0(WORKING.DIR,"/human_counts_edgeR_",names(edgeR.pairwise.degenes)[i],".tsv"),
		  row.names = T,
		  col.names = NA,
		  quote = F,
		  sep = "\t")
}
```

```{R, eval = F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "151943 genes excluded with CPM cutoff"
[1] "1 DE genes identified using edgeR"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "154985 genes excluded with CPM cutoff"
[1] "0 DE genes identified using edgeR"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "142229 genes excluded with CPM cutoff"
[1] "0 DE genes identified using edgeR"
```

#### Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff

All subsequent analyses will be done using the list of genes determined to be differentially expressed by edgeR run in longitudinal mode.

##### Create sample legend for PCA and hierarchical clustering plots

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,5], y=seq(1,length(groups[,2]),1), group = groups[,1]))+
    geom_point(aes(color = groups[,5],shape=groups[,5]), size = 4)+
    scale_shape_manual(values = as.numeric(as.character(groups[!duplicated(groups[,3]),4])))+
    scale_color_manual(values = as.character(unique(groups[,3])))+
    guides(shape = guide_legend(title = "Samples", title.position = "top",nrow=4),
	     colour = guide_legend(title = "Samples", title.position = "top",nrow=4))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/human_pca_hc_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_pca_hc_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/human_pca_hc_samplekey.png)

##### Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff

```{R, fig.height=5, fig.width=8}
dendrogram <- as.data.frame(t(scale(t(log2(tpm.keep+1)))))

result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100)

structure <- get_dendro_structure(result)
dendro.data <- get_dendro_data(result)
bootstrap.positions <- get_dendro_bootstraps(dendro.data)
  
points.df <- as.data.frame(cbind(seq(1,length(structure),1),
					   structure))
dendrogroups <- groups[,2][result$hclust$order]
dendrocol <- groups[,3][result$hclust$order]
dendroshape <- groups[,4][result$hclust$order]
# dendrosize <- colSums(counts.keep)[result$hclust$order]
dendrosize <- 1

dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
  geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize, shape = dendroshape))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  scale_color_manual(values = levels(groups[,3]))+
  labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
  guides(colour = guide_legend(ncol = 2), size = F, shape = F)+
  #scale_x_discrete(limits = as.character(dendrolabels))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
	  axis.text.y = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank())

for(i in 1:length(result$edges$bp)){
  text <- round(result$edges$bp[i] * 100,0)
  dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[i,1] + 0.35, y=bootstrap.positions[i,2] + 0.01, size = 2)
}
pdf(paste0(WORKING.DIR,"/plots/human_tpm_kept_dendrogram.pdf"),
    height=5,
    width=8)
print(dendrogram.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_tpm_kept_dendrogram.png"),
    height=5,
    width=8,
	units = "in",res=300)
print(dendrogram.plot)
dev.off()

print(dendrogram.plot)
```

![image](/images/human_tpm_kept_dendrogram.png)

##### Conduct a PCA on the TPM values of all genes that passed CPM cutoff

```{R,fig.height=5,fig.width=5}
pca.df <- t(scale(t(log2(tpm.keep + 1))))
pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.keep) - 1)

pca.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = groups[,2],size = 1, shape = dendroshape))+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
	 x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
	 y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  guides(color = F,size = F, shape = F)+
  # guides(colour = guide_legend(ncol = 2))+
  scale_color_manual(values = levels(groups[,3]))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/human_tpm_kept_pca.pdf"),
    height=5,
    width=5)
print(pca.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_tpm_kept_pca.png"),
    height=5,
    width=5,
	units = "in",res=300)
print(pca.plot)
dev.off()

print(pca.plot)
```

![image](/images/human_tpm_kept_pca.png)

#### Divide differentially expressed genes into expression modules

##### Find soft power value for WGCNA

```{R}
tpm.de <- tpm[rownames(tpm) %in% rownames(edgeR.longitudinal.degenes),]

wgcna <- as.data.frame(t(tpm.de))
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(wgcna, powerVector = powers, verbose = 5)
softpower <- find_soft_power(sft)

text.color <- rep("black",length(sft$fitIndices[,1]))
text.color[which(sft$fitIndices[,1] == softpower)] <- "red"

scale_independence.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Scale Independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
  theme_bw()

mean_connectivity.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Mean Connectivity", x = "soft threshold (power)", y = "mean connectivity")+
  theme_bw()

wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
lay <- rbind(c(1,2))

pdf(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_soft_power_plot.pdf"),
    width = 10, 
    height = 5)
grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_soft_power_plot.png"),
    width = 10, 
    height = 5,
	units = "in",res=300)
grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)
dev.off()

grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)

```

![image](/images/human_tpm_de_wgcna_soft_power_plot.png)

##### Identify expression modules 

```{R}
adjacency <- adjacency(wgcna, power = softpower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");

minModuleSize <- 1
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
				     deepSplit = 2, pamRespectsDendro = FALSE,
				     minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(wgcna, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
METree = hclust(as.dist(MEDiss), method = "average")

MEDissThres = 0.25

pdf(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_merge_eigengenes_plot.pdf"),
    width = 8, 
    height = 5)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

png(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_merge_eigengenes_plot.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![image](/images/human_tpm_de_wgcna_merge_eigengenes_plot.png)

##### Merge similar expression modules

```{R}
merge = mergeCloseModules(wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

tpm.de.wgcna <- eigengene_invert_id(tpm.de, mergedColors, mergedMEs)

write.table(tpm.de.wgcna,
		paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_modules.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")

pdf(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_eigengene_dendrogram.pdf"),
    width = 8, 
    height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_eigengene_dendrogram.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

![image](/images/human_tpm_de_wgcna_eigengene_dendrogram.png)

##### Plot WGCNA expression modules as a heatmap

```{R}
tpm.de.wgcna <- wgcna_heatmap_reorder(tpm.de.wgcna)

log2tpm.de <- log2(tpm.de.wgcna[,1:(ncol(tpm.de.wgcna) - 2)] + 1)
zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
rowcol1 <- tpm.de.wgcna$module
rowcol2 <- unlist(lapply(tpm.de.wgcna$invert,function(x){if(x == F){return("grey")}else{return("black")}}))
colcol <- as.character(groups[,3])

rowsep <- get_heatmap_separators(rowcol1)
colsep <- get_heatmap_separators(colcol)

rev(sort(table(paste0(tpm.de.wgcna$module))))
```

```{R, eval = F}
         cyan       darkred         black    lightgreen        purple   greenyellow           tan   lightyellow          grey     darkgreen 
          599           250            32            26            11            11            10             3             3             3 
darkturquoise      darkgrey 
            2             2 
```

```{R}
rev(sort(table(paste0(tpm.de.wgcna$module,tpm.de.wgcna$invert))))
```

```{R, eval = F}
         cyanFALSE           cyanTRUE       darkredFALSE        darkredTRUE         blackFALSE    lightgreenFALSE           tanFALSE 
               350                249                192                 58                 31                 22                 10 
       purpleFALSE   greenyellowFALSE    greenyellowTRUE         purpleTRUE     lightgreenTRUE   lightyellowFALSE     darkgreenFALSE 
                 7                  6                  5                  4                  4                  3                  3 
         greyFALSE darkturquoiseFALSE      darkgreyFALSE           greyTRUE          blackTRUE 
                 2                  2                  2                  1                  1 
```


###### Create sample legend for WGCNA heatmap

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,5], y=seq(1,length(groups[,2]),1), group = groups[,1]))+
    geom_line(aes(color = groups[,5]), size = 4)+
    scale_color_manual(values = as.character(unique(groups[,3])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=4))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/human_hm_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_hm_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/human_hm_samplekey.png)

###### Create z-score log2TPM legend for WGCNA heatmap

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
			     colours=hmcol,
			     breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/human_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/human_hm_zscorelog2tpmkey.png)

###### Use module assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_zscorelog2tpm_module_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/human_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png)

###### Use inverse assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/human_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/human_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png)

#### Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs

```{R}
for(i in 1:length(unique(rowcol1))){
  if(nrow(tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == F,]) >= 10){
    tpm.de.wgcna.subset <- tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == F,]
    edgeR.longitudinal.degenes.subset <- edgeR.longitudinal.degenes[match(rownames(tpm.de.wgcna.subset),rownames(edgeR.longitudinal.degenes)),]
    
    write.table(edgeR.longitudinal.degenes.subset,
		    paste0(WORKING.DIR,"/human_edgeR_longitudinal_wgcna_",unique(rowcol1)[i],"F.tsv"),
		    row.names = T,
		    col.names = NA,
		    quote = F,
		    sep = "\t")
  }
  if(nrow(tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == T,]) >= 10){
    tpm.de.wgcna.subset <- tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == T,]
    edgeR.longitudinal.degenes.subset <- edgeR.longitudinal.degenes[match(rownames(tpm.de.wgcna.subset),rownames(edgeR.longitudinal.degenes)),]
    write.table(edgeR.longitudinal.degenes.subset,
		    paste0(WORKING.DIR,"/human_edgeR_longitudinal_wgcna_",unique(rowcol1)[i],"T.tsv"),
		    row.names = T,
		    col.names = NA,
		    quote = F,
		    sep = "\t")
  }
}
```

## Combine cyan F + darkred T and cyan T + darkred F module gene lists for WGCNA analyses due to their similar expression patterns

The individual files were concatenated together and the second header was manually removed.

##### Commands
```{bash, eval = F}
cat "$WORKING_DIR"/human_edgeR_longitudinal_wgcna_cyanF.tsv "$WORKING_DIR"/human_edgeR_longitudinal_wgcna_darkredT.tsv > "$WORKING_DIR"/human_edgeR_longitudinal_wgcna_cyanF_darkredT.tsv
cat "$WORKING_DIR"/human_edgeR_longitudinal_wgcna_cyanT.tsv "$WORKING_DIR"/human_edgeR_longitudinal_wgcna_darkredF.tsv > "$WORKING_DIR"/human_edgeR_longitudinal_wgcna_cyanT_darkredF.tsv
```

## Assess whether excluding 24 h samples for humans allows other clusters to be resolved

The 24 h enriched terms could reflect the cells dying rather than anything biologically significant with regards to the cultures.

### Set R inputs
```{R}
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/human_groups.tsv"
SRR2SAMPLE_MAP.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/srr_sample.map"
WORKING.DIR="Z:/EBMAL/mchung_dir/EHPYL"
```

### Load R functions

```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}

get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}

get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
							    ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}
find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(dy) < 0.05)[1]
  if(softpower == 1){
    softpower <- which(abs(dy) < 0.05)[2]
  }
  return(softpower)
}

eigengene_invert_id <- function(tpm.de, mergedColors, mergedMEs){
  tpm.de.wgcna <- tpm.de
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de[i,]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
	tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}

wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
								 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
								  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
								  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}

get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}
```

### Load packages and view sessionInfo
```{R}
library(dendextend)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(pvclust)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C				   LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WGCNA_1.68			fastcluster_1.1.25	    dynamicTreeCut_1.63-1	 vegan_2.5-5		    
 [5] lattice_0.20-35		 permute_0.9-5		   pvclust_2.0-0		   gridExtra_2.3		  
 [9] gplots_3.0.1.1		  ggplot2_3.2.0		   ggdendro_0.1-20		 FactoMineR_1.42		
[13] edgeR_3.24.3		    limma_3.38.3		    DESeq2_1.22.2		   SummarizedExperiment_1.12.0
[17] DelayedArray_0.8.0	    BiocParallel_1.16.6	   matrixStats_0.54.0	    Biobase_2.42.0		 
[21] GenomicRanges_1.34.0	  GenomeInfoDb_1.18.2	   IRanges_2.16.0		  S4Vectors_0.20.1	     
[25] BiocGenerics_0.28.0	   dendextend_1.12.0	    

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1	 htmlTable_1.13.1	 XVector_0.22.0	   base64enc_0.1-3	  rstudioapi_0.10	  bit64_0.9-7	     
 [7] mvtnorm_1.0-11	   AnnotationDbi_1.44.0   codetools_0.2-15	 splines_3.5.0	    leaps_3.0		  doParallel_1.0.15     
[13] impute_1.56.0	    robustbase_0.93-5	geneplotter_1.60.0     knitr_1.23		 zeallot_0.1.0	    Formula_1.2-3	   
[19] annotate_1.60.1	  cluster_2.0.7-1	  GO.db_3.7.0		rrcov_1.4-7		compiler_3.5.0	   backports_1.1.4	 
[25] assertthat_0.2.1	 Matrix_1.2-14	    lazyeval_0.2.2	   acepack_1.4.1	    htmltools_0.3.6	  tools_3.5.0	     
[31] gtable_0.3.0	     glue_1.3.1		 GenomeInfoDbData_1.2.0 dplyr_0.8.3		Rcpp_1.0.2		 vctrs_0.2.0	     
[37] gdata_2.18.0	     preprocessCore_1.44.0  nlme_3.1-137	     iterators_1.0.12	 xfun_0.8		   stringr_1.4.0	   
[43] gtools_3.8.1	     XML_3.98-1.20	    DEoptimR_1.0-8	   zlibbioc_1.28.0	  MASS_7.3-51.4	    scales_1.0.0	    
[49] RColorBrewer_1.1-2     memoise_1.1.0	    rpart_4.1-13	     latticeExtra_0.6-28    stringi_1.4.3	    RSQLite_2.1.2	   
[55] genefilter_1.64.0	pcaPP_1.9-73	     foreach_1.4.7	    checkmate_1.9.4	  caTools_1.17.1.2	 rlang_0.4.0	     
[61] pkgconfig_2.0.2	  bitops_1.0-6	     purrr_0.3.2		htmlwidgets_1.3	  bit_1.1-14		 tidyselect_0.2.5	
[67] robust_0.4-18.1	  magrittr_1.5	     R6_2.4.0		   Hmisc_4.2-0		fit.models_0.5-14	DBI_1.0.0		 
[73] pillar_1.4.2	     foreign_0.8-70	   withr_2.1.2		mgcv_1.8-23		survival_2.41-3	  scatterplot3d_0.3-41  
[79] RCurl_1.95-4.12	  nnet_7.3-12		tibble_2.1.3	     crayon_1.3.4	     KernSmooth_2.23-15     viridis_0.5.1	   
[85] locfit_1.5-9.1	   grid_3.5.0		 data.table_1.12.2	blob_1.2.0		 digest_0.6.20	    flashClust_1.01-2     
[91] xtable_1.8-4	     munsell_0.5.0	    viridisLite_0.3.0     
```

### Create counts data frame for human genes

```{R}
groups <- read.delim(GROUPS.PATH, header = F)
srr2sample_map <- read.delim(SRR2SAMPLE_MAP.PATH, header = F)
kallisto_output <- list.files(paste0(WORKING.DIR,"/kallisto/"), recursive = T, full.names = T, pattern = "abundance.tsv")

srr2sample_map <- srr2sample_map[srr2sample_map[,3] == "polyA",]

counts <- as.data.frame(matrix(0,
					  nrow = nrow(read.delim(kallisto_output[1],header = T)),
					  ncol = length(unique(srr2sample_map[,2]))))
rownames(counts) <- read.delim(kallisto_output[1],header = T)[,1]
colnames(counts) <- unique(srr2sample_map[,2])

for(i in 1:ncol(counts)){
  srr2sample_map.subset <- srr2sample_map[srr2sample_map[,2] == colnames(counts)[i],]
  for(j in 1:nrow(srr2sample_map.subset)){
    counts.subset <- read.delim(grep(srr2sample_map.subset[j,1],kallisto_output,value=T),header = T)
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),4]
  }
}
counts <- counts[grep("ENST",rownames(counts)),match(groups[,1],colnames(counts))]
```

### Remove 24 h samples from analysis

```{R}
counts <- counts[,!(grepl("24h",colnames(counts)))]
groups <- groups[!(grepl("24h",groups[,1])),]
```

### Create TPM data frame for human genes

```{R}
genelength <- read.delim(kallisto_output[1],header = T)[match(rownames(counts),read.delim(kallisto_output[1],header = T)[,1]),3]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

dim(tpm)
```

```{R, eval = F}
[1] 188753     15
```

### Set group levels

```{R}
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
groups[,5] <- factor(groups[,5], levels = unique(groups[,5]))
```

### Identify differentially expressed genes longitudinally

edgeR is run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

```{R}
FDRcutoff=0.05
cpm.cutoff=5

y <- DGEList(counts = counts, group = groups[,2])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,2]))
keep.df <- as.data.frame(table(keep))
print(paste0(keep.df[keep.df[,1] == F,2]," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,2])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))
qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
edgeR.longitudinal.degenes <- qlf$table[qlf$table$padj < FDRcutoff,]
print(paste0(nrow(edgeR.longitudinal.degenes)," DE genes identified using edgeR longitudinal"))

write.table(edgeR.longitudinal.degenes,
		paste0(WORKING.DIR,"/human_no24h_counts_edgeR_longitudinal.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")

counts.keep <- counts[keep,]
tpm.keep <- tpm[keep,]
```

```{R, eval = F}
[1] "159446 genes excluded with CPM cutoff"
[1] "130 DE genes identified using edgeR longitudinal
```

### Divide differentially expressed genes into expression modules

#### Find soft power value for WGCNA

```{R}
tpm.de <- tpm[rownames(tpm) %in% rownames(edgeR.longitudinal.degenes),]

wgcna <- as.data.frame(t(tpm.de))
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(wgcna, powerVector = powers, verbose = 5)
softpower <- find_soft_power(sft)

text.color <- rep("black",length(sft$fitIndices[,1]))
text.color[which(sft$fitIndices[,1] == softpower)] <- "red"

scale_independence.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Scale Independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
  theme_bw()

mean_connectivity.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Mean Connectivity", x = "soft threshold (power)", y = "mean connectivity")+
  theme_bw()

wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
lay <- rbind(c(1,2))

pdf(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_soft_power_plot.pdf"),
    width = 10, 
    height = 5)
grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_soft_power_plot.png"),
    width = 10, 
    height = 5,
	units = "in",res=300)
grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)
dev.off()

grid.arrange(grobs = wgcna_soft_power_plots,
		 widths = c(5,5),
		 heights = c(5),
		 layout_matrix = lay)

```

![image](/images/human_no24h_tpm_de_wgcna_soft_power_plot.png)

#### Identify expression modules 

```{R}
adjacency <- adjacency(wgcna, power = softpower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");

minModuleSize <- 1
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
				     deepSplit = 2, pamRespectsDendro = FALSE,
				     minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(wgcna, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
METree = hclust(as.dist(MEDiss), method = "average")

MEDissThres = 0.25

pdf(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_merge_eigengenes_plot.pdf"),
    width = 8, 
    height = 5)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

png(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_merge_eigengenes_plot.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![image](/images/human_no24h_tpm_de_wgcna_merge_eigengenes_plot.png)

#### Merge similar expression modules

```{R}
merge = mergeCloseModules(wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

tpm.de.wgcna <- eigengene_invert_id(tpm.de, mergedColors, mergedMEs)

write.table(tpm.de.wgcna,
		paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_modules.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")

pdf(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_eigengene_dendrogram.pdf"),
    width = 8, 
    height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_eigengene_dendrogram.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
			  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

![image](/images/human_no24h_tpm_de_wgcna_eigengene_dendrogram.png)

#### Plot WGCNA expression modules as a heatmap

```{R}
tpm.de.wgcna <- wgcna_heatmap_reorder(tpm.de.wgcna)

log2tpm.de <- log2(tpm.de.wgcna[,1:(ncol(tpm.de.wgcna) - 2)] + 1)
zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
rowcol1 <- tpm.de.wgcna$module
rowcol2 <- unlist(lapply(tpm.de.wgcna$invert,function(x){if(x == F){return("grey")}else{return("black")}}))
colcol <- as.character(groups[,3])

rowsep <- get_heatmap_separators(rowcol1)
colsep <- get_heatmap_separators(colcol)
```

##### Create sample legend for WGCNA heatmap

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,5], y=seq(1,length(groups[,2]),1), group = groups[,1]))+
    geom_line(aes(color = groups[,5]), size = 4)+
    scale_color_manual(values = as.character(unique(groups[,3])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=4))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/human_no24h_hm_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_no24h_hm_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/human_no24h_hm_samplekey.png)

##### Create z-score log2TPM legend for WGCNA heatmap

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
			     colours=hmcol,
			     breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/human_no24h_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_no24h_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/human_no24h_hm_zscorelog2tpmkey.png)

##### Use module assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_zscorelog2tpm_module_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol1,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/human_no24h_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png)

##### Use inverse assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/human_no24h_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
		  col=hmcol,
		  trace="none",
		  labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
		  Rowv = F,
		  Colv = F,
		  RowSideColors=rowcol2,
		  ColSideColors=colcol,
		  lhei = c(2,8),
		  breaks = seq(-3,3,by=.5),
		  rowsep = rowsep,
		  colsep = colsep,
		  dendrogram = "none")
```

![image](/images/human_no24h_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png)

### Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs

```{R}
for(i in 1:length(unique(rowcol1))){
  if(nrow(tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == F,]) >= 10){
    tpm.de.wgcna.subset <- tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == F,]
    edgeR.longitudinal.degenes.subset <- edgeR.longitudinal.degenes[match(rownames(tpm.de.wgcna.subset),rownames(edgeR.longitudinal.degenes)),]
    
    write.table(edgeR.longitudinal.degenes.subset,
		    paste0(WORKING.DIR,"/human_no24h_edgeR_longitudinal_wgcna_",unique(rowcol1)[i],"F.tsv"),
		    row.names = T,
		    col.names = NA,
		    quote = F,
		    sep = "\t")
  }
  if(nrow(tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == T,]) >= 10){
    tpm.de.wgcna.subset <- tpm.de.wgcna[tpm.de.wgcna$module == unique(rowcol1)[i] & tpm.de.wgcna$invert == T,]
    edgeR.longitudinal.degenes.subset <- edgeR.longitudinal.degenes[match(rownames(tpm.de.wgcna.subset),rownames(edgeR.longitudinal.degenes)),]
    write.table(edgeR.longitudinal.degenes.subset,
		    paste0(WORKING.DIR,"/human_no24h_edgeR_longitudinal_wgcna_",unique(rowcol1)[i],"T.tsv"),
		    row.names = T,
		    col.names = NA,
		    quote = F,
		    sep = "\t")
  }
}
```

# Conduct in vivo human RNA-Seq analysis

## Quantify human genes directly from FASTQ files

##### Inputs
```{bash, eval = F}
FASTQ_DIR=/local/aberdeen2ro/ESTAD
NUC_TRANSCRIPT_FNA="$WORKING_DIR"/references/combined_hsapiensGRCh38_hpylori26695.cds.fna
OUTPUT_DIR="$WORKING_DIR"/EHPYL/invivo/kallisto
THREADS=4
```

##### Commands
```{bash, eval = F}
for FASTQ in $(find "$FASTQ_DIR"/*R/ILLUMINA_DATA/ -name *"fastq.gz" | sed "s/_R[12].fastq.gz//g" | sort -n | uniq)
do

  SAMPLE="$(echo "$(dirname "$FASTQ")" | sed "s/\\/ILLUMINA_DATA.*//g" | sed "s/.*\\///g")"
  mkdir "$OUTPUT_DIR"/"$SAMPLE"
  
  echo -e "qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N kallisto -wd "$OUTPUT_DIR"/"$SAMPLE" -b y "$KALLISTO_BIN_DIR"/kallisto quant -t "$THREADS" -i "$NUC_TRANSCRIPT_FNA".kallisto.index -o "$OUTPUT_DIR"/"$SAMPLE"/"$(basename "$FASTQ")" "$FASTQ"_R1.fastq.gz "$FASTQ"_R2.fastq.gz"
  
done
```

## Find genes in tumor and metaplasia samples that are most up- and down-regulated relative to the tumor adjacent sample

### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/EHPYL/"
```

### Load R packages and view sessionInfo

```{R}
library(edgeR)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C				   LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] edgeR_3.24.3 limma_3.38.3

loaded via a namespace (and not attached):
[1] compiler_3.5.0  tools_3.5.0     yaml_2.2.0	Rcpp_1.0.2	grid_3.5.0	locfit_1.5-9.1  knitr_1.23	xfun_0.8	 
[9] lattice_0.20-35
```

### Create counts data frame

```{R}
counts.dir <- paste0(WORKING.DIR,"/invivo/kallisto/")
counts.files <- list.files(counts.dir, recursive = T, pattern = "abundance.tsv")

samples <- unique(gsub("\\/.*","",counts.files))

counts <- as.data.frame(matrix(0,
					 nrow=nrow(read.delim(paste0(counts.dir,"/",counts.files[1]))),
					 ncol=length(samples)))
rownames(counts) <- read.delim(paste0(counts.dir,"/",counts.files[1]))[,1]
colnames(counts) <- samples


for(i in 1:length(counts.files)){
  counts.file <- read.delim(paste0(counts.dir,"/",counts.files[i]))
  sample <- gsub("\\/.*","",counts.files[i])
  
  counts[,which(colnames(counts) == sample)] <- counts[,which(colnames(counts) == sample)] + counts.file[match(rownames(counts),counts.file[,1]),4]
}

write.table(counts,
		paste0(WORKING.DIR,"/invivo/counts.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")
```

### Check how many read counts were quantified from each organism

#### H. pylori

```{R}
colSums(counts[grep("NC_000915.1", rownames(counts)),])
```

```{R, eval = F}
  3CG_046T_R 3CG_051T_A_R   3CG_051T_R     4GB011_R 
	    10	   2487	   2739		4 
```

#### Human

```{R}
colSums(counts[grep("NC_000915.1", rownames(counts), invert = T),])
```

```{R, eval = F}
  3CG_046T_R 3CG_051T_A_R   3CG_051T_R     4GB011_R 
   158337944    158666244    145880207    161222830 
```

### Create TPM data frame

```{R}
genelength <- read.delim(paste0(counts.dir,"/",counts.files[1]))[,3]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(tpm,
		paste0(WORKING.DIR,"/invivo/tpm.tsv"),
		row.names = T,
		col.names = T,
		quote = F,
		sep = "\t")
```

### Calculate the log2TPM ratios between the two tumor samples and the metaplasia sample versus the tumor adjacent sample

The 190,198 transcripts in the analysis were pared down based on two criteria. Transcripts for downstream analyses must have a CPM value in all samples greater than 5 reads in the lowest sequenced sample of each pair.

```{R}
edgeR.pairwise.degenes <- list("3CG_046T_R_vs_3CG_051T_A_R"=c(), 
					 "3CG_051T_R_vs_3CG_051T_A_R"=c(),
					 "4GB011_R_vs_3CG_051T_A_R"=c())

for(i in 1:length(edgeR.pairwise.degenes)){
  group1 <- substr(names(edgeR.pairwise.degenes)[i],1,regexpr("_vs_",names(edgeR.pairwise.degenes)[i])-1)
  group2 <- substr(names(edgeR.pairwise.degenes)[i],regexpr("_vs_",names(edgeR.pairwise.degenes)[i])+4,nchar(names(edgeR.pairwise.degenes)[i]))
  
  groups.pairwise <- c(group1,group2)
  groups.pairwise <- factor(groups.pairwise,levels=unique(groups.pairwise))
  
  counts.pairwise <- counts[,colnames(counts) %in% groups.pairwise]
  
  cpm.cutoff <- 5/min(colSums(counts.pairwise)) * 1000000
  
  y <- DGEList(counts = counts.pairwise, group = groups.pairwise)
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= cpm.cutoff) >= ncol(counts.pairwise)
  
  counts.pairwise <- counts.pairwise[keep,]
  tpm.pairwise <- tpm[keep,colnames(tpm) %in% groups.pairwise]
  log2ratiotpm <- as.data.frame(log2(tpm.pairwise[,1]/tpm.pairwise[,2]))
  rownames(log2ratiotpm) <- rownames(counts)[keep]
  colnames(log2ratiotpm) <- "log2ratiotpm"
  
  edgeR.pairwise.degenes[[i]] <- log2ratiotpm
}
```

### For each comparison, find the top 3000 up- and down-regulated genes relative to the tumor adjacent sample

IPA allows only 3000 genes for core analysis. 

```{R}
for(i in 1:length(edgeR.pairwise.degenes)){
  log2ratiotpm.upreg <- edgeR.pairwise.degenes[[i]] 
  log2ratiotpm.upreg <- log2ratiotpm.upreg[rev(order(log2ratiotpm.upreg[,1]))[1:3000],,drop=F]
  
  log2ratiotpm.downreg <- edgeR.pairwise.degenes[[i]] 
  log2ratiotpm.downreg <- log2ratiotpm.downreg[order(log2ratiotpm.downreg[,1])[1:3000],,drop=F]
  
  write.table(log2ratiotpm.upreg,
		  paste0(WORKING.DIR,"/invivo/",names(edgeR.pairwise.degenes)[i],"_log2ratiotpm_upreg.tsv"),
		  row.names = T,
		  col.names = T,
		  quote = F,
		  sep = "\t")
  
  write.table(log2ratiotpm.downreg,
		  paste0(WORKING.DIR,"/invivo/",names(edgeR.pairwise.degenes)[i],"_log2ratiotpm_downreg.tsv"),
		  row.names = T,
		  col.names = T,
		  quote = F,
		  sep = "\t")
}
```

## Assess which in vivo human samples are most similar to the in vitro human samples

### Set R inputs

```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/EHPYL/"
INVITRO_COUNTS.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/human_counts.tsv"
INVIVO_COUNTS.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/invivo/counts.tsv"
```

### Load packages and view sessionInfo

```{R}
library(cowplot)
library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(pvclust)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C				   LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggrepel_0.8.1   pvclust_2.0-0   ggplot2_3.2.0   FactoMineR_1.42 cowplot_1.0.0  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2	     rstudioapi_0.10	cluster_2.0.7-1	knitr_1.23	     magrittr_1.5	   MASS_7.3-51.4	 
 [7] leaps_3.0		scatterplot3d_0.3-41 tidyselect_0.2.5     munsell_0.5.0	  lattice_0.20-35	colorspace_1.4-1    
[13] R6_2.4.0		 rlang_0.4.0	    dplyr_0.8.3	    tools_3.5.0	    grid_3.5.0	     gtable_0.3.0	  
[19] xfun_0.8		 withr_2.1.2	    yaml_2.2.0	     lazyeval_0.2.2	 assertthat_0.2.1     tibble_2.1.3	  
[25] crayon_1.3.4	   purrr_0.3.2	    glue_1.3.1	     labeling_0.3	   compiler_3.5.0	 pillar_1.4.2	  
[31] scales_1.0.0	   flashClust_1.01-2    pkgconfig_2.0.2     
```

### Combine counts data frame for in vivo and in vitro human data sets

```{R}
invitro_counts <- read.delim(INVITRO_COUNTS.PATH, row.names = 1)
invivo_counts <- read.delim(INVIVO_COUNTS.PATH, row.names = 1, check.names = F)

counts <- as.data.frame(cbind(invitro_counts,
					invivo_counts[match(rownames(invitro_counts),rownames(invivo_counts)),]))
```

### Create TPM data frame for in vivo and in vitro human data sets

```{R}
counts.files <- list.files(paste0(WORKING.DIR,"/invivo/kallisto/", recursive = T, pattern = "abundance.tsv")
genelength <- read.delim(paste0(WORKING.DIR,"/invivo/kallisto/",counts.files[1]))
genelength <- genelength[match(rownames(counts),genelength[,1]),3]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

counts <- counts[rowSums(counts) != 0,]
tpm <- tpm[rowSums(tpm) != 0,]
```

### Conduct a hierarchical cluster analysis on the TPM values across all in vivo and in vitro human samples

```{R, fig.height = 5, fig.width=8}
dendrogram <- as.data.frame(log10(tpm+1))

result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100)

pdf(paste0(WORKING.DIR,"/plots/human_invivo_v_invitro_dendrogram.pdf"),
    width = 6, 
    height = 4)
plot(result)
dev.off()

png(paste0(WORKING.DIR,"/plots/human_invivo_v_invitro_dendrogram.png"),
    width = 6, 
    height = 4,
	units = "in",res=300)
plot(result)
dev.off()

plot(result)
```

![image](/images/human_invivo_v_invitro_dendrogram.png)


### Conduct a PCA on the TPM values across all in vivo and in vitro human samples

```{R}
pca.df <- as.data.frame(log10(tpm+1))
pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm) - 1)

pca.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2]))+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
	 x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
	 y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  theme_bw()

pca_with_label.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2]))+
  geom_text_repel(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2],label = colnames(pca.df)),size=3)+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
	 x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
	 y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/human_invivo_v_invitro_pca.pdf"),
    width = 10, 
    height = 5)
plot_grid(plotlist=list(pca.plot, pca_with_label.plot),
	    nrow = 1,labels=c("A","B"))
dev.off()

png(paste0(WORKING.DIR,"/plots/human_invivo_v_invitro_pca.png"),
    width = 10, 
    height = 5,
    units = "in",res=300)
plot_grid(plotlist=list(pca.plot, pca_with_label.plot),
	    nrow = 1,labels=c("A","B"))
dev.off()

plot_grid(plotlist=list(pca.plot, pca_with_label.plot),
	    nrow = 1,labels=c("A","B"))
```

![image](/images/human_invivo_v_invitro_pca.png)

