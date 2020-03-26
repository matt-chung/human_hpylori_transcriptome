
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
- [Conduct in vitro H. pylori and human RNA-Seq analysis](#conduct-in-vitro-h-pylori-and-human-rna-seq-analysis)
	- [Create SRA ID and sample map file](#create-sra-id-and-sample-map-file)
	- [Download FASTQs from SRA](#download-fastqs-from-sra)
	- [Find](#find)
	- [Quantify transcipts](#quantify-transcipts)
		- [H. pylori](#h-pylori)
			- [Align FASTQs with splicing to combined human and H. pylori reference](#align-fastqs-with-splicing-to-combined-human-and-h-pylori-reference)
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
		- [Human](#human-2)
			- [Set R inputs](#set-r-inputs-1)
			- [Load R functions](#load-r-functions)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo)
			- [Create counts data frame for human genes](#create-counts-data-frame-for-human-genes)
			- [Create TPM data frame for human genes](#create-tpm-data-frame-for-human-genes)
			- [Set group levels](#set-group-levels)
			- [Conduct saturation analysis for human genes](#conduct-saturation-analysis-for-human-genes)
			- [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally)
			- [Identify differentially expressed genes with pairwise comparisons](#identify-differentially-expressed-genes-with-pairwise-comparisons)
			- [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff)
			- [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules)
			- [Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs](#create-a-list-of-genes-and-their-logfc-and-fdr-values-for-each-wgcna-module-as-ipa-inputs)
	- [Identify differentially expressed H. pylori genes across the data set](#identify-differentially-expressed-h-pylori-genes-across-the-data-set)
		- [Set R inputs](#set-r-inputs-2)
		- [Load R functions](#load-r-functions-1)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-1)
		- [Create counts data frame for H. pylori genes](#create-counts-data-frame-for-h-pylori-genes)
		- [Exclude genes with 0 unique positions](#exclude-genes-with-0-unique-positions)
		- [Create TPM data frame for H. pylori genes](#create-tpm-data-frame-for-h-pylori-genes)
		- [Set group levels](#set-group-levels-1)
		- [Conduct saturation analysis for human genes](#conduct-saturation-analysis-for-human-genes-1)
		- [Exclude low count samples](#exclude-low-count-samples)
		- [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally-1)
			- [DESeq2](#deseq2)
			- [edgeR](#edger)
		- [Identify differentially expressed genes with pairwise comparisons](#identify-differentially-expressed-genes-with-pairwise-comparisons-1)
			- [DESeq2](#deseq2-1)
			- [edgeR](#edger-1)
		- [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff-1)
			- [Create sample legend for PCA and hierarchical clustering plots](#create-sample-legend-for-pca-and-hierarchical-clustering-plots)
			- [Conduct a PCA on the TPM values of all genes that passed CPM cutoff](#conduct-a-pca-on-the-tpm-values-of-all-genes-that-passed-cpm-cutoff)
			- [Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff](#conduct-a-hierarchical-cluster-analysis-on-the-tpm-values-of-all-genes-that-passed-cpm-cutoff)
		- [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules-1)
			- [Find soft power value for WGCNA](#find-soft-power-value-for-wgcna)
			- [Identify expression modules](#identify-expression-modules)
			- [Merge similar expression modules](#merge-similar-expression-modules)
			- [Plot WGCNA expression modules as a heatmap](#plot-wgcna-expression-modules-as-a-heatmap)
		- [Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs](#create-a-list-of-genes-and-their-logfc-and-fdr-values-for-each-wgcna-module-as-ipa-inputs-1)
		- [Plot a heatmap of the module containing the cag pathogenicity genes](#plot-a-heatmap-of-the-module-containing-the-cag-pathogenicity-genes)

<!-- /MarkdownTOC -->

# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software
```{bash, eval = F}
JULIA_DEPOT_PATH=/home/mattchung/.julia
PYTHON_LIB_PATH=/usr/local/packages/python-3.5/lib

JULIA_BIN_DIR=/usr/local/bin
R_BIN_DIR=/usr/local/packages/r-3.6.0/bin

BOWTIE2_BIN_DIR=/usr/local/packages/bowtie2-2.3.4.3
FADU_BIN_DIR=/home/mattchung/scripts/FADU
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
INTERPROSCAN_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/interproscan-5.34-73.0/
KALLISTO_BIN_DIR=/usr/local/packages/kallisto-0.45.0
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
```

## Directories
```{bash, eval = F}
SCRIPTS_DIR=/home/mattchung/scripts/
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/EHPYL
```

## Create directories
```{bash, eval = F}
mkdir -p "$WORKING_DIR"/bam
mkdir -p "$WORKING_DIR"/fadu
mkdir -p "$WORKING_DIR"/kallisto
mkdir -p "$WORKING_DIR"/plots
mkdir -p "$WORKING_DIR"/reads
mkdir -p "$WORKING_DIR"/references
```

## Set up reference files

### Download H. pylori and human reference files

```{bash, eval = F}
wget -O "$WORKING_DIR"/references/hpylori26995.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/525/GCF_000008525.1_ASM852v1/GCF_000008525.1_ASM852v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/hpylori26995.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/525/GCF_000008525.1_ASM852v1/GCF_000008525.1_ASM852v1_genomic.gff.gz
gunzip "$WORKING_DIR"/references/hpylori26995.fna.gz
gunzip "$WORKING_DIR"/references/hpylori26995.gff.gz

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
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.2.0 reshape_0.8.8

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2       rstudioapi_0.10  magrittr_1.5     tidyselect_0.2.5 munsell_0.5.0    colorspace_1.4-1 R6_2.4.0         rlang_0.4.0     
 [9] stringr_1.4.0    plyr_1.8.4       dplyr_0.8.3      tools_3.5.0      grid_3.5.0       gtable_0.3.0     withr_2.1.2      yaml_2.2.0      
[17] lazyeval_0.2.2   assertthat_0.2.1 tibble_2.1.3     crayon_1.3.4     purrr_0.3.2      reshape2_1.4.3   glue_1.3.1       labeling_0.3    
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

# Conduct in vitro H. pylori and human RNA-Seq analysis

## Create SRA ID and sample map file
```{bash, eval = F}
vim "$WORKING_DIR"/srr_sample.map
```

```{bash, eval = F}
SRR5423603      N87_a
SRR5423604      N87_a
SRR5422947      N87_b
SRR5422948      N87_b
SRR5423599      N87_c
SRR5422749      N87_Hp_cag_2h_a
SRR5422857      N87_Hp_cag_2h_a
SRR5422858      N87_Hp_cag_2h_a
SRR5422746      N87_Hp_cag_2h_b
SRR5423069      N87_Hp_cag_2h_b
SRR5423070      N87_Hp_cag_2h_b
SRR5423636      N87_Hp_cag_2h_c
SRR5422747      N87_Hp_cag_4h_a
SRR5422856      N87_Hp_cag_4h_a
SRR5422748      N87_Hp_cag_4h_b
SRR5423058      N87_Hp_cag_4h_b
SRR5422946      N87_Hp_cag_4h_c
SRR5422750      N87_Hp_cag_24h_a
SRR5422751      N87_Hp_cag_24h_a
SRR5422752      N87_Hp_cag_24h_a
SRR5423659      N87_Hp_cag_24h_a
SRR5422753      N87_Hp_cag_24h_b
SRR5422754      N87_Hp_cag_24h_b
SRR5422755      N87_Hp_cag_24h_b
SRR5423632      N87_Hp_cag_24h_b
SRR5423658      N87_Hp_cag_24h_c
SRR5422756      N87_Hp_ko_2h_a
SRR5423638      N87_Hp_ko_2h_a
SRR5422833      N87_Hp_ko_2h_b
SRR5423059      N87_Hp_ko_2h_b
SRR5422928      N87_Hp_ko_2h_c
SRR5422831      N87_Hp_ko_4h_a
SRR5423655      N87_Hp_ko_4h_a
SRR5422832      N87_Hp_ko_4h_b
SRR5422927      N87_Hp_ko_4h_b
SRR5422945      N87_Hp_ko_4h_c
SRR5422834      N87_Hp_ko_24h_b
SRR5422835      N87_Hp_ko_24h_b
SRR5422836      N87_Hp_ko_24h_b
SRR5424793      N87_Hp_ko_24h_b
SRR5422837      N87_Hp_ko_24h_c
SRR5422838      N87_Hp_ko_24h_c
SRR5422839      N87_Hp_ko_24h_c
SRR5424794      N87_Hp_ko_24h_c
SRR5422043      Hp_cag_0h_a
SRR5422045      Hp_cag_0h_b
SRR5422044      Hp_cag_0h_c
SRR5422048      Hp_cag_2h_a
SRR5422046      Hp_cag_2h_b
SRR5422047      Hp_cag_2h_c
SRR5422049      Hp_cag_4h_a
SRR5422050      Hp_cag_4h_b
SRR5422741      Hp_cag_4h_c
SRR5422742      Hp_cag_24h_a
SRR5422743      Hp_cag_24h_b
SRR5422744      Hp_cag_24h_c
SRR5422840      Hp_ko_0h_a
SRR5422841      Hp_ko_0h_b
SRR5422842      Hp_ko_0h_c
SRR5422843      Hp_ko_0h_c
SRR5422844      Hp_ko_0h_c
SRR5422845      Hp_ko_2h_a
SRR5422846      Hp_ko_2h_b
SRR5422847      Hp_ko_2h_c
SRR5422848      Hp_ko_4h_a
SRR5422849      Hp_ko_4h_b
SRR5422850      Hp_ko_4h_c
SRR5422851      Hp_ko_24h_a
SRR5422852      Hp_ko_24h_b
SRR5422853      Hp_ko_24h_c
SRR5422854      Hp_ko_24h_c
SRR5422855      Hp_ko_24h_c
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
	qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" -O "$READS_DIR"
done
```

## Find 


## Quantify transcipts

### H. pylori

#### Align FASTQs with splicing to combined human and H. pylori reference

##### Inputs
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
	echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$REF_FNA" -1 "$READS_DIR"/"$SRR"_1.fastq -2 "$READS_DIR"/"$SRR"_2.fastq | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
done < "$SRR_ID_LIST"
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

```{R}
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
 [1] WGCNA_1.68                  fastcluster_1.1.25          dynamicTreeCut_1.63-1      
 [4] vegan_2.5-5                 lattice_0.20-35             permute_0.9-5              
 [7] pvclust_2.0-0               gridExtra_2.3               gplots_3.0.1.1             
[10] ggplot2_3.2.0               ggdendro_0.1-20             FactoMineR_1.42            
[13] edgeR_3.24.3                limma_3.38.3                DESeq2_1.22.2              
[16] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6        
[19] matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0       
[22] GenomeInfoDb_1.18.2         IRanges_2.16.0              S4Vectors_0.20.1           
[25] BiocGenerics_0.28.0         dendextend_1.12.0          

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1       htmlTable_1.13.1       XVector_0.22.0         base64enc_0.1-3       
 [5] rstudioapi_0.10        bit64_0.9-7            mvtnorm_1.0-11         AnnotationDbi_1.44.0  
 [9] codetools_0.2-15       splines_3.5.0          leaps_3.0              doParallel_1.0.15     
[13] impute_1.56.0          robustbase_0.93-5      geneplotter_1.60.0     knitr_1.23            
[17] zeallot_0.1.0          Formula_1.2-3          annotate_1.60.1        cluster_2.0.7-1       
[21] GO.db_3.7.0            rrcov_1.4-7            compiler_3.5.0         backports_1.1.4       
[25] assertthat_0.2.1       Matrix_1.2-14          lazyeval_0.2.2         acepack_1.4.1         
[29] htmltools_0.3.6        tools_3.5.0            gtable_0.3.0           glue_1.3.1            
[33] GenomeInfoDbData_1.2.0 dplyr_0.8.3            Rcpp_1.0.2             vctrs_0.2.0           
[37] gdata_2.18.0           preprocessCore_1.44.0  nlme_3.1-137           iterators_1.0.12      
[41] xfun_0.8               stringr_1.4.0          gtools_3.8.1           XML_3.98-1.20         
[45] DEoptimR_1.0-8         zlibbioc_1.28.0        MASS_7.3-51.4          scales_1.0.0          
[49] RColorBrewer_1.1-2     yaml_2.2.0             memoise_1.1.0          rpart_4.1-13          
[53] latticeExtra_0.6-28    stringi_1.4.3          RSQLite_2.1.2          genefilter_1.64.0     
[57] pcaPP_1.9-73           foreach_1.4.7          checkmate_1.9.4        caTools_1.17.1.2      
[61] rlang_0.4.0            pkgconfig_2.0.2        bitops_1.0-6           purrr_0.3.2           
[65] htmlwidgets_1.3        labeling_0.3           bit_1.1-14             tidyselect_0.2.5      
[69] robust_0.4-18.1        magrittr_1.5           R6_2.4.0               fit.models_0.5-14     
[73] Hmisc_4.2-0            DBI_1.0.0              pillar_1.4.2           foreign_0.8-70        
[77] withr_2.1.2            mgcv_1.8-23            survival_2.41-3        scatterplot3d_0.3-41  
[81] RCurl_1.95-4.12        nnet_7.3-12            tibble_2.1.3           crayon_1.3.4          
[85] KernSmooth_2.23-15     viridis_0.5.1          locfit_1.5-9.1         grid_3.5.0            
[89] data.table_1.12.2      blob_1.2.0             digest_0.6.20          flashClust_1.01-2     
[93] xtable_1.8-4           munsell_0.5.0          viridisLite_0.3.0
```

#### Create counts data frame for human genes

```{R}
groups <- read.delim(GROUPS.PATH, header = F)
srr2sample_map <- read.delim(SRR2SAMPLE_MAP.PATH, header = F)
kallisto_output <- list.files(paste0(WORKING.DIR,"/kallisto/"), recursive = T, full.names = T, pattern = "abundance.tsv")

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
counts <- counts[grep("ENST",rownames(counts)),match(groups1[,1],colnames(counts))]

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

png(paste0(WORKING.DIR,"/plots/human_rarefication_plot.pdf"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

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

```{R, eval - F}
[1] "69453 genes excluded with CPM cutoff"
[1] "1604 DE genes identified using DESeq2 longitudinal"
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
[1] "919 DE genes identified using edgeR longitudinal
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

```{R, eval - F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "100628 genes excluded with CPM cutoff"
[1] "0 DE genes identified using DESeq2"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "108194 genes excluded with CPM cutoff"
[1] "0 DE genes identified using DESeq2"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "108990 genes excluded with CPM cutoff"
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

```{R, eval - F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "151943 genes excluded with CPM cutoff"
[1] "1 DE genes identified using edgeR"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "154985 genes excluded with CPM cutoff"
[1] "0 DE genes identified using edgeR"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "142229 genes excluded with CPM cutoff"
[1] "2 DE genes identified using edgeR"
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

## Identify differentially expressed H. pylori genes across the data set

### Set R inputs
```{R}
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/hpylori_groups.tsv"
SRR2SAMPLE_MAP.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/srr_sample.map"
GFF3_MAP.PATH <- "Z:/EBMAL/mchung_dir/EHPYL/References/hpylori26995.gff.map"
WORKING.DIR="Z:/EBMAL/mchung_dir/EHPYL/"
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

```{R}
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
 [1] vegan_2.5-5                 lattice_0.20-35             permute_0.9-5              
 [4] pvclust_2.0-0               gridExtra_2.3               gplots_3.0.1.1             
 [7] ggplot2_3.2.0               ggdendro_0.1-20             FactoMineR_1.42            
[10] edgeR_3.24.3                limma_3.38.3                DESeq2_1.22.2              
[13] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6        
[16] matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0       
[19] GenomeInfoDb_1.18.2         IRanges_2.16.0              S4Vectors_0.20.1           
[22] BiocGenerics_0.28.0         dendextend_1.12.0          

loaded via a namespace (and not attached):
 [1] nlme_3.1-137           bitops_1.0-6           bit64_0.9-7            RColorBrewer_1.1-2    
 [5] tools_3.5.0            backports_1.1.4        R6_2.4.0               KernSmooth_2.23-15    
 [9] rpart_4.1-13           mgcv_1.8-23            Hmisc_4.2-0            DBI_1.0.0             
[13] lazyeval_0.2.2         colorspace_1.4-1       nnet_7.3-12            withr_2.1.2           
[17] tidyselect_0.2.5       bit_1.1-14             compiler_3.5.0         htmlTable_1.13.1      
[21] flashClust_1.01-2      caTools_1.17.1.2       scales_1.0.0           checkmate_1.9.4       
[25] genefilter_1.64.0      stringr_1.4.0          digest_0.6.20          foreign_0.8-70        
[29] XVector_0.22.0         base64enc_0.1-3        pkgconfig_2.0.2        htmltools_0.3.6       
[33] htmlwidgets_1.3        rlang_0.4.0            rstudioapi_0.10        RSQLite_2.1.2         
[37] gtools_3.8.1           acepack_1.4.1          dplyr_0.8.3            RCurl_1.95-4.12       
[41] magrittr_1.5           GenomeInfoDbData_1.2.0 Formula_1.2-3          leaps_3.0             
[45] Matrix_1.2-14          Rcpp_1.0.2             munsell_0.5.0          viridis_0.5.1         
[49] scatterplot3d_0.3-41   stringi_1.4.3          yaml_2.2.0             MASS_7.3-51.4         
[53] zlibbioc_1.28.0        grid_3.5.0             blob_1.2.0             gdata_2.18.0          
[57] crayon_1.3.4           splines_3.5.0          annotate_1.60.1        locfit_1.5-9.1        
[61] zeallot_0.1.0          knitr_1.23             pillar_1.4.2           geneplotter_1.60.0    
[65] XML_3.98-1.20          glue_1.3.1             latticeExtra_0.6-28    data.table_1.12.2     
[69] vctrs_0.2.0            gtable_0.3.0           purrr_0.3.2            assertthat_0.2.1      
[73] xfun_0.8               xtable_1.8-4           survival_2.41-3        viridisLite_0.3.0     
[77] tibble_2.1.3           AnnotationDbi_1.44.0   memoise_1.1.0          cluster_2.0.7-1   
```

### Create counts data frame for H. pylori genes

```{R}
groups <- read.delim(GROUPS.PATH, header = F)
srr2sample_map <- read.delim(SRR2SAMPLE_MAP.PATH, header = F)
fadu_output <- list.files(paste0(WORKING.DIR,"/fadu/"), recursive = T, full.names = T, pattern = "counts.txt")

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
N87_Hp_cag_24h_c  N87_Hp_cag_4h_c   N87_Hp_ko_4h_c  N87_Hp_cag_2h_c   N87_Hp_ko_2h_c     Hp_cag_24h_c 
         3071.71         16520.82         31907.22         34886.26         60839.53       5139048.19 
     Hp_cag_4h_a     Hp_cag_24h_b       Hp_ko_2h_c      Hp_ko_24h_b     Hp_cag_24h_a      Hp_ko_24h_a 
      6026149.54       7117275.66       7516840.17       7732774.68       8322814.38       8382443.90 
      Hp_ko_0h_b      Hp_cag_2h_b       Hp_ko_4h_b       Hp_ko_4h_a      Hp_cag_4h_c       Hp_ko_0h_a 
      9189260.42       9398914.86       9554572.16       9637179.54      11127451.44      11153378.47 
      Hp_ko_4h_c      Hp_cag_0h_c      Hp_cag_0h_b       Hp_ko_2h_a      Hp_cag_2h_c      Hp_ko_24h_c 
     11689215.36      11795833.81      12161552.73      12337556.79      12589538.55      12856729.85 
     Hp_cag_0h_a      Hp_cag_2h_a       Hp_ko_2h_b       Hp_ko_0h_c      Hp_cag_4h_b   N87_Hp_ko_4h_b 
     13300948.58      13487127.37      13681967.48      14458854.88      15245199.89      16447977.90 
N87_Hp_cag_24h_b   N87_Hp_ko_4h_a  N87_Hp_cag_4h_b  N87_Hp_cag_2h_b  N87_Hp_cag_2h_a  N87_Hp_cag_4h_a 
     18112648.30      19388385.65      22037640.22      27078238.33      30870027.59      31353003.53 
N87_Hp_cag_24h_a  N87_Hp_ko_24h_b   N87_Hp_ko_2h_a   N87_Hp_ko_2h_b  N87_Hp_ko_24h_c 
     34309415.00      37945751.47      43733495.35      44622923.10      59278742.54 
```

### Exclude genes with 0 unique positions

```{R}
genelength <- read.delim(fadu_output[1],header = T)[match(rownames(counts),read.delim(fadu_output[1],header = T)[,1]),2]

counts <- counts[genelength != 0,]
genelength <- genelength[genelength != 0]
```

### Create TPM data frame for H. pylori genes
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
[1] 1445   41
```

### Set group levels

```{R}
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
groups[,5] <- factor(groups[,5], levels = unique(groups[,5]))
```

### Conduct saturation analysis for human genes

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

png(paste0(WORKING.DIR,"/plots/hpylori_rarefication_plot.pdf"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

### Exclude low count samples

From the rarefaction analysis, the following samples were excluded (<100k gene counts): "N87_Hp_cag_2h_c", "N87_Hp_cag_4h_c", "N87_Hp_cag_24h_c", "N87_Hp_ko_2h_c", "N87_Hp_ko_4h_c"

```{R}
exclude.samples <- c("N87_Hp_cag_2h_c","N87_Hp_cag_4h_c","N87_Hp_cag_24h_c","N87_Hp_ko_2h_c","N87_Hp_ko_4h_c")
groups <- groups[!(groups[,1] %in% exclude.samples),]

groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
groups[,5] <- factor(groups[,5], levels = unique(groups[,5]))

counts <- counts[,colnames(counts) %in% groups[,1]]
tpm <- tpm[,colnames(tpm) %in% groups[,1]]
```

### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

#### DESeq2

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

```{R, eval - F}
[1] "0 genes excluded with CPM cutoff"
[1] "1444 DE genes identified using DESeq2 longitudinal"
```

#### edgeR

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

### Identify differentially expressed genes with pairwise comparisons

#### DESeq2

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

```{R, eval - F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "0 genes excluded with CPM cutoff"
[1] "814 DE genes identified using DESeq2"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "2 genes excluded with CPM cutoff"
[1] "1099 DE genes identified using DESeq2"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "0 genes excluded with CPM cutoff"
[1] "994 DE genes identified using DESeq2"

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

#### edgeR

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

```{R, eval - F}
[1] "N87_Hp_cag_2h_vs_N87_Hp_ko_2h"
[1] "2 genes excluded with CPM cutoff"
[1] "818 DE genes identified using edgeR"

[1] "N87_Hp_cag_4h_vs_N87_Hp_ko_4h"
[1] "2 genes excluded with CPM cutoff"
[1] "1131 DE genes identified using edgeR"

[1] "N87_Hp_cag_24h_vs_N87_Hp_ko_24h"
[1] "2 genes excluded with CPM cutoff"
[1] "972 DE genes identified using edgeR"

[1] "Hp_cag_0h_vs_Hp_ko_0h"
[1] "0 genes excluded with CPM cutoff"
[1] "943 DE genes identified using edgeR"

[1] "Hp_cag_2h_vs_Hp_ko_2h"
[1] "0 genes excluded with CPM cutoff"
[1] "1009 DE genes identified using edgeR"

[1] "Hp_cag_4h_vs_Hp_ko_4h"
[1] "2 genes excluded with CPM cutoff"
[1] "838 DE genes identified using edgeR"

[1] "Hp_cag_24h_vs_Hp_ko_24h"
[1] "1 genes excluded with CPM cutoff"
[1] "1062 DE genes identified using edgeR"
```
### Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff

#### Create sample legend for PCA and hierarchical clustering plots

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
#### Conduct a PCA on the TPM values of all genes that passed CPM cutoff

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
#### Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff

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

#### Identify expression modules 

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

#### Merge similar expression modules

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
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=5))+
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

##### Use module assigments as the row color bar

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
              rowsep = rowsep1,
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

##### Use inverse assigments as the row color bar

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

### Create a list of genes and their logFC and FDR values for each WGCNA module as IPA inputs

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

### Plot a heatmap of the module containing the cag pathogenicity genes

##### Create heatmap row labels containing gene names and encoded products 

```{R}
module <- "darkslateblue"

zscore.log2tpm.cag <- zscore.log2tpm.de[which(rowcol1 == module),]

gff3_map <- read.delim(GFF3_MAP.PATH, comment.char = "#", header = T)
gff3_map$old_locus_tag <- gsub("%2C.*","",gff3_map$old_locus_tag)

labrow <- paste0(rownames(zscore.log2tpm.cag),": ", gff3_map$product[gff3_map$old_locus_tag %in% rownames(zscore.log2tpm.cag)])
```

##### Use module assigments as the row color bar

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

##### Use inverse assigments as the row color bar

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
