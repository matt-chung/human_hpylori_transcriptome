
<!-- MarkdownTOC -->

- Set software and directory paths
	- Software
	- Directories
	- Create directories
	- Set up reference files
		- Download H. pylori and human reference files
		- Create combined H. pylori and human reference
- Confirm cag KO in H. pylori
	- Create list of SRAs to download
			- Commands
	- Download H. pylori genomic FASTQs from the SRA
			- Inputs
			- Commands
	- Align H. pylori genomic FASTQs to H. pylori reference
			- Inputs
			- Commands
	- Assess depth over H. pylori genome
			- Inputs
			- Commands
	- Visualize depth over H. pylori genome
		- Set R inputs
		- Load R packages and view sessionInfo
		- Construct data frame of the depth at H. pylori genome from positions 500,000-600,000
		- Calculate average depth across the H. pylori knockout region at positions 547,165-584,551
		- Create plot data frame
		- Plot depth in H. pylori genome from positions 500,000-600,000
- Conduct in vitro H. pylori and human RNA-Seq analysis
	- Create list of SRAs to download
			- Inputs
			- Commands
	- Calculate gene counts for H. pylori genes
		- Align FASTQs with splicing to combined human and H. pylori reference
			- Inputs
			- Commands
		- Find reads that mapped to H. pylori
			- Inputs
			- Commands
		- Subset FASTQs to only contain H. pylori-mapping reads
			- Inputs
			- Commands
		- Align subset FASTQs with no splicing to combined human and H. pylori reference
			- Inputs
			- Commands
		- Sort BAM files
			- Inputs
			- Commands
		- Index BAM files
			- Inputs
			- Commands
		- Quantify H. pylori genes from BAM files
			- Inputs
			- Commands
	- Calculate gene counts for human genes
		- Quantify human genes directly from FASTQ files
			- Inputs
			- Commands

<!-- /MarkdownTOC -->

# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software
```{bash, eval = F}
JULIA_DEPOT_PATH=/home/mattchung/.julia

JULIA_BIN_DIR=/usr/local/bin

BOWTIE2_BIN_DIR=/usr/local/packages/bowtie2-2.3.4.3
FADU_BIN_DIR=/home/mattchung/scripts/FADU
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
KALLISTO_BIN_DIR=/usr/local/packages/kallisto-0.45.0
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
```

## Directories
```{bash, eval = F}
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

## Create list of SRAs to download
```{bash, eval = F}
vim  "$WORKING_DIR"/PRJNA378649.srr.list
```

```{bash, eval = F}
SRR5422043
SRR5422044
SRR5422045
SRR5422046
SRR5422047
SRR5422048
SRR5422049
SRR5422050
SRR5422741
SRR5422742
SRR5422743
SRR5422744
SRR5422746
SRR5422747
SRR5422748
SRR5422749
SRR5422750
SRR5422751
SRR5422752
SRR5422753
SRR5422754
SRR5422755
SRR5422756
SRR5422831
SRR5422832
SRR5422833
SRR5422834
SRR5422835
SRR5422836
SRR5422837
SRR5422838
SRR5422839
SRR5422840
SRR5422841
SRR5422842
SRR5422843
SRR5422844
SRR5422845
SRR5422846
SRR5422847
SRR5422848
SRR5422849
SRR5422850
SRR5422851
SRR5422852
SRR5422853
SRR5422854
SRR5422855
SRR5422856
SRR5422857
SRR5422858
SRR5422927
SRR5422928
SRR5422945
SRR5422946
SRR5422947
SRR5422948
SRR5423058
SRR5423059
SRR5423069
SRR5423070
SRR5423599
SRR5423603
SRR5423604
SRR5423632
SRR5423636
SRR5423638
SRR5423655
SRR5423658
SRR5423659
SRR5424793
SRR5424794
SRR5424795
SRR5424796
SRR5424797
SRR5424826
```

##### Inputs
```{bash, eval = F}
READS_DIR="$WORKING_DIR"/reads
SRR_ID_LIST="$WORKING_DIR"/PRJNA378649.srr.list
```

##### Commands
```{bash, eval = F}
while read SRR_ID
do
	qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" -O "$READS_DIR"
done < "$SRR_ID_LIST"
```

## Calculate gene counts for H. pylori genes

### Align FASTQs with splicing to combined human and H. pylori reference

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

### Find reads that mapped to H. pylori

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

### Subset FASTQs to only contain H. pylori-mapping reads

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

### Align subset FASTQs with no splicing to combined human and H. pylori reference

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

### Sort BAM files

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

### Index BAM files

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

### Quantify H. pylori genes from BAM files

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

## Calculate gene counts for human genes
### Quantify human genes directly from FASTQ files
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

