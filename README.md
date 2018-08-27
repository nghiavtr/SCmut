# SCmut: A method to detect cell-specific mutation from single-cell RNA-sequencing

## 1. Introduction
SCmut is a novel and robust statistical method for cell-specific somatic mutation detection from single-cell RNA-sequencing. SCmut requies RNA-sequencing data of single cells and bulk-cell DNA-sequencing (e.g whole exome sequencing - WES) of matched samples (tumor and normal). If the DNA-sequencing data are not available, the list of somatic mutations can be used.

Software requirements for SCmut:
- Java 1.8 or higher
- R 3.2 or higher
- Samtools 1.3 or higher
- Picard 2.3 or higher
- VarScan 2.3.7 or higher
- GATKAnalysisTK 3.6 or higher

Annotation reference: Several common annotations and databases are required for alignment, preprocessing and variant calling.  For example, we can use the annotation b37 of the reference genome and known variant sites from phase I of 1000 Genomes Project and dbSNP-138 from broadinstitute.org.

## 2. Pre-processing for RNA-seq and DNA-seq data
### 2.1 Alignment
Start from FASTQ files, all data samples are input for read alignment to get BAM files. Since RNA-seq and DNA-seq have different biological data structures, tools for the read alignment should be selected properly. For example, we use Tophat for RNA-seq and BWA for DNA-seq.

### 2.2 Preparation
In order to reduce the number of false positives, we use standard BAM processing steps suggested from GATK pipeline. Assume that we do processing for a input BAM file input.bam, then for user's convenience, we aim to present the following steps in the copy-and-paste manner. 
```sh
bam_fn=input.bam
```
The workflow can be used to both DNAseq and RNAseq, but an extra step is required for RNA-seq data. We need to specify the type of data ($seqType) in advance:
```sh
seqType="RNA"
```
Assume that we already downloaded the reference genome and known variant sites from phase I of 1000 Genomes Project and dbSNP-138 from broadinstitute.org. We use annotation b37 version 2.8 for this example.
```sh
genomeFasta_b37="human_g1k_v37.fasta"
refKnown1="Mills_and_1000G_gold_standard.indels.b37.vcf"
refKnown2="1000G_phase1.indels.b37.vcf"
dbsnp="dbsnp_138.b37.vcf"
cosmic="b37_cosmic_v54_120711.vcf"
```

#### 2.3 Remove duplicate reads
```sh
rmdup_fn=$(echo $(basename $bam_fn)$"_rmDup.bam")
rmdup_metric_fn=$(echo $(basename $bam_fn)$"_rmDup.metric")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar picard.jar MarkDuplicates I=$bam_fn O=$rmdup_fn REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$rmdup_metric_fn
```

#### 2.4 Assign a read group ID
```sh
read_group_id=$(echo $(basename $bam_fn))
readgroup_fn=$(echo $(basename $bam_fn)$"_rg.bam")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar picard.jar AddOrReplaceReadGroups I=$rmdup_fn O=$readgroup_fn SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=$read_group_id RGLB=RNA RGPL=Illumina RGPU=illumina RGSM=$read_group_id
```

#### 2.5 Reorder a bam file
```sh
reorder_fn=$(echo $(basename $bam_fn)$"_reorder.bam")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar picard.jar ReorderSam VALIDATION_STRINGENCY=SILENT I=$readgroup_fn O=$reorder_fn REFERENCE=$genomeFasta_b37
samtools index $reorder_fn
```

#### 2.6 Extra step for for RNA-seq: Split'N'Trim and reassign mapping qualities
```sh
Realignment_inFn=$reorder_fn
if [ "$seqType" == "RNA" ]; then
  split_fn=$(echo $(basename $bam_fn)$"_split.bam")
  java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomeFasta_b37 -I $reorder_fn -o $split_fn -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
  Realignment_inFn=$split_fn  
fi
```
Thus, for RNA-seq data we get input of RealingerTargetCreator step from the output of SplitNCigarReads (Realignment_inFn=$split_fn), for DNA-seq data we use the output from Reorder step (Realignment_inFn=$reorder_fn). 

#### 2.7 GATK - RealingerTargetCreator
```sh
RealignerTargetCreator_fn=$(echo $(basename $bam_fn)$"_forIndelRealinger.intervals")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $genomeFasta_b37 -I $Realignment_inFn -known $refKnown1 -known $refKnown2 -o $RealignerTargetCreator_fn -U ALLOW_SEQ_DICT_INCOMPATIBILITY
```
#### 2.8 GATK - IndelRealinger
```sh
IndelRealinger_fn=$(echo $(basename $bam_fn)$"_realigned.bam")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar GenomeAnalysisTK.jar -T IndelRealigner -R $genomeFasta_b37 -I $Realignment_inFn -targetIntervals $RealignerTargetCreator_fn -o $IndelRealinger_fn -known $refKnown1 -known $refKnown2
```
#### 2.9 GATK - BaseRecalibrator
```sh
    BaseRecalibrator_fn=$(echo $(basename $bam_fn)$"_recal.grp")
    java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I $IndelRealinger_fn -R $genomeFasta_b37 -knownSites $refKnown1 -knownSites $refKnown2 -knownSites $dbsnp -o $BaseRecalibrator_fn
```
#### 2.10 GATK - PrintReads
```sh
PrintReads_fn=$(echo $(basename $bam_fn)$"_recalibrated.bam")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar GenomeAnalysisTK.jar -T PrintReads -R $genomeFasta_b37 -I $IndelRealinger_fn -BQSR $BaseRecalibrator_fn -o $PrintReads_fn
```

## 3. Somatic mutation detection from matched samples
If the DNA-sequencing data of matched samples are available, we can discover the somatic mutations. Assume that after preprocessing we obtain two BAM files of germline (normal) sample (normal.bam) and tumor sample (tumor.bam) from DNA-seq data. Any somatic mutation detection tools such as Mutect or VarScan can be used to discover somatic mutations from the matched samples. The following example is for Mutect:
```sh
DNA_g_fn="normal.bam"
DNA_t_fn="tumor.bam"
mutation_fn="SomaticMutation.vcf")
java -Xmx16g -Djava.io.tmpdir=`pwd`/tmp -jar muTect.jar --analysis_type MuTect --reference_sequence $genomeFasta_b37 --cosmic $cosmic --dbsnp $dbsnp --input_file:normal $DNA_g_fn --input_file:tumor $DNA_t_fn --out $mutation_fn
```

## 4. Variant calling of multiple files from both RNA-seq and DNA-seq data
We list the names of processed BAM files of RNA single cells and (if available) the two DNA bulk-cell samples ($DNA_t_fn and $DNA_g_fn) in a variable $fileList. Then, the samtools and varscan2 are used to call variants of all samples simultaneously.

```sh
snv_fn="output.snp.vcf"
samtools mpileup -f $genomeFasta_b37 $fileList | java -jar varscan2.jar mpileup2snp --min-coverage 5  --min-avg-qual 15 --min-var-freq 0.01 --p-value 1 > $snv_fn
```

The data of variants are collected for the next step to detect cell-specific mutations.

## 5. Cell-specific mutation detection
In this section, we introduce how to use SCmut by an example from a public sc-RNAseq dataset [1]. The data and source codes are available in the SCmut project site (https://github.com/nghiavtr/SCmut).

```R
source("SCmut.R")

# load data
load("example.RData")

### The input data include
# raFull: count from variant allele
# rrFull: count from reference allele
# mut.sites: a list of somatic mutations, possibly discovered from the BC-WES (step 3) or other sources
# cellType: labels of cells (tumor/non-tumor)
###

# germ counts
germ = ncol(raFull)
galt = raFull[,germ]
gn = rrFull[,germ] + raFull[,germ]
germstat = cbind(alt=galt, total=gn)

# observed-mutation sites
ncell = ncol(rrFull)-2
x0.obs = c(rrFull[mut.sites,1:ncell])
   x.obs = x0.obs[!is.na(x0.obs)]
y0.obs = c(raFull[mut.sites,1:ncell])
   y.obs = y0.obs[!is.na(x0.obs)]
nread.obs = x.obs+y.obs
vaf.obs = y.obs/nread.obs

# run fdr2d
set.seed(2017)
fdr = scfdr(rrFull[,1:ncell], raFull[,1:ncell],  mut.sites, germstat)

# plots
par(mar=c(5,5,4,2)+0.1)
contour(fdr$x, fdr$y, fdr$fdr.xy, levels=seq(0.1,1,len=10), xlab='scRNA total reads', ylab='scRNA VAF',cex.axis=2.0, cex.main=2.0, cex.lab=2.0)
# get cell types
tum = cellType$index=='Tumor'
tum.mat = matrix(rep(tum,length(mut.sites)),nrow=length(mut.sites), byrow=TRUE)
tum.indic = c(tum.mat)[!is.na(x0.obs)]
points(nread.obs[tum.indic], vaf.obs[tum.indic], pch=16, col='red',cex=1.0, lwd=2)
points(nread.obs[!tum.indic], vaf.obs[!tum.indic], pch=6, col='blue',cex=0.6, lwd=1)
# plot the cell-specific mutations  
library(gplots)
library(RColorBrewer)
mycol=rev(brewer.pal(n=10, name="RdBu"))
# signif mutated cells: fdr<0.2
signif = c(fdr$ifdr<0.2)[!is.na(x0.obs)]
points(nread.obs[signif], vaf.obs[signif], pch=0, cex=2.0, lwd=2,col=mycol[7])
# signif mutated cells: fdr<0.05
signif = c(fdr$ifdr<0.05)[!is.na(x0.obs)]
points(nread.obs[signif], vaf.obs[signif], pch=0, cex=2.0, lwd=2,col=mycol[10])
```

## 6. License
SCmut uses GNU General Public License GPL-3.

## 7. References
1. Chung, W. et al. Single-cell RNA-seq enables comprehensive tumour and immune cell profiling in primary breast cancer. Nat Commun 8, 15081 (2017).

(update later)
