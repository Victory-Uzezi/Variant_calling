# Variant_calling
CHRONOS PIPELINE FOR VARIANT CALLING AND ANNOTATION

Softwares and the sequences needed for variant calling:

-bwa
-samtools
-bcftools
-vcfutils.pl

We will be working with a reference genome and aligning the reads genome to the reference in order to search if there is a variant in our reads genome.
The genome is a tumor and normal tissue with a reference genome. Genome with tumor tissues are indicated as either "T" or "tumor" and genome with normal tissues are indicated as "N" or "normal" in this pipeline.
The location to check for this variants are Chromosome 5, 12, and 17.

Below are the steps for the variant calling pipeline:

1. create directories for results and reference

code: $ mkdir reference
      $ mkdir results

2. We need to prepare our reads for alignment to our reference:

* the reads consists of normal tissues for forward and reverse strands and a tumor tissues for forward and reverse strands
Reads genome:
https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

* reference genome:
https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

3. Get the reads and reference into the reference directory

i. wget your reference genome as an output fasta.gz using the (-O) and unzip the fasta.gz to a fasta file

$ wget -O reference/reference.fasta.gz hg19.fa.gz
$ gunzip reference/reference.fasta.gz
	=reference.fasta

example:
$ wget -O reference/hg19_tissue.fasta.gz https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
$ gunzip reference/hg19_tissue.fasta.gz

ii. wget your reads into the reference directory as a single fastq file and unzip. 
note: assuming your reads are normal_1.fastq.gz normal_2.fastq.gz tumor_1.fastq.gz tumor_2.fastq.gz and they are consisted of both forward and reverse strands.
 
code:
$ wget normal_1.fastq.gz normal_2.fastq.gz tumor_1.fastq.gz tumor_2.fastq.gz
$ gunzip *fastq.gz

##make directory in the reference directory and send your reads into it
code:
$ mkdir reads.fastq
$ mv normal_1.fastq normal_2.fastq tumor_1.fastq tumor_2.fastq reads.fastq

example:
wget reference/https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
gunzip *fastq.gz
mkdir reads.fastq
mv SLGFSK-N_231335_r1_chr5_12_17.fastq SLGFSK-N_231335_r2_chr5_12_17.fastq SLGFSK-T_231336_r1_chr5_12_17.fastq SLGFSK-T_231336_r2_chr5_12_17.fastq reads.fastq/

4. make directories in the results directory for our bam, sam, bcf, and vcf outputs using one line of code (-p).

code:
$ mkdir -p results/sam results/bam results/bcf results/vcf

5. the next step is to index our reference genome using BWA (burrow-wheeler aligner). This allow it to more efficiently search the genome during sequence alignment

code:
$ bwa index reference/reference.fasta

example:
bwa index reference/hg19_tissue.fasta

while the indexing is created, you will see it unpacking as thus:




6. the next is to align the reads to the reference genome using the BWA MEM command and send the result as a sam file to the results directory. This alignment of reads to reference tells you which precise location in the genome each base pair in each sequencing read comes from

code:
$ cd reference
$ bwa mem reference.fasta reads_fastq/normal_1.fastq reads_fastq/normal_2.fastq >../results/sam/normal_aligned.sam
$ bwa mem reference.fasta reads_fastq/tumor_1.fastq reads_fastq/tumor_2.fastq >../results/sam/tumor_aligned.sam

example:
bwa mem hg19_tissue.fasta reads.fastq/SLGFSK-N_231335_r1_chr5_12_17.fastq reads.fastq/SLGFSK-N_231335_r2_chr5_12_17.fastq > ../results/sam/SLGFSK-N_normal_aligned.sam
bwa mem hg19_tissue.fasta reads.fastq/SLGFSK-T_231336_r1_chr5_12_17.fastq reads.fastq/SLGFSK-T_231336_r2_chr5_12_17.fastq > ../results/sam/SLGFSK-Tumor_aligned.sam


tips: a variety of new alignment tools (Langmead et al., 2009; Li et al., 2008) have been designed to realize efficient read mapping against large reference sequences, including the human genome.
 A common alignment format that supports all sequence types and aligners creates a well-defined interface between alignment and downstream analyses, including variant detection, genotyping and assembly.
The Sequence Alignment/Map (SAM) format is designed to achieve this goal. The SAM format consists of one header section and one alignment section. The lines in the header section start with character ‘@’, and lines in the alignment section do not. All lines are TAB delimited.

In SAM, each alignment line has 11 mandatory fields and a variable number of optional fields. The optional fields are presented as key-value pairs in the format of TAG:TYPE:VALUE. They store extra information from the platform or aligner. 
To improve the performance, we designed a companion format Binary Alignment/Map (BAM), which is the binary representation of SAM and keeps exactly the same information as SAM. 

7. Convert sam to bam format using samtools with the view command and flag S (-S) indicating that the input file is in sam format and to output it as bam format using flag b (-b)
Note: The samtools view command is the most versatile tool in the samtools package. It’s main function, not surprisingly, is to allow you to convert the binary (i.e., easy for the computer to read and process) alignments in the BAM file view to text-based SAM alignments that are easy for humans to read and process.
code:
$ samtools view -S -b results/sam/normal_aligned.sam > bam/normal_aligned.bam
$ samtools view -S -b results/sam/tumor_aligned.sam > bam/tumor_aligned.bam

Example:
samtools view -S -b sam/SLGFSK_normal_aligned.sam > bam/SLGFSK_normal_aligned.bam
samtools view -S -b sam/SLGFSK-Tumor_aligned.sam > bam/SLGFSK-Tumor_aligned.bam


8. Sort your bam files. This is necessary because doing anything meaningful such as calling variants or visualizing alignments in IGV) requires that the BAM is further manipulated. It must be sorted such that the alignments occur in “genome order”. That is, ordered positionally based upon their alignment coordinates on each chromosome.
code:
$ samtools sort bam/normal_aligned.bam -o bam/normal_sorted.bam
$ samtools sort bam/tumor_aligned.bam -o bam/normal_sorted.bam

example:
samtools sort bam/SLGFSK_normal_aligned.bam -o bam/SLGFSK_normal_sorted.bam
samtools sort bam/SLGFSK-Tumor_aligned.bam -o bam/SLGFSK-Tumor_sorted.bam


9. View the first 5 of your sorted bam files
code:
$ samtools view bam/normal_sorted.bam | head -n 5
$ samtools view bam/tumor_sorted.bam | head -n 5

example:
samtools view bam/SLGFSK_normal_sorted.bam | head -n 5
samtools view bam/SLGFSK-Tumor_sorted.bam | head -n 5

10. you can get the statistics of your sorted bam files using:
code:
$ samtools flagstat bam/normal_sorted.bam
$ samtools flagstat bam/tumor_sorted.bam

example:
samtools flagstat bam/SLGFSK_normal_sorted.bam


VARIANT CALLING
Variant calling entails identifying single nucleotide polymorphisms (SNPs) and small insertions and deletion (indels) from next generation sequencing data.

insert pictures for indicating reads to variant calling

step 1: count the read coverage with bcftools and use the command mpileup, 
-O b to generate a bcf format 
-o to output it
-f path to the reference genome and the sorted bam file.
 
code:
$ bcftools mpileup -O b -o bcf/normal_raw.bcf -f ../ reference/hg19_tissue.fasta bam/normal_sorted.bam
$ bcftools mpileup -O b -o bcf/tumor_raw.bcf -f ../reference/hg19_tissue.fasta bam/tumor_sorted.bam

alternatively, you can achieve same with the following commands:
$ samtools mpileup -g -f /reference/hg19_tissue.fasta bam/bam/normal_sorted.bam > bcf/normal_raw.bcf
$ samtools mpileup -g -f /reference/hg19_tissue.fasta bam/bam/tumor_sorted.bam > bcf/tumor_raw.bcf

-g generate genotype likelihoods in bcf format
-f fasta reference file

example:
bcftools mpileup -O b -o results/bcf/SLGFSK_normal_raw.bcf -f ../reference/hg19_tissue.fasta bam/SLGFSK_normal_sorted.bam
bcftools mpileup -O b -o results/bcf/SLGFSK-Tumor_raw.bcf -f ../reference/hg19_tissue.fasta bam/SLGFSK-Tumor_sorted.bam


step 2: detect the single nucleotide variants (SNVs) using 
-bcftools
-call 
-O b output type in bcf
-v for indicating the output of only the variant sites, 
-c use the original calling method

code:
$ bcftools call -O b -vc results/bcf/SLGFSK_normal_raw.bcf > bcf/SLGFSK_normal_var.bcf
$ bcftools call -O b -vc results/bcf/SLGFSK_tumor_raw.bcf > bcf/SLGFSK_tumor_var.bcf

example:
bcftools call -O b -vc results/bcf/SLGFSK_normal_raw.bcf > bcf/SLGFSK_normal_var.bcf 
bcftools call -O b -vc results/bcf/SLGFSK-Tumor_raw.bcf > bcf/SLGFSK-Tumor_var.bcf


step 3: filter SNVs with vcfutils.pl
code:
$ bcftools view results/bcf/SLGFSK_normal_var.bcf | vcfutils.pl varFilter - > vcf/SLGFSK_normal_final.vcf
$ bcftools view results/bcf/SLGFSK-Tumor_var.bcf | vcfutils.pl varFilter - > vcf/SLGFSK-Tumor_final.vcf

example:
bcftools view results/bcf/SLGFSK_normal_var.bcf | vcfutils.pl varFilter - > vcf/SLGFSK_normal_final.vcf
bcftools view results/bcf/SLGFSK-Tumor_var.bcf | vcfutils.pl varFilter - > vcf/SLGFSK-Tumor_final.vcf

note: vcf file is made up of;
1. meta-information which starts the vcf and it is indicated by ##
2. header ( CHROM, POS, ID, REF, QUAL, FILTER etc)
3. variants 

example:
less -S vcf/SLGFSK_normal_final.vcf


