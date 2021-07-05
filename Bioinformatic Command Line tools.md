# Bioinformatic Command Line tools

## FASTQ

### Generate index and dictionary

```bash
#Generate Index
samtools faidx Homo_sapiens.GRCh38.fa
# Generate Dict
picard CreateSequenceDictionary -R Homo_sapiens.GRCh38.fa
```

## BAM files

mote at http://www.htslib.org/workflow/

### view BAM file

```bash
samtools view -b file.bam | less
```

### generate index

```bahs
samtools index *.bam --threads 8
```

### Check Read Group

```bash
samtools view -H file1.bam | grep '^@RG'
samtools stats --split RG <file1.bam>
```

### Split Bam By Read Group

Read documentation for info to include in the read groups

```bash
parallel samtools split $FILE_DIR/input.bam -f $OUT_DIR/"%*_%#.bam"
```

### filter only specified UMI from BAM file (scRNAseq)

`LC_ALL=C` and `-F` (fixed search) are required for speed improvement. Do not add `-i` (ignore case) because it will significantely affect performances

```bash
# Code adapted from https://kb.10xgenomics.com/hc/en-us/articles/360022448251-Is-there-way-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
for i in *.splitNtrim.bam
do
# Save the header lines
	samtools view -H $i > ${i//splitNtrim.bam/filtered.sam}
  # Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
 	samtools view $i | LC_ALL=C grep -F -f $WD/processed/filtered_umi/${i//splitNtrim.bam/barcodes.txt} >> ${i//splitNtrim.bam/filtered.sam}
	# Convert filtered.sam to BAM format
 	samtools view -b ${i//splitNtrim.bam/filtered.sam} > $OUT_DIR/${i//splitNtrim.bam/filtered.bam}
 	rm -f ${i//splitNtrim.bam/filtered.sam}
done
```

**Parallelized version**

```bash
find ./ -name "*.splitNtrim.bam" | sed "s/.*\///" | \
parallel "samtools view -H {} > {=s/splitNtrim.bam/filtered.sam/=} && \
		  samtools view {} | LC_ALL=C grep -F -f $WD/processed/filtered_umi/{=s/splitNtrim.bam/barcodes.txt/=} >> {=s/splitNtrim.bam/filtered.sam/=} && \
		  samtools view -b {=s/splitNtrim.bam/filtered.sam/=} > $OUT_DIR/{=s/splitNtrim.bam/filtered.bam/=} && \
		  rm -f {=s/splitNtrim.bam/filtered.sam/=}"
```

the **filter.txt** should be like that:

```
CB:Z:AAACCTGCAACACGCC
CB:Z:AAACCTGTCGTTTGCC
CB:Z:AAAGATGCACGCTTTC
...
```

### Convert SAM to BAM, filter by -Flag, sort and generate index

```bash
# Convert to binary filtering for flag 3844: read unmapped, not primary alignment,
# read fails platform/vendor quality checks, read is PCR or optical duplicate, supplementary alignment
find ./processed/sam -name "*.sam" | parallel 'samtools view -F 3844 -S -b "{}" > ./processed/bam/"{/.}".bam' &&
# Sort by genome position
find ./processed/bam -name "*.bam" | parallel 'samtools sort "{}" -o "{.}".sorted.bam' &&
# Remove unsorted files
#rm ./processed/bam/*R1.bam
# Generate index file
find ./processed/bam -name "*.sorted.bam" | parallel samtools index {}

```

### Convert BAM to CRAM

```bash
samtools view -T GRCh38.fasta -C -o file.cram file.bam
```

### Convert BAM to BED and add offset (for ATACseq footprint analysis)

```shell
bedtools bamtobed -i bowtie_dup_rm.bam > my.bed
# Shift the forward reads 4bp and reverse reads 5bp:
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' my.bed > my_shifted.bed

```

### generate positive and negative coverage visualization from BAM file

```bash
samtools view -h -b reads.sorted.bam | \
bedtools bamtobed -ed -i - | \
bedtools genomecov -split -strand "+" -i - -g hg19.chromosomeSize -bg -trackline -trackopts 'name="File_name" visibility=2 color=255,30,30' > reads.bedgraph

samtools view -h -b reads.sorted.bam | \
bedtools bamtobed -ed -i - | \
bedtools genomecov -split -strand "-" -i - -g hg19.chromosomeSize -bg | awk '{OFS="\t"}{print $1,$2,$3,$4*-1;}' >> reads.bedgraph

# You need to resort the file since we simply uppended the positive and negative reads one belove the other
igvtools sort reads.bedgraph reads.sorted.bedgraph

# Compress to binary version to increase loading speed (also applies some smoothing)
igvtools toTDF reads.sorted.bedgraph reads.tdf hg19
```

## Cut adapters: cutadapt tool

```bash
# Code for illumina TruSeq Bulk RNAseq processing

# cutadapt can't work in multicore on python 3.8 on macOS so you need a conda enviroment with python 3.7 and cutadapt installed
conda activate cutadapt-env

# Adapter sequence obtained from  https://gist.github.com/photocyte/3edd9401d0b13476e60f8b104c2575f8
# Cut TruSeq Adapter, Index 5 from R1
# Cut Illumina Single End PCR Primer 1 (reverse complemented) from R2
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  --minimum-length 20 \
  --cores 0 \
  --output $WD/variant_call/processed/trimmed_fastq/SW480_bulk_ACAGTG_L007_R1_001.trimmed.fastq.gz \
  --paired-output $WD/variant_call/processed/trimmed_fastq/SW480_bulk_ACAGTG_L007_R2_001.trimmed.fastq.gz \
  $FILE_DIR/SW480_bulk_ACAGTG_L007_R1_001.fastq.gz \
  $FILE_DIR/SW480_bulk_ACAGTG_L007_R2_001.fastq.gz |& tee $WD/variant_call/logs/cutadapt.log
```

## STAR

### Genome build

https://groups.google.com/g/rna-star/c/h9oh10UlvhI/m/BfSPGivUHmsJ

> James is right, using large enough --sjdbOverhang is safer and should not generally cause any problems with reads of varying length.
>
> If your reads are very short, <50b, then I would strongly recommend using optimum --sjdbOverhang=mateLength-1
>
> By mate length I mean the length of one of the ends of the read, i.e. it's 100 for 2x100b PE or 1x100b SE.
>
> For longer reads you can simply use generic --sjdbOverhang 100.
>
> It is a bit confusing because of the way I named this parameter. --sjdbOverhang <Noverhang> is only used at the genome generation step  for constructing the reference sequence out of the annotations.
>
> Basically, the Noverhang exonic bases from the donor site and Noverhang exonic bases from the acceptor site are spliced together for each of the junctions, and these spliced sequences are added to the genome sequence.
>
> At the mapping stage, the reads are aligned to both genomic and splice sequences simultaneously. If a read maps to one of spliced sequences and crosses the "junction" in the middle of it, the coordinates of two pspliced pieces are translated back to genomic space and added to the collection of mapped pieces, which are then all "stitched" together to form the final alignment. Since in the process of "maximal mapped length" search the read is split into pieces of no longer than --seedSearchStartLmax (=50 by default) bases, even if the read (mate) is longer than --sjdbOverhang, it can still be mapped to the spliced reference, as long as --sjdbOverhang > --seedSearchStartLmax.

### Star alignment

Requires to be executed as a bash script to make the <(gunzip -c file.fastq.gz) trick work

```bash
#!/bin/bash
STAR --runThreadN 8 --genomeDir $GENOME_DIR \
  --readFilesIn <(gunzip -c $WD/file1_R1.fastq.gz) <(gunzip -c $WD/file1_R2.fastq.gz) \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix $WD/output/bam/file1
```

#### Two pass alignment 

Two pass is recomended for improved quantification of novel splice mappings

1. An easy way is adding ```--twopassMode Basic``` as option.
2. The ideal and manual way is to run a first single pass alignment for all the sample, collect the *SJ.out.tab* files from all samples and use pass all of them for a 2nd single pass alignment with the option ```--sjdbFileChrStartEnd /path/to/sj1.tab /path/to/sj2.tab ...``` 



## Run programs that requires different version of Java

https://gatk.broadinstitute.org/hc/en-us/articles/360035891211-Need-to-run-programs-that-require-different-versions-of-Java

Allways install JDK version and not base because is more complete and some tools requires it

```bash
/usr/libexec/java_home -v 1.7.0_79 --exec java -jar GenomeAnalysisTK.jar -T ...
```

of course, this is not supporting the standard GATK script call

## VCF

#### bgzip and generate index

```bash
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz
```

### Remove all reads with filter different from PASS 

```bash
bcftools view -f PASS sample_filtered.bam > sample_PASS.bam
# or with VCF tools
# use -c or --stdout to redirect the output to stdoutput
vcftools --vcf sample_filtered.bam -c --remove-filtered-all > sample_PASS.bam
# you can use out if you want to simply specify a prefix
vcftools --vcf sample_filtered.bam --out pass_ --remove-filtered-all
```

### Split multiallelic variants into multiline biallelic variants

From: https://www.biostars.org/p/270693/

* `-Ov` specifies the **O**utput to be *uncompressed **v**cf*

* `-m-any` or `--multiallelics -any` specifies to reduce (-) multiallelig variants into biallelic for both SNPs and InDels (note the space in the second version)

```bash
bcftools norm -Ov --multiallelics-any MyFile.vcf > MyFile.Split.vcf
```

 Also get into the habit of checking your ref alleles against an existing reference genome (this also left-aligns indels):

```bash
bcftools norm --multiallelics-any MyFile.vcf | bcftools norm -Ov --check-ref w -f hg38.fa -o MyFile.Split.RefCheck.vcf
```

*Note: the normalization needs to be applied after splitting multiallelic variants because is not possible to simoultaniously normalized an SNV and an Indel:*

```
REF		ALT
A			C,*
TA		TC,T
```

