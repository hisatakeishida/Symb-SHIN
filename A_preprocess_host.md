## Sample information
- We analyzed coral hologenome data of _Acropora kenti_ (Cooke et al., 2020), sampled across the inshore central Great Barrier Reef (GBR)
- These data were generously provided by Dr Ira Cooke and Dr Jia Zhang from James Cook University (JCU)
-  In 2015, 148 adult colonies were sampled from five locations in inshore central GBR
-  These included Magnetic Island, Dunk Island, Fitzroy Island, Pandora Reef, and Pelorus Island (see the map below)
-  Magnetic Island, Dunk Island, and Pandora Reef are known to experience high riverine influence
-  Fitzroy Island and Pelorus Island experience low riverine influence
-  Low-coverage whole genome sequencing (3X of coral genome per sample) was performed using Illumina HiSeq 2500 platform (2*100bp)

## Required softwares 
- FASTQC v0.12.1 (https://github.com/s-andrews/FastQC)
- MultiQC v1.19 (https://github.com/MultiQC/MultiQC)
- Trimmomatic v0.39 (https://github.com/usadellab/Trimmomatic)
- bwa v.0.7.17 (https://github.com/lh3/bwa)
- samtools v.1.19.2 (https://github.com/samtools/samtools)
- bedtools v.2.31.1 (https://github.com/arq5x/bedtools2)
- GATK4
- ANGSD
- PCAngsd
- pigz 

## 1. QC and reprocessing
- Quality check with FASTQC and MultiQC for read quality and adapter contamination
- QC using Trimmomatic (4bp sliding window, minimum phred score quality of 20, minimum read legth of 50 bp, Illuminaclip option in palindrome mode)

## 2. Map post-QC reads to coral genome and get coral and non-coral reads (i.e. reads that did not map against host coral reference genome; avaiable at http://aten.reefgenomics.org/download/) 
```
cd post_qc_raw_reads_merged
# index reference genome 
bwa index aten_final_0.11.fasta.gz
$ process each read pairs
for infile in *_R1.fastq.gz
do
     base=$(basename ${infile} _R1.fastq.gz)
     bwa mem aten_final_0.11.fasta.gz -t 24 ${infile} ${base}_R2.fastq.gz -o 0_coralgenome_mapping/${base}.sam

     # get coral reads using -F4 flag
     samtools view 0_coralgenome_mapping/${base}.sam -F4  -b -@ 96 | samtools sort > coral_bams/${base}.bam -@ 96

     # get non-coral reads using -f12 -F256 flags
     samtools view 0_coralgenome_mapping/${base}.sam -f12 -F256  -b -@ 96 | samtools sort > non_coral_bams/${base}.bam -@ 96

     bedtools bamtofastq -i coral_bams/${base}.bam -fq coral_reads/${base}_1.fq -fq2 coral_reads/${base}_2.fq
     bedtools bamtofastq -i non_coral_bams/${base}.bam -fq noncoral_reads/${base}_1.fq -fq2 noncoral_reads/${base}_2.fq
done

pigz coral_reads/*.fq
pigz noncoral_reads/*.fq

```
## 3. Coral genotyping using GATK, ANGSD, PCANGSD
- Coral genotype can be identified based on genotype likelihood
- See https://ecoevorxiv.org/repository/view/6714/ for more details
- See https://github.com/iracooke/atenuis_wgs_pub for Cooke et al. (2020) and https://github.com/bakeronit/acropora_digitifera_wgs for Zhang et al. (2022) 

![FIGURE 1462  copy](https://github.com/hisatakeishida/Symb-SHIN/assets/95674651/8a91bc78-c762-49e5-9099-c32623fc09f9)






