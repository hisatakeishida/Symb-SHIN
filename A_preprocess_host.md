## Sample information
- We analyzed coral hologenome data of _Acropora kenti_ (Cooke et al., 2020), sampled across the inshore central Great Barrier Reef (GBR)
- These data were generously provided by Dr Ira Cooke and Dr Jia Zhang from James Cook University (JCU)
-  In 2015, 148 adult colonies were sampled from five locations in inshore central GBR
-  These included Magnetic Island, Dunk Island, Fitzroy Island, Pandora Reef, and Pelorus Island (see the map below)
-  Magnetic Island, Dunk Island, and Pandora Reef are known to experience high riverine influence
-  Fitzroy Island and Pelorus Island experience low riverine influence
-  Low-coverage whole genome sequencing (3X of coral genome per sample) was performed using Illumina HiSeq 2500 platform (2*100bp)

## Required softwares 
- bwa
- samtools
- bedtools
- GATK
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
- See https://github.com/iracooke/atenuis_wgs_pub for Cooke et al. (2020) and https://github.com/bakeronit/acropora_digitifera_wgs for Zhang et al. (2022) 

## 4. Get putative Symbiodiniaceae reads by mapping non-coral reads to Symbiodiniaceae genomes including 
- Clade A: Symbiodinium microadriaticum CCMP2467 (Nand et al., 2021)
- Clade B: Breviolum minutum Mf1.05b (Shoguchi et al., 2013)
- Clade C: Cladocopium goreaui SCF055-01 (Chen et al., 2022)
- Clade D: Durusdinium trenchii CCMP2556 (Dougan et al., 2023)
- Clade E: Effernium voratum RCC1521 (Shah et al., 2023a)
- Clade F: Fugacium kawagutii CCMP2468 (Li et al., 2019)

```
cd noncoral_reads/
for infile in *_1.fq.gz
do
     base=$(basename ${infile} _1.fq.gz)
     bwa mem Smicro_CCMP2467_v2_Nand_et_al_2021_HiC_Smic1_1N.Bmin_Chen_et_al_2020.Cgor_v2_Chen_et_al_2022.Dtre_CCMP2556_v1_Dougan_et_al_2022.Evor_RCC1521_v2_Shah_et_al_XXXX.FkawHic_v3_Li_et_al_2020.fa.gz -t 24 ${infile} ${base}_2.fq.gz -o 0_symbgenome_mapping/${base}.sam
     samtools view 0_symbgenome_mapping/${base}.sam -F 4 -b -@ 24 | samtools sort > symb_map_bams/${base}.bam -@ 24
     samtools fasta symb_map_bams/${base}.bam -@ 24 > symb_reads_fa/${base}.fa
done
```
![FIGURE 1462  copy](https://github.com/hisatakeishida/Symb-SHIN/assets/95674651/8a91bc78-c762-49e5-9099-c32623fc09f9)






