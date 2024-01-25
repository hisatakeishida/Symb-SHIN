![image](https://github.com/hisatakeishida/Symb-SHIN/assets/95674651/a225b2c2-050a-4d34-8acf-c5ac66fe221d)# Perform read-based k-mer based alignment-free analysis on hologenome dataset

## Table of Contents
0. [Input data](#input)
1. [k-mer based alignment-free sample distance calculation](#d2s)

## 0. Input data <a name="input"></a>
- Non-coral reads for each hologenome sample
- Required softwares 
     - bwa v.0.7.17 (https://github.com/lh3/bwa)
     - samtools v.1.19.2 (https://github.com/samtools/samtools)
     - d2ssect (https://github.com/bakeronit/d2ssect)
 
## 1. k-mer based alignment-free sample distance calculation <a name="d2s"></a>
- Map non-coral reads to 6 genomes of Symbiodiniaceae (each representing their genera) and ITS2 from SymPortal database to recover putative Symbiodiniaceae reads
- Genomes we used include
  - _Symbiodinium_: _Symbiodinium microadriaticum_ CCMP2467 (Nand et al., 2021)
  - _Breviolum_: _Breviolum minutum_ Mf1.05b (Shoguchi et al., 2013) 
  - _Cladocopium_: _Cladocopium proliferum_ SCF055-01 (Chen et al., 2022) 
  - _Durusdinium_: _Durusdinium trenchii_ CCMP2556 (Dougan et al ., 2023)
  - _Effrenium_: _Effrenium voratum_ RCC1521 (Shah et al., 2023a)
  - _Fugacium_: _Fugacium kawagutii_ CCMP2468 (Li et al., 2019)

```
bwa index Symb_ref_genome.ITS2.fa

cd noncoral_reads/

for infile in *_1.fq.gz
do
    base=$(basename ${infile} _1.fq.gz)
    bwa mem Symb_ref_genome.ITS2.fa ${infile} ${base}_2.fq.gz -o symb_mapping/${base}.sam    
    samtools view symb_mapping/${base}.sam -F 4 -b | samtools sort > symb_mapping/${base}.bam 
    samtools fasta symb_mapping/${base}.bam > symb_mapping/${base}.fa
done
```

- Quantify metagenomic variation of Symbiodiniaceaen communities across samples using putative extracted Symbiodiniaceaen sequences 
- d2ssect calculates an alignment-free distance between samples based on D2S statistic

```
cd  symb_mapping/

for f in *.fa ;do
     jellyfish count -m 21 -s 10000000 $f -o ${f%.fa}.jf
done

d2ssect -l *.jf -f *.fa -o symb_community_distance.txt 
```
- d2ssect generates D2s-derived distance matrix, representing metagenomic variation of Symbiodiniaceaen communities across samples, and this can be analysed using different dimensionality reduction methods such as Non-metric multidimensional scaling (NMDS) and Distance-based redundancy analysis (db-RDA) or tree-building (e.g. hierarchical clustering and Neighbour joining) 
- Used in Figure 2C
- See https://github.com/bakeronit/acropora_digitifera_wgs for application of this method
