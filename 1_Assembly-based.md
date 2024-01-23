# Perform metagenomic assembly to recover contigs or MAGs from hologenome dataset

## Table of Contents
0. [Input data](#input)
1. [Metagenomic assembly](#MA)
2. [Recover markers in contigs](#generecov)
3. [Recover MAGs](#magrecov)

## 0. Input data <a name="input"></a>
- Non-coral reads for each hologenome sample
- ## Required softwares 
     - metaSPAdes
     - Blastn

## 1. Metagenomic assembly <a name="MA"></a>
- Perfomring metagenomic assembly (metaSPAdes) of non-coral reads from each hologenome sample

```
cd noncoral_reads/
for infile in *_1.fq.gz
do
     base=$(basename ${infile} _1.fq.gz)
     spades.py -1 ${infile} -2 ${base}_2.fq.gz --memory 1000 --meta -t 64 -o noncoral_assembly_tmp/${base}
     cp noncoral_assembly_tmp/${base}/contigs.fasta noncoral_assembly/${base}_contigs.fasta
     rm -r noncoral_assembly_tmp/${base}
done
```

## 2. Recover markers in contigs <a name="generecov"></a>
- We used Blastn to identify contigs that contain reference markers of interests
- We searched for contigs that contain full-length mitochondrial cytochrome b (mtCOB) sequence of _Cladocopium goreaui_ RT152 (KF206028.1)

```
cd noncoral_assembly/
for infile in *_contigs.fasta
do
     base=$(basename ${infile} _contigs.fasta)
     makeblastdb -in ${infile} -dbtype nucl -parse_seqids 
     blastn -query C_gor_RT152_mtcob.fa -db ${infile} -outfmt 6 -perc_identity 100 -evalue 1e-50 -out mtcob_blast_hit/${base}
     sort -k3,3g -k12,12gr -k11,11g  mtcob_blast_hit/${base} | sort -u -k1,1 --merge > mtcob_blast_hit/${base}.sorted
done

# Put top hits from all samples into one text file
python blast_sum.py

- 39 out of 64 hologenome assesmbly had contigs that included full-length (914 bp) mitochondrial cytochrome b (mtCOB) sequence of _Cladocopium goreaui_ RT152 (KF206028.1)
```

## 3. Recover MAGs <a name="magrecov"></a>
- This is a work in progress
- MAGs of symbionts can be recovered, but requires additional pre-filtering to minimize the impact of contamination to assembly
- Samples with similar community composition of microbial taxa of interests can be identified (based on all the approaches we discuss in this paper), prior to assembly 
- This will facilitate strain-aware assembly
- Here, we provided an workflow to recover MAGs of important baterial symbionts of coral, _Endozoicomonas_ sp, as an example for such analysis
- Required softwares 
     - Aviary
     - d2ssect
### 3.1 Sample clustering 
- Check **[2.1_Read-based_alignment.md](2.1_Read-based_alignment.md)** and **[2.2_Read-based_alignment-free.md)** for different methods to identify corals with similar microbial composition
  
### 3.2 Metagenomic assembly and recovery 

### 3.3 Alignment-free phylogeny  











