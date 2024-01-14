# Read-based (with alignment) analysis 

## Table of Contents
0. [Input data](#input)
1. [Generating consensus marker sequence](#consensus)
2. [GraftM](#graftm)
3. [Alignment to marker/genome/mitogenome](#alignmentgenome)

## 0. Input data <a name="input"></a>
- Non-coral reads for each hologenome sample
- Softwares (Samtools, ) 

## 1. Generating consensus marker sequence <a name="consensus"></a>
- Recoverying consensus marker sequence for each sample
- See https://github.com/institut-de-genomique/TaraPacific_Pocillopora-transcriptomic for workflow details 
- We used psbAncr sequences below as a reference to construct consensus psbAncr sequence of each sample:
  - _Cladocopium proliferum_ SCF_055_C1_Atenuis_MagneticIsland (OQ359937.1)
  - _Cladocopium vulgare_ A03_85_C1_Acropora (OQ359929.1)
  - _Cladocopium latusorum_ Aust03_41 (MW819757.1)
  - _Cladocopium pacificum_ Phuket07_238 (MW861711.1)
  - _Cladocopium madreporum_ Pal18_ORT1_18_C40_Acropora (OQ359918.1)
  - _Cladocopium patulum_ Zan07_10_C3u_Physogyra (OQ359922.1)
  - _Cladocopium sodalum_ A03_81_C3K_Acropora (OQ359896.1)
  - _Cladocopium goreaui_ RT152 (KF572162.1)

- Multiple sequence alignment (MSA) of psbA sequence set, including consensus psbAncr sequence for each sample and other psbA reference marker sequences (provided below), was generated using MAFFT v7.471 in mafft-linsi mode (Katoh and Standley 2013)
- MSA was trimmed using trimAl v1.4.rev15 with-automated1 (Capella-Gutiérrez et al. 2009)
- Maximum-likelihood phylogenetic tree was inferred from the trimmed MSA using IQ-TREE v2.1.3 (Nguyen et al. 2015)
- Additional markers we compared include
  - _Cladocopium proliferum_
    - SCF_055 (OQ359937)
    - Aten-MI-1 (MW691104.1)
  - _Cladocopium vulgare_ (OQ359929.1)
  - Cladocopium.latusorum MW819757.1
  - Cladocopium.pacificum MW861711.1
  - Cladocopium.madreporum OQ359918.1
  - Cladocopium.patulum OQ359922.1
  - Cladocopium.sodalum OQ359896.1

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
```

## 3. Recover MAGs <a name="magrecov"></a>
- This is a work in progress
- MAGs of symbionts can be recovered, but requires additional pre-filtering to minimize the impact of contamination to assembly
- Samples with similar community composition of microbial taxa of interests can be identified (based on all the approaches we discuss in this paper), prior to assembly 
- This will facilitate computationally intensive strain-aware assembly
- Here, we provided an workflow to recover MAGs of important baterial symbionts of coral, _Endozoicomonas_ sp, as an example for such analysis

### 3.1 Sample clustering 
### 3.2 Metagenomic assembly and recovery 
### 3.3 Alignment-free phylogeny 















