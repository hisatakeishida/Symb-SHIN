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
    - SCF_055 (OQ359937), Aten-MI-1 (MW691104.1), Aten-MI-2 (MW691103.1), A03_50 (KF572189.1), Amil-MI (MW691105.1), Aten-WSY (MW691106.1), A02_71 (KF572195.1), JPB08_92_C1_Siderastrea (OQ359935.1), JPB08_77_C1_Siderastrea (OQ359936.1), A03_50_C1_Acropora_tenuis (OQ359938.1), A02_71_C1_Galaxea_fasicularis (OQ359939.1), Zam03_70_C1_Astropora (OQ359940.1)
  - _Cladocopium vulgare_ (OQ359929.1)
  - _Cladocopium latusorum_ (MW819757.1)
  - _Cladocopium pacificum_ (MW861711.1)
  - _Cladocopium madreporum_ (OQ359918.1)
  - _Cladocopium patulum_ (OQ359922.1)
  - _Cladocopium sodalum_
    - A03_75_C3_Coelastrea OQ359899.1), A03_235_C3_Acropora_sarmentosa (OQ359900.1), A03_97_C3_Favites (OQ359901.1), Zan07_314_C3_Acropora (OQ359906.1), Zan07_379_C3_Acropora (OQ359907.1), Zan07_67_C3_Acropora (OQ359908.1), A02_20_C3_Acropora (OQ359902.1), A02_10_C3_Acropora (OQ359903.1), HI07_11_C3_Acropora (OQ359904.1), Pal16_ORT1_40A_C3_Acropora (OQ359905.1), A03_82_C3K_Acropora (OQ359894.1), A03_81_C3K_Acropora (OQ359896.1), A03_327_C3_Echinophyllia_mammiformis (OQ359897.1)
  - _Cladocopium goreaui_
    - RT152 (KF572162.1), RT113 (KF572161.1)
![example (20)-01](https://github.com/hisatakeishida/Symb-SHIN/assets/95674651/0d58ab85-c604-4d19-a5c6-c52784155a87)

## 2. GraftM <a name="graftm"></a>
- Recovering markers in reads using taxonomic sequence identifier (GraftM)
- We used GraftM with custom ITS2 HMM profile (input:1327 ITS2 seqs from SymPortal, removed 121 ITS2 seqs as duplicate, continued analysis with 1210 ITS2 seqs)

```
# creating GraftM with custom ITS2 HMM profile 
graftM create --output ITS2_graftm.gpkg --sequences symportal_ITS2.fa --taxonomy ITS2_taxonomy.txt

#ITS2_taxonomy.txt looks like
C1	k__Dinophyceae,p__Suessiales,c__Symbiodiniaceae,o__Cladocopium,f__C1,g__C1
C1b	k__Dinophyceae,p__Suessiales,c__Symbiodiniaceae,o__Cladocopium,f__C1,g__C1b
C1au	k__Dinophyceae,p__Suessiales,c__Symbiodiniaceae,o__Cladocopium,f__C1,g__C1au

cd noncoral_reads/
for infile in *_1.fq.gz
do
   base=$(basename ${infile} _1.fq.gz)
   graftM graft --forward ${infile} --reverse ${base}_2.fq.gz --graftm_package ITS2_graftm.gpkg --input_sequence_type nucleotide --output_directory graftm_result/${base}
done
```

cd /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/paired_noncoral_reads

# # bwa index /scratch/project/ncmas_d85/ishida/symb_refseq/aligned_trimmed/Cladocopium_psbA_frontier_alignment_trimmed.fa 
# # bwa index /scratch/project/ncmas_d85/ishida/symb_refseq/aligned_trimmed/symportal_ITS2_alignment_trimmed.fa 
bwa index /scratch/project/ncmas_d85/ishida/symb_refseq/psbA-test.fa

for infile in *_1.fq.gz
do
    base=$(basename ${infile} _1.fq.gz)
    bwa mem /scratch/project/ncmas_d85/ishida/symb_refseq/psbA-test.fa -t 24 ${infile} ${base}_2.fq.gz -o /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/psba_sam/${base}.sam
    samtools view /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/psba_sam/${base}.sam -F 4 -b -@ 24 | samtools sort > /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/psba_bam/${base}.bam -@ 24
    samtools index -@ 12 /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/psba_bam/${base}.bam
    rm /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/psba_sam/${base}.sam
    # bwa mem /scratch/project/ncmas_d85/ishida/symb_refseq/aligned_trimmed/symportal_ITS2_alignment_trimmed.fa  -t 24 ${infile} ${base}_2.fq.gz -o /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/ITS2_all_sam/${base}.sam
    # samtools view /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/ITS2_all_sam/${base}.sam -F 4 -b -@ 24 | samtools sort > /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/ITS2_all_bam/${base}.bam -@ 24

    # bwa mem /scratch/project/ncmas_d85/ishida/symb_refseq/aligned_trimmed/symportal_Cladocopium.ITS2_alignment_trimmed.fa  -t 24 ${infile} ${base}_2.fq.gz -o /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/ITS2_C_sam/${base}.sam
    # samtools view /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/ITS2_C_sam/${base}.sam -F 4 -b -@ 24 | samtools sort > /scratch/project/ncmas_d85/ishida/aten_cooke_analysis/marker_alignment/ITS2_C_bam/${base}.bam -@ 24
done









