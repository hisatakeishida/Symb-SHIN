# Read-based (with alignment) analysis 

## Table of Contents
0. [Input data](#input)
1. [Generating consensus marker sequence](#consensus)
2. [GraftM](#graftm)
3. [Alignment to marker/genome/mitogenome](#alignmentgenome)

## 0. Input data <a name="input"></a>
- Non-coral reads for each hologenome sample
- ## Required softwares 
     - Samtools v1.19.2 (https://github.com/samtools/samtools) 
     - bwa v0.7.17 (https://github.com/lh3/bwa)
     - GraftM v0.14.0 (https://github.com/geronimp/graftM)
     - MAFFT v7.471 (https://github.com/GSLBiotech/mafft)
     - seqmagic (https://github.com/fhcrc/seqmagick)
     - MrBayes v3.2.7a (https://nbisweden.github.io/MrBayes/)
     - beagle-lib v4.0.0 (https://github.com/beagle-dev/beagle-lib?tab=readme-ov-file)
     - perl  (https://github.com/Perl/perl5)

## 1. Generating consensus marker sequence <a name="consensus"></a>
- Recoverying consensus marker sequence for each sample 
- Adopting methods from https://github.com/institut-de-genomique/TaraPacific_Pocillopora-transcriptomic (BAMFilteration.pl and MpileuptoConsensus.pl were adopted) 
- We used psbAncr sequences below as a reference (seed) to construct consensus psbAncr sequence of each sample (psbA-seed.fa):
  - _Cladocopium proliferum_ SCF_055_C1_Atenuis_MagneticIsland (OQ359937.1)
  - _Cladocopium vulgare_ A03_85_C1_Acropora (OQ359929.1)
  - _Cladocopium latusorum_ Aust03_41 (MW819757.1)
  - _Cladocopium pacificum_ Phuket07_238 (MW861711.1)
  - _Cladocopium madreporum_ Pal18_ORT1_18_C40_Acropora (OQ359918.1)
  - _Cladocopium patulum_ Zan07_10_C3u_Physogyra (OQ359922.1)
  - _Cladocopium sodalum_ A03_81_C3K_Acropora (OQ359896.1)
  - _Cladocopium goreaui_ RT152 (KF572162.1)

```
bwa index psbA-seed.fa
cd noncoral_reads/
#read mapping
for infile in *_1.fq.gz
do
    base=$(basename ${infile} _1.fq.gz)
    bwa mem psbA-seed.fa ${infile} ${base}_2.fq.gz | samtools view -b -@ 4 -F 4 - -o psba_bam/${base}.bam
    # keep only full-length reads aligned with more than 90% of identity
    perl BAMFilteration.pl -in psba_bam/${base}.bam -out psba_bam/${base}.filtered100-90.bam -minPCaligned 100 -minPCidentity 90    
done

#identify reference with highest number of reads mapped 
for infile in *.filtered100-90.bam
do
    echo -ne "$i\t" ;samtools view $i | awk '{print $3}' | sort | uniq -c | sort -nrk1,1 |head -1 | awk '{print $2"\t"$1}'
done > Summary_mapping_Cladocopium_psbA_100-90.tab

# sort bam files 
cat Summary_mapping_Cladocopium_psbA_100-90.tab | while read a b c ;do samtools sort -o psba_bam/${a%.filtered100-90.bam}.sorted.bam $a; done

# index sorted bam files 
for i in *.sorted.bam
do
    samtools index $i
done

# Generate text pileup output
cat Summary_mapping_Cladocopium_psbA_100-90.tab | while read a b c
do
    samtools mpileup -Aa -f psbA-seed.fa -r $b ${a%.bam}.sorted.bam -o ${a%.bam}.sorted.mpileup
done

# building consensus (minimum coverage of 2) 
for i in *.filtered100-90.sort.mpileup
do
    perl MpileuptoConsensus.pl -mpileup $i -Mincov 2 -out ${i}.consensus.fa
done

# keep consensus sequences in one fasta file 
for i in *.filtered100-90.sort.mpileup.consensus.fa
do
    echo -e ">${i%_*}"; tail -n +2 $i
done > All_psbA-Cladocopium100-90_psbA.aln.mpileup.consensus.fa
```

- Multiple sequence alignment (MSA) of a set of psbA sequences, including consensus psbAncr sequence for each sample and other psbA reference marker sequences (provided below), was generated using MAFFT in mafft-linsi mode (Katoh and Standley 2013)
- Bayesian phylogenetic tree was inferred from the MSA using MrBayes
- Additional markers we used include
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

- MSA using mafft 
```
mafft --maxiterate 1000 --localpair --thread 24 psbA_consenus_and_ref.fa > psbA_consenus_and_ref_alignment.fa

# convert fasta to nexus
seqmagick convert --output-format nexus --alphabet dna psbA_consenus_and_ref_alignment.fa psbA_consenus_and_ref_alignment.nexus
```

- build bayesian phylogeney
```
mpirun -np 1 MrbayesCommands.nex
```

- MrbayesCommands.nex
```
begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute psbA_consenus_and_ref_alignment.nexus;
   lset nst=6 rates=invgamma;
   prset brlenspr=clock:birthdeath;
   mcmc nruns=1 ngen=1000000 samplefreq=100 printfreq=100 diagnfreq=1000 file=psbA_consenus_and_ref_alignment.phylogeny.nexus;
   sump;
   sumt;
end;
```
- Resulting tree was investigated using TVBOT (https://www.chiplot.online/tvbot.html) and visualized in Figure 2B

## 2. GraftM <a name="graftm"></a>
- Recovering markers in reads using taxonomic sequence identifier (GraftM)
- We used GraftM with custom ITS2 HMM profile (input:1327 ITS2 seqs from SymPortal, removed 121 ITS2 seqs as duplicate, continued analysis with 1210 ITS2 seqs)
- Benchmark is required for accuracy (Beta version available at **[ ITS2_graftm_final.gpkg](ITS2_graftm_final.gpkg)**)
- Results are visualised as a heatmap in Figure 2B

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

python graftM_result_summary.py
```

## 3. Alignment to marker/genome/mitogenome <a name="alignmentgenome"></a>
- Use the number of reads uniquely mapped to markers as a proxy of abundance

```
bwa index psbA-Cgor_Cpro_markers.fa

for infile in *_1.fq.gz
do
    base=$(basename ${infile} _1.fq.gz)
    bwa mem psbA-Cgor_Cpro_marker -t 24 ${infile} ${base}_2.fq.gz -o psba_sam_cgor_cpro/${base}.sam
    samtools view psba_sam_cgor_cpro/${base}.sam -F 4 -b -@ 24 | samtools sort > psba_bam_cgor_cpro/${base}.bam -@ 24
    samtools index -@ 12 psba_bam_cgor_cpro/${base}.bam
    rm psba_sam_cgor_cpro/${base}.sam
done
```

- See https://github.com/institut-de-genomique/TaraPacific_Pocillopora-transcriptomic and https://github.com/iracooke/atenuis_wgs_pub for applications of read-based alignment-methods with genome or mitogenome
