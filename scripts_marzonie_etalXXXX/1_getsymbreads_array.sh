#!/bin/bash --login
#SBATCH --job-name="bwa-symb0"    
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1     
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=2:00:00		
#SBATCH --account=xxx	
#SBATCH --partition=general	       	
#SBATCH -e bwa-symb0%A_%a.e    
#SBATCH -o bwa-symb0%A_%a.o
#SBATCH --array=1-51

echo "#=================== JOB INFO ===================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "#================================================#"

# define in, out, filtering threshold for coverm filter
aln_id_length=80
INDIR=02_unmappedFQs
OUTDIR=03_symbFAs_${aln_id_length}_${aln_id_length}
mkdir -p ${OUTDIR}

# get noncoral_R1 and noncoral_R2 in array 
FILENAME=`ls ${INDIR}/*_noncoral_1.fastq.gz | sort | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`
name=`basename ${FILENAME} _noncoral_1.fastq.gz`
R1=${INDIR}/${name}_noncoral_1.fastq.gz
R2=${INDIR}/${name}_noncoral_2.fastq.gz
ls $R1 $R2
echo $name
echo $FILENAME

# define symbiodiniaceae reference
# here we are using composite reference that includes
# - Symbiodinium microadriaticum CCMP2467 (Hi-C ver)
# - Breviolum minutum 
# - Cladocopium proliferum SCF055
# - Cladocopium infistilum rt203
# - Cladocopium C15 
# - Cladocopium C92
# - Durusdinium trenchii CCMP2556
# - Effrenium voratium RCC1521
# - Fugacium kawagutii 
REF=Smic1_1N_CCMP2467_v2_Nand2021.Bmin_Chen2020.Cpro_v2_Chen2022.Dtre_CCMP2556_v1_Dougan2022.Evor_RCC1521_v2_Shah.Fkaw_v3_Li2020.Cinf_rt203_GP2024.C15_Robbins2019.C103_v1_1_Chen2020.fa.gz

# module/software required 
module load bwa/0.7.17-gcccore-11.3.0 
module load samtools/1.16.1-gcc-11.3.0

# read mapping 
bwa mem ${REF} -t ${SLURM_CPUS_PER_TASK} ${R1} ${R2} | samtools view -1 -bS - > /scratch/temp/${SLURM_JOB_ID}/${name}.bam

# filter primary alignments (i.e. mapped reads) by alignment identitiy and aligned length, here we using 80% and 80% 
conda activate 70f45da0230f1c4fa9331adb70eda1b9_ 
coverm filter -b /scratch/temp/${SLURM_JOB_ID}/${name}.bam -o /scratch/temp/${SLURM_JOB_ID}/${name}.filtered.bam \
--min-read-percent-aln_id_length $aln_id_length --min-read-aligned-percent $aln_id_length --threads ${SLURM_CPUS_PER_TASK}

# bam to FASTA and gzip
samtools fasta /scratch/temp/${SLURM_JOB_ID}/${name}.filtered.bam -@ ${SLURM_CPUS_PER_TASK} > /scratch/temp/${SLURM_JOB_ID}/${name}.${aln_id_length}-${aln_id_length}filtered.fa
pigz /scratch/temp/${SLURM_JOB_ID}/${name}.${aln_id_length}-${aln_id_length}filtered.fa
cp /scratch/temp/${SLURM_JOB_ID}/${name}.${aln_id_length}-${aln_id_length}filtered.fa.gz ${OUTDIR}


