#!/bin/bash --login
#SBATCH --job-name="bwa-noncoral"    
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1     
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --time=12:00:00		
#SBATCH --account=xxx	
#SBATCH --partition=general	       	
#SBATCH -e bwa-noncoral%A_%a.e    
#SBATCH -o bwa-noncoral%A_%a.o
#SBATCH --array=1-51

echo "#=================== JOB INFO ===================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "#================================================#"

# deine in and out
INDIR=01_rawFQs
OUTDIR=02_unmappedFQs
mkdir -p ${OUTDIR}

# get reads in array 
FILENAME=`ls ${INDIR}/*_1.fastq.gz | sort | tail -n +${SLURM_ARRAY_TASK_ID} | head -1`
name=`basename ${FILENAME} _1.fastq.gz `
R1=${INDIR}/${name}_1.fastq.gz 
R2=${INDIR}/${name}_2.fastq.gz 
echo $FILENAME
echo $name
ls $R1 $R2

# define reference 
REF=/scratch/project_mnt/S0026/ishida/symb/Magena/00_genref/Magena_fourcorals.genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/634/125/GCA_014634125.1_Agem_1.0/GCA_014634125.1_Agem_1.0_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/669/915/GCF_036669915.1_ASM3666991v2/GCF_036669915.1_ASM3666991v2_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/949/126/865/GCA_949126865.1_jaMonCapi2.1/GCA_949126865.1_jaMonCapi2.1_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/212/065/GCA_964212065.1_jaIsoPali11.1/GCA_964212065.1_jaIsoPali11.1_genomic.fna.gz 

# bwa index $REF if you haven't

# module/software required 
module load bwa/0.7.17-gcccore-11.3.0 
module load samtools/1.16.1-gcc-11.3.0
module load bedtools/2.30.0-gcc-11.3.0

# read mapping 
bwa mem ${REF} -t ${SLURM_CPUS_PER_TASK} ${R1} ${R2} -o /scratch/temp/${SLURM_JOB_ID}/${name}.sam

# remove reads that mapped (i.e. retaining only the non-coral fraction)
samtools view /scratch/temp/${SLURM_JOB_ID}/${name}.sam -f 12 -F 256 -b -@ ${SLURM_CPUS_PER_TASK} | samtools sort > /scratch/temp/${SLURM_JOB_ID}/${name}.bam -@ ${SLURM_CPUS_PER_TASK}

# bamfiles to paired R1 and R2
bedtools bamtofastq -i /scratch/temp/${SLURM_JOB_ID}/${name}.bam -fq /scratch/temp/${SLURM_JOB_ID}/${name}_noncoral_1.fastq -fq2 /scratch/temp/${SLURM_JOB_ID}/${name}_noncoral_2.fastq

# gzip 
pigz /scratch/temp/${SLURM_JOB_ID}/${name}_noncoral_1.fastq
pigz /scratch/temp/${SLURM_JOB_ID}/${name}_noncoral_2.fastq
cp /scratch/temp/${SLURM_JOB_ID}/${name}_noncoral_1.fastq.gz ${OUTDIR}
cp /scratch/temp/${SLURM_JOB_ID}/${name}_noncoral_2.fastq.gz ${OUTDIR}

