#!/bin/bash --login
#SBATCH --job-name="magena-d2s_jk"    
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1     
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G			
#SBATCH --time=3:00:00		
#SBATCH --account=xxxx	
#SBATCH --partition=general	       	
#SBATCH -e magena-d2s_jk%A_%a.e
#SBATCH -o magena-d2s_jk%A_%a.o
#SBATCH --array=1-100

# repeat 100 times in array 

echo "#=================== JOB INFO ===================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "#================================================#"

# define portion ${PORTION} to be removed and size ${CHUNK_SIZE} for random sampling  
PORTION=40
CHUNK_SIZE=100

# define in, out, filtering threshold used 
iden_t=80
s_size=40000000
INDIR=05_symbFAseqs_${iden_t}_${iden_t}
OUTDIR=06_jkdistmat_${iden_t}_${iden_t}_${PORTION}_${CHUNK_SIZE}
OUT=${iden_t}_${iden_t}_magena_all.d2s.k21.batch__protion_${PORTION}__chunksize_${CHUNK_SIZE}__s_size_${s_size}__${SLURM_ARRAY_TASK_ID}.txt
mkdir -p $OUTDIR

# need conda env with python 
conda activate python_env
mkdir -p  /scratch/temp/${SLURM_JOB_ID}/batch_${SLURM_ARRAY_TASK_ID}

# jackknife.py from https://github.com/chanlab-genomics/alignment-free-tools/tree/main/jackknife/jackknife.py for jackknifing 
for i in ${INDIR}/*
do
    ls $i
    name=$(basename $i .fa)
    echo $name 
    python 00_jackknife.py --input_paths ${i} --output_path /scratch/temp/${SLURM_JOB_ID}/batch_${SLURM_ARRAY_TASK_ID} --portion $PORTION --chunk_size $CHUNK_SIZE --threads ${SLURM_CPUS_PER_TASK} 
done

conda activate d2ssect
# k-mer enumeration with jellyfish with k=21 
for f in /scratch/temp/${SLURM_JOB_ID}/batch_${SLURM_ARRAY_TASK_ID}/*.fa ;do 
   name=$(basename $f .fa)
   echo $name 
   echo $f

   jellyfish count -t ${SLURM_CPUS_PER_TASK} -m 21 -s "${s_size}" $f -o /scratch/temp/${SLURM_JOB_ID}/batch_${SLURM_ARRAY_TASK_ID}/${name}.jf 
done 

# running d2ssect to generate distance matrix. See https://github.com/bakeronit/d2ssect
d2ssect -l /scratch/temp/${SLURM_JOB_ID}/batch_${SLURM_ARRAY_TASK_ID}/*.jf -f /scratch/temp/${SLURM_JOB_ID}/batch_${SLURM_ARRAY_TASK_ID}/*.fa -o /scratch/temp/${SLURM_JOB_ID}/${OUT} -t ${SLURM_CPUS_PER_TASK}
cp /scratch/temp/${SLURM_JOB_ID}/${OUT} ${OUTDIR}
