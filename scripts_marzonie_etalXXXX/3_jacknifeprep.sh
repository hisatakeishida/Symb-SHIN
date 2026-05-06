#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20G
#SBATCH --job-name=d2s-80_80
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=xxxx
#SBATCH -o d2s-80_80.o
#SBATCH -e d2s-80_80.e

echo "#=================== JOB INFO ===================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "#================================================#"

# define in, out, filtering threshold used for coverm filter
iden_t=80
INDIR=03_symbFAs_${iden_t}_${iden_t}
OUTDIR=05_symbFAseqs_${iden_t}_${iden_t}
threshold=40000000 #needs to be changed according to the smallest datasets you have.
s_size=50000000
OUT=40Mbp_${iden_t}_${iden_t}.d2s.k21.txt

# module/software required 
conda activate seqkit_2.5.1
module load seqtk/1.5-gcc-14.2.0

# downsampling all datasets to $threshdold (i.e. 40Mbp) remove bias from differential coverage across datasets 
for i in ${INDIR}/*.gz
do
    ls $i
    name=$(basename $i .${iden_t}-${iden_t}filtered.fa.gz)
    echo $name 

    len=$(seqkit stats -T $i -j ${SLURM_CPUS_PER_TASK} | tail -n 1  | cut -f 5)

    if [ "$len" -gt "$threshold" ]; then
        ratio=$(echo "scale=5; $threshold / $len" | bc -l)

        formatted_ratio=$(printf "%.5f" $ratio)
        echo "Ratio: $formatted_ratio"

        seqkit sample $i --proportion ${formatted_ratio} --threads ${SLURM_CPUS_PER_TASK} | seqtk seq -a > /scratch/temp/${SLURM_JOB_ID}/${name}.fa
    else
        echo "Length is smaller than threshold, ignoring."
    fi
    # break
done

cp /scratch/temp/${SLURM_JOB_ID}/*.fa $OUTDIR