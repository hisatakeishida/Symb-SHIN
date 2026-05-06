#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --job-name=d2s-hclust
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=xxxx
#SBATCH -o d2s-hclust.o
#SBATCH -e d2s-hclust.e

echo "#=================== JOB INFO ===================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "#================================================#"

# define in, out, filtering thresholds
iden_t=80
INDIR=06_jkdistmat_${iden_t}_${iden_t}
OUTDIR=07_HCLUST_${iden_t}_${iden_t}
cp ${INDIR}/* /scratch/temp/${SLURM_JOB_ID}

ref_matrix=40Mbp_80_80_magena_all.d2s.k21.txt #reference distance matrix from full dataset 
ref_hclust=/scratch/temp/${SLURM_JOB_ID}/ref_hclust.awk
outfilename=${iden_t}_${iden_t}_hclust_pseudoreplicates_ref_bs.unroot.awk
outfilename_avg=${iden_t}_${iden_t}_hclust_pseudoreplicates_avg_bs.unroot.awk
mkdir -p /scratch/temp/${SLURM_JOB_ID}/hclust_replicates

# module/software required 
module load r/4.4.0-gfbf-2023a

# generate original dendrogram 
Rscript 01_hclust.R ${ref_matrix} ${ref_hclust}

# generate dendrogram for each pseudo-replicate 
for i in /scratch/temp/${SLURM_JOB_ID}/*.txt
do
    ls $i 
    name=$(basename $i .txt)
    echo $name 
    outpath=/scratch/temp/${SLURM_JOB_ID}/hclust_replicates/${name}.awk

    echo $outpath
    Rscript hclust.R ${i} ${outpath} 
done

# modified from https://github.com/chanlab-genomics/alignment-free-tools/blob/main/jackknife/Jackknife.r 
Rscript 01_Jackknife_param.r \
${ref_hclust} \
/scratch/temp/${SLURM_JOB_ID}/hclust_replicates \
/scratch/temp/${SLURM_JOB_ID}/${outfilename}

Rscript 01_avg_tree.R \
${ref_hclust} \
/scratch/temp/${SLURM_JOB_ID}/hclust_replicates \
/scratch/temp/${SLURM_JOB_ID}/${outfilename_avg}

cp /scratch/temp/${SLURM_JOB_ID}/${outfilename} ${OUTDIR}
cp /scratch/temp/${SLURM_JOB_ID}/${outfilename_avg} ${OUTDIR}

