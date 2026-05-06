#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --job-name=d2s-magena_all
#SBATCH --time=120:00:00
#SBATCH --partition=general
#SBATCH --account=a_ace
#SBATCH -o d2s-magena_all.o
#SBATCH -e d2s-magena_all.e

echo "#=================== JOB INFO ===================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "#================================================#"

# ############################################################

NJ_method=NJ
ref_matrix=/QRISdata/Q8193/01_Marzonie_CoralSea2024/05_unmapped/06_d2s/60Mbp_80_80_magena_all.d2s.k21.v1.txt
ref_NJ=/scratch/temp/${SLURM_JOB_ID}/ref_nj.awk
outfilename=${NJ_method}_pseudoreplicates_ref_bs.awk
outfilename_avg=${NJ_method}_pseudoreplicates_avg_bs.awk

module load r/4.4.0-gfbf-2023a

Rscript /scratch/project_mnt/S0026/ishida/0_scripts/batch/SYMB_MAGENA/nj_tree.R ${ref_matrix} ${ref_NJ}

INDIR=/QRISdata/Q8193/01_Marzonie_CoralSea2024/05_unmapped/08_d2s_jack/2_60mbp_distmat
cp ${INDIR}/* /scratch/temp/${SLURM_JOB_ID}

mkdir -p /scratch/temp/${SLURM_JOB_ID}/NJ_replicates

for i in /scratch/temp/${SLURM_JOB_ID}/*.txt
do
    ls $i 
    name=$(basename $i .txt)
    echo $name 
    outpath=/scratch/temp/${SLURM_JOB_ID}/NJ_replicates/${name}.awk

    echo $outpath
    Rscript /scratch/project_mnt/S0026/ishida/0_scripts/batch/SYMB_MAGENA/nj_tree.R ${i} ${outpath} 
    # break
done

Rscript /scratch/project_mnt/S0026/ishida/0_scripts/batch/SYMB_MAGENA/Jackknife_param.r \
${ref_NJ} \
/scratch/temp/${SLURM_JOB_ID}/NJ_replicates \
/scratch/temp/${SLURM_JOB_ID}/${outfilename}

Rscript /scratch/project_mnt/S0026/ishida/0_scripts/batch/SYMB_MAGENA/avg_tree.R \
${ref_NJ} \
/scratch/temp/${SLURM_JOB_ID}/NJ_replicates \
/scratch/temp/${SLURM_JOB_ID}/${outfilename_avg}


cp /scratch/temp/${SLURM_JOB_ID}/${outfilename} /scratch/project_mnt/S0026/ishida/symb/Magena/7_D2S

cp /scratch/temp/${SLURM_JOB_ID}/${outfilename_avg} /scratch/project_mnt/S0026/ishida/symb/Magena/7_D2S

