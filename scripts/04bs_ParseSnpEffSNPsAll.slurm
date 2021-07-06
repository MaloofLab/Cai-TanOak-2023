#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=R_annoate # Job name
#SBATCH --array=1-10
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem-per-cpu=12000 # in MB default is 2GB per cpu
#SBATCH --time=3-00:00:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=Tan_anno_%A_%a.out # STDOUT
#SBATCH --error=Tan_anno_%A_%a.err # STDERR
#SBATCH --mail-user=jnmaloof@ucdavis.edu #
#SBATCH --mail-type=ALL #

# This will be run once for a single process

/bin/hostname

start=`date +%s`

module load R/4.1.0

cd /share/malooflab/Julin/git/TanOak/scripts

# Identify each array run

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

firstLast=`sed "${SLURM_ARRAY_TASK_ID}q;d" ../input/contig_splits.txt`

Rscript --no-save --no-restore 04bs_ParseSnpEffSNPsAll.R ${firstLast}

end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion
