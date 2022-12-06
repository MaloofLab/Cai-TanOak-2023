#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=ortho_align # Job name
#SBATCH --array=1-8
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem-per-cpu=4000 # in MB default is 2GB per cpu
#SBATCH --time=01:00:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=ortho_align_%A_%a.out # STDOUT
#SBATCH --error=ortho_align_%A_%a.err # STDERR
#SBATCH --mail-user=prmagalang@ucdavis.edu #
#SBATCH --mail-type=ALL #

/bin/hostname

start=`date +%s`

module purge
module load # something java

cd /share/malooflab/Julin/git/TanOak/scripts # change file path later

# Identify each array run
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# set up starting indices
start_index=$(( ($SLURM_ARRAY_TASK_ID - 1) * 1000 ))

# want to run 1000 alignments per job except for the last
job_length=1000
if [[ $SLURM_ARRAY_TASK_ID = 8 ]]; then
    job_length=638
fi

# take a subset of fasta files
all_files=(*.fasta)
files=("{$all_files[@]:start_index:job_length}")

let i=1
total="${#files[@]}"

# align each fasta file
for file in "${files[@]}"; do
    echo "Aligning: $i/$total"

    out_nt=
    out_aa=
    java -jar ~/macse_v2.o5.jar -prog alignSequences -seq $file -out_NT $out_nt -out_AA $out_aa

    let i++
done

# stats
end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion

