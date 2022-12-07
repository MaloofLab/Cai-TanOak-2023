#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=ortho_align_blastp # Job name
#SBATCH --array=1-9
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem-per-cpu=4000 # in MB default is 2GB per cpu
#SBATCH --time=03:00:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=ortho_align_%A_%a.out # STDOUT
#SBATCH --error=ortho_align_%A_%a.err # STDERR
#SBATCH --mail-user=prmagalang@ucdavis.edu #
#SBATCH --mail-type=ALL #

/bin/hostname

start=`date +%s`

module load java/jdk11.0

cd /share/malooflab/Julin/git/TanOak/output/ortho/blastp_ortholog 

# Identify each array run
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# make directories if needed
mkdir -p aligned_seqs/aa
mkdir -p aligned_seqs/nt

# set up starting indices
start_index=$(( ($SLURM_ARRAY_TASK_ID - 1) * 1000 ))

# want to run 1000 alignments per job except for the last
job_length=1000
if [[ $SLURM_ARRAY_TASK_ID = 9 ]]; then
    job_length=496
fi

# take a subset of fasta files
all_files=(raw_seqs/*.fasta)
files=("${all_files[@]:start_index:job_length}")

# align each fasta file
for file in "${files[@]}"; do

    file_nopath=$(echo $file | cut -d '/' -f 2)
    file_noext="${file_nopath%.*}"

    out_nt="aligned_seqs/nt/${file_noext}.nt.fasta"
    out_aa="aligned_seqs/aa/${file_noext}.aa.fasta"
    java -jar /share/malooflab/Paulo/tools/macse_v2.05.jar -prog alignSequences -seq $file -out_NT $out_nt -out_AA $out_aa

done

# stats
end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion

