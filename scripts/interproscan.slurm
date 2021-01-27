#!/bin/bash
#SBATCH --partition=production # partition to submit to
#SBATCH --job-name=interproscan # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=16 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem-per-cpu=2000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=2-00:30:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=ips_%A.out # STDOUT
#SBATCH --error=ips_%A.err # STDERR
#SBATCH --mail-user=jnmaloof@ucdavis.edu #
#SBATCH --mail-type=ALL #

# This will be run once for a single process

/bin/hostname

start=`date +%s`

# Load modules

module load interproscan/5.30-69.0
# set up envrionment

export IPS_CONFIG_DIR=/share/malooflab/Julin/git/TanOak/ipr_config

# Change directory

cd /share/malooflab/Julin/git/TanOak

interproscan.sh --input input/PanGenome_ONT_pilon.all.maker.proteins.fasta \
  --output-file-base output/PanGenome_ONT_V1 \
  --goterms \
  --iprlookup \
  --pathways \
  --formats TSV \
  --cpu 15
  

end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion
