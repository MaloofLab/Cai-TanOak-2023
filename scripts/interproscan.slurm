#!/bin/bash
#SBATCH --partition=development # partition to submit to
#SBATCH --job-name=interproscan # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem-per-cpu=8000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-00:10:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=ips_%A.out # STDOUT
#SBATCH --error=ips_%A.err # STDERR
#SBATCH --mail-user=jnmaloof@ucdavis.edu #
#SBATCH --mail-type=ALL #

# This will be run once for a single process

/bin/hostname

start=`date +%s`

# Load modules

module load interproscan/5.30-69.0

# Change directory

cd /share/malooflab/Julin/git/TanOak

interproscan --input input/PanGenome_ONT_pilon.all.maker.proteins.fasta.gz \
  --output-file-base output/PanGenome_ONT_V1 \
  --goterms \
  --iprlookup \
  --pathways
  

end=`date +%s`
runtime=$((end-start))
echo $runtime seconds to completion
