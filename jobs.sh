#!/bin/sh
#SBATCH --job-name=AcinetobacterGenomes
#SBATCH --array=1-6041%40

echo $HOSTNAME
source ~/VirtualEnv/SOMseq/bin/activate
module load Python/3.8.2
roortdir=`/c7/sc`

sleep 1
taxid=`awk '(NR=='"$SLURM_ARRAY_TASK_ID"'){print}' /work/ifilella/bytaxid/TssM/acinetobacters.txids`

