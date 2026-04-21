#!/bin/bash
#SBATCH --job-name=sara_hmm
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=5
#SBATCH --output=slurm_outputs/array_job_out_%A.txt
#SBATCH --error=slurm_outputs/array_job_err_%A.txt

module load matlab
srun matlab -nodisplay -batch "triton_dual_from_group_sensor"

seff $SLURM_JOBID
