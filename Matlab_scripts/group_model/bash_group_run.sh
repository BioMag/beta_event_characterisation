#!/bin/bash
#SBATCH --job-name=sara_hmm
#SBATCH --time=80:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_outputs/array_job_out_%A.txt
#SBATCH --error=slurm_outputs/array_job_err_%A.txt

module load matlab
srun matlab -nodisplay -batch "make_group_model"

seff $SLURM_JOBID
