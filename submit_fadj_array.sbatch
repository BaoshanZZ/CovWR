#!/bin/bash
#SBATCH --job-name=wrcov_fadj_n1000
#SBATCH --array=1-14%14
#SBATCH --cpus-per-task=27
#SBATCH --mem-per-cpu=3G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/home/bz91/WRCov/slurm_logs/fadj_n1000/%A_%a.out
#SBATCH --error=/hpc/home/bz91/WRCov/slurm_logs/fadj_n1000/%A_%a.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /hpc/home/bz91
mkdir -p /hpc/home/bz91/WRCov/slurm_logs/fadj_n1000

Rscript WRCov/R/run_scenario_batch.R \
  --scenario-id="${SLURM_ARRAY_TASK_ID}" \
  --cores="${SLURM_CPUS_PER_TASK}" \
  --seed=20260709 \
  --output-dir=/hpc/home/bz91/WRCov/simulation_output_fadj_n1000
