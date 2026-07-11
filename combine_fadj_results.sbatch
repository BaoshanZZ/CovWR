#!/bin/bash
#SBATCH --job-name=wrcov_fadj_n1000_combine
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/hpc/home/bz91/WRCov/slurm_logs/fadj_n1000/combine_%j.out
#SBATCH --error=/hpc/home/bz91/WRCov/slurm_logs/fadj_n1000/combine_%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /hpc/home/bz91
mkdir -p /hpc/home/bz91/WRCov/slurm_logs/fadj_n1000

Rscript WRCov/R/combine_batch_results.R \
  --output-dir=/hpc/home/bz91/WRCov/simulation_output_fadj_n1000
