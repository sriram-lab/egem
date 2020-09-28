#!/bin/bash
#@Author: Scott Campit

SBATCH --job-name=PredMetFromGCP_1
SBATCH --mail-user=scampit@umich.edu
SBATCH --mail-type=BEGIN,END,FAIL
SBATCH --cpus-per-task=1
SBATCH --nodes=1
SBATCH --ntasks-per-node=1
SBATCH --mem-per-cpu=1000m 
SBATCH --time=10:00
SBATCH --account=test
SBATCH --partition=standard
SBATCH --output=/home/%u/%x-%j.log
