#!/bin/bash
#SBATCH --job-name=nonlinear_regress1
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=./output.txt
#SBATCH --error=./error.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=16
#SBATCH --time=03-00:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

source /nfs/turbo/umms-csriram/scampit/PyEnvs/ml/bin/activate
python3 /nfs/turbo/umms-csriram/scampit/Software/egem/srv/python/ml/regression.py
