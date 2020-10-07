#!/bin/bash
#SBATCH --job-name=nonlinear_regress1
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=03-00:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

source /nfs/turbo/umms-csriram/scampit/PyEnvs/ml/bin/activate
python3 /nfs/turbo/umms-csriram/scampit/Software/egem/srv/python/ml/nonlinear_regression.py
