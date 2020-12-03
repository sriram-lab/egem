#!/bin/bash
#SBATCH --job-name=met2gcp_regression
#SBATCH --mail-user=scampit@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=./output.log
#SBATCH --error=./error.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1g
#SBATCH --cpus-per-task=16
#SBATCH --time=03-00:00:00
#SBATCH --account=lsa1
#SBATCH --partition=standard

module load python3
source /nfs/turbo/umms-csriram/scampit/Envs/python/ml/bin/activate
pip install --upgrade pip --user
pip3 install -r /nfs/turbo/umms-csriram/scampit/Software/egem/srv/python/ml/requirements.txt --user
python3 /nfs/turbo/umms-csriram/scampit/Software/egem/srv/python/ml/regression.py
