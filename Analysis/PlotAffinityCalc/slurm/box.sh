#!/bin/tcsh
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=box_py_April_14_Run_1_single_pion
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --chdir=/work/clas12/users/rojokell/MCLundAnalysis
#SBATCH --output=/work/clas12/users/rojokell/MCLundAnalysis/Analysis/slurm/output/%x.out
#SBATCH --error=/work/clas12/users/rojokell/MCLundAnalysis/Analysis/slurm/error/%x.err
source /u/home/rojokell/.cshrc
source /work/clas12/users/rojokell/MCLundAnalysis/venv/bin/activate.csh
python3 ./Analysis/BoxAffinity/Box.py
