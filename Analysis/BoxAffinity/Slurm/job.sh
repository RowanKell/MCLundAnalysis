#!/bin/tcsh
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=AffinityCalc
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --chdir=/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis
#SBATCH --output=/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Slurm/output/%x.out
#SBATCH --error=/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Slurm/error/%x.err
source /w/hallb-scshelf2102/clas12/users/rojokell/venv/bin/activate.csh
cd Analysis/BoxAffinity
python3 main.py --no-useArgs
