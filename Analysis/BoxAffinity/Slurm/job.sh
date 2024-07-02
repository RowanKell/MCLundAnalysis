#!/bin/tcsh
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=AffinityCalc
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --chdir=/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis
#SBATCH --output=/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Slurm/output/%x.out
#SBATCH --error=/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Slurm/error/%x.err
source /w/hallb-scshelf2102/clas12/users/rojokell/venv/bin/activate.csh
cd Analysis/BoxAffinity
python3 main.py --outRootName /root_files/May24_all_events.root --plotName driver_single_pion_all_events.pdf --xlsxFileName xlsx/May24_all_events --fileFromLundAnalysis Files_Spring_24/May_24/file_0_all_events.root --useDriver True --useArgs True 
