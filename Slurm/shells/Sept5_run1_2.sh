#!/bin/tcsh
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=Sept5_run1_2
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --chdir=/work/clas12/users/rojokell/MCLundAnalysis
#SBATCH --output=/work/clas12/users/rojokell/MCLundAnalysis/Slurm/output/%x-%j-%N.out
#SBATCH --error=/work/clas12/users/rojokell/MCLundAnalysis/Slurm/output/%x-%j-%N.err
echo /work/clas12/users/rojokell/MCLundAnalysis
source /group/clas12/packages/setup.csh
module load clas12/pro
set CLAS12ROOT=/w/hallb-scshelf2102/clas12/users/rojokell/clas12root
set CCDB_HOME=/group/clas12/packages/clas12root/1.7.8.c/ccdb
cd /work/clas12/users/rojokell/MCLundAnalysis/
clas12root LundAnalysis.C\(\"/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV//45nA_job_3051_2.hipo\",\"/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Slurm/Sept_5//Sept5_run1_2.root\"\)
