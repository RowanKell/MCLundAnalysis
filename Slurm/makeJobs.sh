#!/bin/bash

workdir="/work/clas12/users/rojokell/MCLundAnalysis"
hipodir="/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/"
#hipodir="/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v0/nobkg_10604MeV/" #THERE ARE FILES HERE, NOWHERE ELSE IN torus-1, nothing in torus+1 either
outputdir="/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Slurm/March_2/run_1/"
rootname="file_"
processdir="/work/clas12/users/rojokell/MCLundAnalysis/"
processcodename="LundAnalysis.C"
runJobs="${workdir}/Slurm/runJobs.sh"
touch $runJobs
chmod +x $runJobs
echo " " > $runJobs

i=0

for hipofile in "$hipodir"/*
do
    file="${workdir}/Slurm/shells/${rootname}${i}.sh"
    touch $file
    echo "#!/bin/tcsh" > $file
    echo "#SBATCH --account=clas12" >> $file
    echo "#SBATCH --partition=production" >> $file
    echo "#SBATCH --mem-per-cpu=4000" >> $file
    echo "#SBATCH --job-name=${rootname}${i}" >> $file
    echo "#SBATCH --cpus-per-task=4" >> $file
    echo "#SBATCH --time=24:00:00" >> $file
    echo "#SBATCH --chdir=${workdir}" >> $file
    echo "#SBATCH --output=${workdir}/Slurm/output/%x-%j-%N.out" >> $file
    echo "#SBATCH --error=${workdir}/Slurm/output/%x-%j-%N.err" >> $file
    echo "echo ${workdir}" >> $file
    echo "source /group/clas12/packages/setup.csh" >> $file
    echo "module load clas12/pro" >> $file
    echo "set CLAS12ROOT=/w/hallb-scshelf2102/clas12/users/rojokell/clas12root" >> $file
    echo "set CCDB_HOME=${CLAS12ROOT}/ccdb" >> $file
#    echo "source ${CCDB_HOME}/environment.csh" >> $file
    echo "cd ${processdir}" >> $file    
    echo "clas12root ${processcodename}\\(\\\"${hipofile}\\\",\\\"${outputdir}/${rootname}${i}.root\\\"\\)" >> $file   
    echo "sbatch shells/${rootname}${i}.sh" >> $runJobs
    i=$((i+1))
done
