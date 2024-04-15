#!/bin/bash
current_date=$(date +"%B_%d")
workdir="/work/clas12/users/rojokell/MCLundAnalysis"
hipodir="/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/"
# hipodir="/lustre19/expphy/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV"
slurm_output="${workdir}/OutputFiles/Slurm_Spring_24"
daydir="${slurm_output}/${current_date}"
#USER SET VALUES
outputdir="${daydir}/Run_1_single_pion/"
single_pion_flag=1 #1 = true, bin by single pion values | 0 = false, bin by dihadron values

out_folder="/work/clas12/users/rojokell/MCLundAnalysis/Slurm/output/output${current_date}"
error_folder="/work/clas12/users/rojokell/MCLundAnalysis/Slurm/error/error${current_date}"

rootname="file_"
processdir="/work/clas12/users/rojokell/MCLundAnalysis/"
processcodename="LundAnalysis.C"
runJobs="${workdir}/Slurm/runJobs.sh"
touch $runJobs
chmod +x $runJobs
echo " " > $runJobs
echo $daydir
i=0

if [ ! -d "$daydir" ]; then
  mkdir "$daydir"
fi

if [ ! -d "$outputdir" ]; then
  mkdir "$outputdir"
fi

if [ ! -d "$out_folder" ]; then
  mkdir "$out_folder"
fi

if [ ! -d "$error_folder" ]; then
  mkdir "$error_folder"
fi

for hipofile in "$hipodir"/*
do
    file="${workdir}/Slurm/shells/${rootname}${i}.sh"
    touch $file
    echo "#!/bin/tcsh" > $file
    echo "#SBATCH --account=clas12" >> $file
    echo "#SBATCH --partition=production" >> $file
    echo "#SBATCH --mem-per-cpu=4000" >> $file
    echo "#SBATCH --job-name=${rootname}${i}" >> $file
    echo "#SBATCH --cpus-per-task=1" >> $file
    echo "#SBATCH --time=24:00:00" >> $file
    echo "#SBATCH --chdir=${workdir}" >> $file
    echo "#SBATCH --output=${out_folder}/%x.out" >> $file
    echo "#SBATCH --error=${error_folder}/%x.err" >> $file
    echo "module purge -f" >> $file
    echo "source /u/home/rojokell/.cshrc" >> $file
    echo "cd ${processdir}" >> $file    
    echo "clas12root ${processcodename}\\(\\\"${hipofile}\\\",\\\"${outputdir}/${rootname}${i}.root\\\",${single_pion_flag}\\)" >> $file   
    echo "sbatch shells/${rootname}${i}.sh" >> $runJobs
    i=$((i+1))
done
