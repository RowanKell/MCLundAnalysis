#!/bin/bash

workdir="/work/clas12/users/rojokell/MCLundAnalysis"
hipodir="/mss/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV"
runJobs="${workdir}/jcache/runJcache.sh"

touch $runJobs
chmod +x $runJobs
echo -n "jcache default " > $runJobs

i=0


for hipofile in "$hipodir"/*
do
    echo -n "${hipofile} " >> $runJobs
    i=$((i+1))
done
