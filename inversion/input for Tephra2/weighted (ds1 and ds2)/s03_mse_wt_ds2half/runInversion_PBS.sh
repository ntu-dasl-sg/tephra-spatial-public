#!/bin/bash

#PBS -N mse_wt
#PBS -j oe
#PBS -V
#PBS -m n
#PBS -M maricarl001@ntu.edu.sg
#PBS -l nodes=1:ppn=12
#PBS -q q12
#PBS

# Load modules
module load python/2.7.4
module load openmpi/1.4.5-gnu

# Usage on OpenPBS
# qsub -t 9-11 T2inversion_openPBS.sh    Batch, mass between 10^9 and 10^11
# qsub runMass.sh                        Single run

#____________________________________________________________
#
# DO NOT CHANGE ANYTHING BELOW THIS LINE
#____________________________________________________________

function sum() {
perl -e "print $1 + $2"
}

# With OpenPBS, got to cd to folder
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

source inversionConfig.conf

chmod 755 ../_scripts/genConfig.py
chmod 755 ../_scripts/genConfigForward.py

# Check the seed
if [ "$SEED" = "-1" ]; then
	RANGE=20
	SEED=$RANDOM
	let "SEED %= $RANGE"
fi

if [ "$BATCH" = "1" ]; then
	echo "BATCH"
	mass1=$PBS_ARRAYID
	mass2=$(($mass1 + deltaMass))

	for ht in `seq $minHt $deltaHt $maxHt`; do
		outDir=mass${mass1}_ht$(($ht/1000))
		mkdir $outDir
		cp $inputFile $outDir/inversionInput.txt
		cd $outDir

		touch h_q_rmse1.dat
		../../_scripts/genConfig.py $ht `sum $ht $incrHt` "1e$mass1" "1e$mass2" $ventE $ventN $ventA $minDiff $maxDiff $eddy $minMedPhi $maxMedPhi $minSigPhi  $maxSigPhi $minAlpha $maxAlpha $minBeta $maxBeta $minFTT $maxFTT $minWindSpeed $maxWindSpeed $minWindDir $maxWindDir $plumeModel $fixedWind $windLevels $colSteps $partSteps $lithicDensity $pumiceDensity $minPhi $maxPhi $fitTest $SEED
		date
         	mpirun  ../../../tephraMSEweighted_inversion tmp.conf ../$inputFile ../$windFile
		date

		rm node_*

		# Write out a configuration file for tephra2
		../../_scripts/genConfigForward.py
		# Run the tephra2 forward model
		mpirun -np 1 -machinefile $PBS_NODEFILE ../../../tephra2-MSEweighted tephra2.conf ../$gridFile wind_levels.out > tephra2.out

		cd ..
	done

elif [ "$BATCH" = "2" ]; then
	echo "SINGLE with variable seed"
	SEED=$PBS_ARRAYID

	outDir=mass${minMass}_ht$(($minHt/1000))_$SEED
	mkdir $outDir

	cp $inputFile $outDir/inversionInput.txt
	cd $outDir

	touch h_q_rmse1.dat
	../../_scripts/genConfig.py $minHt $maxHt "1e$minMass" "1e$maxMass" $ventE $ventN $ventA $minDiff $maxDiff $eddy $minMedPhi $maxMedPhi $minSigPhi  $maxSigPhi $minAlpha $maxAlpha $minBeta $maxBeta $minFTT $maxFTT $minWindSpeed $maxWindSpeed $minWindDir $maxWindDir $plumeModel $fixedWind $windLevels $colSteps $partSteps $lithicDensity $pumiceDensity $minPhi $maxPhi $fitTest $SEED
	date
		mpirun ../../../tephraMSEweighted_inversion tmp.conf ../$inputFile ../$windFile
	date

	rm node_*

	# Write out a configuration file for tephra2
	../../_scripts/genConfigForward.py
	# Run the tephra2 forward model
	mpirun -np 1 -machinefile $PBS_NODEFILE ../../../tephra2-MSEweighted tephra2.conf ../$gridFile wind_levels.out >tephra2.out

	cd ..

else
	echo "SINGLE"
	outDir=mass${minMass}_ht$(($minHt/1000))
	mkdir $outDir
	cp $inputFile $outDir/inversionInput.txt
	cd $outDir

	touch h_q_rmse1.dat

	../../_scripts/genConfig.py $minHt $maxHt "1e$minMass" "1e$maxMass" $ventE $ventN $ventA $minDiff $maxDiff $eddy $minMedPhi $maxMedPhi $minSigPhi  $maxSigPhi $minAlpha $maxAlpha $minBeta $maxBeta $minFTT $maxFTT $minWindSpeed $maxWindSpeed $minWindDir $maxWindDir $plumeModel $fixedWind $windLevels $colSteps $partSteps $lithicDensity $pumiceDensity $minPhi $maxPhi $fitTest $SEED
	date
	mpirun -np $NCPUS -machinefile ../../../tephraMSEweighted_inversion tmp.conf ../$inputFile ../$windFile
	date

	rm node_*

	# Write out a configuration file for tephra2
	../../_scripts/genConfigForward.py
	# Run the tephra2 forward model
	mpirun -np 1 -machinefile $PBS_NODEFILE ../../../tephra2-MSEweighted tephra2.conf ../$gridFile wind_levels.out >tephra2.out

	cd ..
fi
