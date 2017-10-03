#!/bin/bash

cd ../MUTATIONS
for i in *
do
	cd $i
	cd analysis
	cp ../../../ANALYSE/extraction.py .
    cp ../../../ANALYSE/TIMD_integration_prep.py .
	cp ../../../ANALYSE/calc* .
	cp ../../../ANALYSE/FML.sh .
	sh FML.sh $i > RESULT_$i.txt
	cd ../../
done
cd ../ANALYSE
mkdir RESULTS
cd RESULTS

cp ../../MUTATIONS/*/analysis/RESULT* .
more RESULT* > FULL_REPORT.txt
cd ../
cp RESULTS/FULL_REPORT.txt .

