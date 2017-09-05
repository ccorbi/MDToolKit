#!/bin/bash

## first change Job-name to the actual job name with sed
## then change Native to Job-res-native with sed

mkdir STEP1 STEP2 STEP3
rm ./STEP*/*

cp ../complex_run/STEP1/*.out ./STEP1/
cp ../ligand_run/STEP1/*.out ./STEP1/
cp ../complex_run/STEP2/*.out ./STEP2/
cp ../ligand_run/STEP2/*.out ./STEP2/
cp ../complex_run/STEP3/*.out ./STEP3/
cp ../ligand_run/STEP3/*.out ./STEP3/

rm ./STEP*/*equi*
rm ./STEP*/*min*

#sed 's/name-x/GEN-native/g' TIMD_generic-comp_extraction.py > TIMD_GEN-native-comp_extraction.py
#sed 's/name-x/GEN-native/g' TIMD_generic-lig_extraction.py > TIMD_GEN-native-lig_extraction.py

cp extraction.py ./STEP1/
cp extraction.py ./STEP2/
cp extraction.py ./STEP3/

cd ./STEP1
python extraction.py native-comp_prod_v0_l
python extraction.py native-lig_prod_v0_l
mv report.native-comp_prod_v0_l comp_step1.report
mv report.native-lig_prod_v0_l lig_step1.report
cd ../

cd ./STEP2
python extraction.py $1-native-comp_prod_v1_l
python extraction.py  $1-native-lig_prod_v1_l
mv report.$1-native-comp_prod_v1_l comp_step2.report
mv report.$1-native-lig_prod_v1_l lig_step2.report
cd ../

cd ./STEP3
python extraction.py  $1-native-comp_prod_v1_l
python extraction.py  $1-native-lig_prod_v1_l
mv report.$1-native-comp_prod_v1_l comp_step3.report
mv report.$1-native-lig_prod_v1_l lig_step3.report
cd ../

cp ./STEP*/*.report .

python TIMD_integration_prep.py

echo "${1}-native result:"

perl calc_dvdl.pl integ_prep_step1 | tail -1 > step1.out
perl calc_dvdl.pl integ_prep_step2 | tail -1 > step2.out
perl calc_dvdl.pl integ_prep_step3 | tail -1 > step3.out

cat step1.out step2.out step3.out > all_steps.out
more all_steps.out
total=$(perl -0pe 's/:/ /g' all_steps.out | perl -0pe 's/\+/ /g' | awk '{print $3}' | tail -3 | awk '{s+=$1}END{print s}')
echo "TOTAL = $total kcal/mol"
