#!/bin/bash

### first use sed to rename inputname to - for example "r2-bug", inputnative to eg' "r2-native",
### and :mutantres to the mutated res id - in this example - :2

mkdir complex_run
mkdir ligand_run
mkdir complex_run/STEP1
mkdir complex_run/STEP2
mkdir complex_run/STEP3
mkdir ligand_run/STEP1
mkdir ligand_run/STEP2
mkdir ligand_run/STEP3
mkdir analysis

sed 's/name-x/mut_id-native/g' prepare-complex-step1_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nocomp_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./complex_run/STEP1/prepare-complex-step1_script.sh
sed 's/name-x/mut_id-native/g' prepare-complex-step2_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nocomp_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./complex_run/STEP2/prepare-complex-step2_script.sh
sed 's/name-x/mut_id-native/g' prepare-complex-step3_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nocomp_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./complex_run/STEP3/prepare-complex-step3_script.sh
sed 's/name-x/mut_id-native/g' prepare-ligand-step1_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nolig_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./ligand_run/STEP1/prepare-ligand-step1_script.sh
sed 's/name-x/mut_id-native/g' prepare-ligand-step2_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nolig_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./ligand_run/STEP2/prepare-ligand-step2_script.sh
sed 's/name-x/mut_id-native/g' prepare-ligand-step3_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nolig_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./ligand_run/STEP3/prepare-ligand-step3_script.sh


cp mut_id-native* complex_run/STEP1/
cp mut_id-native* complex_run/STEP2/
cp mut_id-native* complex_run/STEP3/
cp mut_id-native* ligand_run/STEP1/
cp mut_id-native* ligand_run/STEP2/
cp mut_id-native* ligand_run/STEP3/

cp native* complex_run/STEP1/
cp native* complex_run/STEP2/
cp native* complex_run/STEP3/
cp native* ligand_run/STEP1/
cp native* ligand_run/STEP2/
cp native* ligand_run/STEP3/

cd complex_run/STEP1
chmod +x *.sh
./prepare-*
cd ../../
cd complex_run/STEP2
chmod +x *.sh
./prepare-*
cd ../../
cd complex_run/STEP3
chmod +x *.sh
./prepare-*
cd ../../
cd ligand_run/STEP1
chmod +x *.sh
./prepare-*
cd ../../
cd ligand_run/STEP2
chmod +x *.sh
./prepare-*
cd ../../
cd ligand_run/STEP3
chmod +x *.sh
./prepare-*
cd ../../

