#!/bin/bash

### first use sed to rename inputname to - for example "r2-bug", inputnative to eg' "r2-native",
### and :mutantres to the mutated res id - in this example - :2

N='2'

sed 's/name-x/mut_id-native/g' extend-complex-step1_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nocomp_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./complex_run/STEP1/extend-complex-step1_script.sh
sed 's/name-x/mut_id-native/g' extend-complex-step2_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nocomp_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./complex_run/STEP2/extend-complex-step2_script.sh
sed 's/name-x/mut_id-native/g' extend-complex-step3_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nocomp_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./complex_run/STEP3/extend-complex-step3_script.sh
sed 's/name-x/mut_id-native/g' extend-ligand-step1_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nolig_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./ligand_run/STEP1/extend-ligand-step1_script.sh
sed 's/name-x/mut_id-native/g' extend-ligand-step2_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nolig_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./ligand_run/STEP2/extend-ligand-step2_script.sh
sed 's/name-x/mut_id-native/g' extend-ligand-step3_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':res_nolig_ions_id/g" | sed "s/mask_noncanon=':res/mask_noncanon=':res_no/g" > ./ligand_run/STEP3/extend-ligand-step3_script.sh



cd complex_run/STEP1
chmod +x *.sh
./extend-complex-step1_script.sh ${N}


cd ../../
cd complex_run/STEP2
chmod +x *.sh
./extend-complex-step2_script.sh ${N}


cd ../../
cd complex_run/STEP3
chmod +x *.sh
./extend-complex-step3_script.sh ${N}




cd ../../
cd ligand_run/STEP1
chmod +x *.sh
./extend-ligand-step1_script.sh ${N}

cd ../../
cd ligand_run/STEP2
chmod +x *.sh

./extend-ligand-step2_script.sh ${N}


cd ../../
cd ligand_run/STEP3
chmod +x *.sh
./extend-ligand-step3_script.sh ${N}



cd ../../



