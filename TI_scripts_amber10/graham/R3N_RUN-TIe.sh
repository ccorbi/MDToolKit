#!/bin/bash

### first use sed to rename inputname to - for example "r2-bug", inputnative to eg' "r2-native",
### and :mutantres to the mutated res id - in this example - :2

N='3'

sed 's/name-x/R3N-native/g' extend-complex-step1_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':3,217/g" | sed "s/mask_noncanon=':res/mask_noncanon=':3/g" > ./complex_run/STEP1/extend-complex-step1_script.sh
sed 's/name-x/R3N-native/g' extend-complex-step2_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':3,217/g" | sed "s/mask_noncanon=':res/mask_noncanon=':3/g" > ./complex_run/STEP2/extend-complex-step2_script.sh
sed 's/name-x/R3N-native/g' extend-complex-step3_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':3,217/g" | sed "s/mask_noncanon=':res/mask_noncanon=':3/g" > ./complex_run/STEP3/extend-complex-step3_script.sh
sed 's/name-x/R3N-native/g' extend-ligand-step1_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':3,20/g" | sed "s/mask_noncanon=':res/mask_noncanon=':3/g" > ./ligand_run/STEP1/extend-ligand-step1_script.sh
sed 's/name-x/R3N-native/g' extend-ligand-step2_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':3,20/g" | sed "s/mask_noncanon=':res/mask_noncanon=':3/g" > ./ligand_run/STEP2/extend-ligand-step2_script.sh
sed 's/name-x/R3N-native/g' extend-ligand-step3_script.sh | sed 's/x-native/native/g' | sed "s/mask_native=':res/mask_native=':3,20/g" | sed "s/mask_noncanon=':res/mask_noncanon=':3/g" > ./ligand_run/STEP3/extend-ligand-step3_script.sh


cp R3N-native* complex_run/STEP1/
cp R3N-native* complex_run/STEP2/
cp R3N-native* complex_run/STEP3/
cp R3N-native* ligand_run/STEP1/
cp R3N-native* ligand_run/STEP2/
cp R3N-native* ligand_run/STEP3/

cp native* complex_run/STEP1/
cp native* complex_run/STEP2/
cp native* complex_run/STEP3/
cp native* ligand_run/STEP1/
cp native* ligand_run/STEP2/
cp native* ligand_run/STEP3/



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

