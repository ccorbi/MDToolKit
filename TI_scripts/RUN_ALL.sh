#!/bin/bash

cd MUTATIONS
for i in *
do
	cd $i
	cp ../../prepare-* .
	basename $i > namex
	perl -0pe 's/\n//g' namex | sed 's:^.\(.*\).$:\1:' > res_number
	res_num_var=$(<res_number)
	mut_id_var=$(<namex)
	comp_ions_id_var=$(more ions_to_mask.native-comp.pdb | sed 's/.$//' | awk '{print ","$0}')
	lig_ions_id_var=$(more ions_to_mask.native-lig.pdb | sed 's/.$//' | awk '{print ","$0}')
	sed 's/res_no/'$res_num_var'/g' ../../GENERIC_RUN-TI.sh | sed 's/mut_id/'$mut_id_var'/g' | sed 's/comp_ions_id/'$comp_ions_id_var'/g' | sed 's/lig_ions_id/'$lig_ions_id_var'/g' > "$mut_id_var"_RUN-TI.sh
	chmod +x "$mut_id_var"_RUN-TI.sh
	./"$mut_id_var"_RUN-TI.sh
	cd ..
done
