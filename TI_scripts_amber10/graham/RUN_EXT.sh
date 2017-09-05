#!/bin/bash

cd MUTATIONS
for i in   K4A R3A  L5A L7A Q10A T8A Q11A  M12A
do
	cd $i
	cp ../../extend-* .
	basename $i > namex
	perl -0pe 's/\n//g' namex | sed 's:^.\(.*\).$:\1:' > res_number
	sed 's/res_no/'$res_num_var'/g' ../../GENERIC_RUN-TIe.sh | sed 's/mut_id/'$mut_id_var'/g' | sed 's/comp_ions_id/'$comp_ions_id_var'/g' | sed 's/lig_ions_id/'$lig_ions_id_var'/g' > "$mut_id_var"_RUN-TIe.sh
	chmod +x "$mut_id_var"_RUN-TIe.sh
	./"$mut_id_var"_RUN-TIe.sh
	cd ..
done
