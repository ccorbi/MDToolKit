#!/bin/bash

### NOTES ###
# 1. mutations to HIP require some manual work to get correct protonation state
# 2. input files: native-comp.pdb & native-lig.pdb need to contain at least 2 x NA & 2 x CL to
# ensure that all charge change setups work.
# 2b. in tleap setup, solvate first, add ions second.

source /home/smona/bin/amber14/amber.sh # change to your amber.sh path

mkdir MUTATIONS-temp
cd MUTATIONS-temp

mkdir H33T H33A  # mutant list, EDIT eg: D23R G33T Y54F

for i in *
do
	cd $i
	cp ../../native-* .
	# parse mutation info 
	basename $i > namex
	head -c 1 namex > name_from
	perl -0pe 's/\n//g' namex | tail -c 1 > name_to
	perl -0pe 's/\n//g' namex | sed 's:^.\(.*\).$:\1:' > res_number

	# convert amino acid names
	for n in name_*
	do
		if grep -q R "$n"; then
		sed 's/./ARG/g' $n > full_$n
		elif grep -q H "$n"; then
		sed 's/./HIP/g' $n > full_$n
		elif grep -q K "$n"; then
		sed 's/./LYS/g' $n > full_$n
		elif grep -q D "$n"; then
		sed 's/./ASP/g' $n > full_$n
		elif grep -q E "$n"; then
		sed 's/./GLU/g' $n > full_$n
		elif grep -q S "$n"; then
		sed 's/./SER/g' $n > full_$n
		elif grep -q T "$n"; then
		sed 's/./THR/g' $n > full_$n
		elif grep -q N "$n"; then
		sed 's/./ASN/g' $n > full_$n
		elif grep -q Q "$n"; then
		sed 's/./GLN/g' $n > full_$n
		elif grep -q C "$n"; then
		sed 's/./CYS/g' $n > full_$n
		elif grep -q G "$n"; then
		sed 's/./GLY/g' $n > full_$n
		elif grep -q P "$n"; then
		sed 's/./PRO/g' $n > full_$n
		elif grep -q A "$n"; then
		sed 's/./ALA/g' $n > full_$n
		elif grep -q V "$n"; then
		sed 's/./VAL/g' $n > full_$n
		elif grep -q I "$n"; then
		sed 's/./ILE/g' $n > full_$n
		elif grep -q L "$n"; then
		sed 's/./LEU/g' $n > full_$n
		elif grep -q M "$n"; then
		sed 's/./MET/g' $n > full_$n
		elif grep -q F "$n"; then
		sed 's/./PHE/g' $n > full_$n
		elif grep -q Y "$n"; then
		sed 's/./TYR/g' $n > full_$n
		elif grep -q W "$n"; then
		sed 's/./TRP/g' $n > full_$n
		fi
	done
	for t in native-*
	do
		# remove non-backbone atoms from mutant residue
		res_num_value=$(<res_number)
		full_name_from_value=$(<full_name_from)	
		full_name_to_value=$(<full_name_to)
		awk '$5 != "'$res_num_value'" || ($5 == "'$res_num_value'" && ($3 == "N" || $3 == "O" || $3 == "C" || $3 == "CA"))' $t > temp1
		
		# rename mutant residue
		echo "at position: $res_num_value"
		echo "mutating residue from $full_name_from_value"
		echo "to $full_name_to_value"
		sed -i -e '/'" $res_num_value "'/ s/'"$full_name_from_value"'/'"$full_name_to_value"'/g' temp1
		grep -v "ATOM    $res_num_value " temp1 > $i-$t
		cp temp1 $i-$t
		cat << EOF > make_$i.leap
source leaprc.ff14SB
loadamberparams frcmod.ionsjc_tip3p

complex = loadpdb native-comp.pdb
setbox complex vdw
savepdb complex native-comp_leap.pdb
saveamberparm complex native-comp.prm native-comp.rst
charge complex

ligand = loadpdb native-lig.pdb
setbox ligand vdw
savepdb ligand native-lig_leap.pdb
saveamberparm ligand native-lig.prm native-lig.rst
charge ligand

complex = loadpdb $i-native-comp.pdb
setbox complex vdw
savepdb complex $i-native-comp_leap.pdb
saveamberparm complex $i-native-comp.prm $i-native-comp.rst
charge complex

ligand = loadpdb $i-native-lig.pdb
setbox ligand vdw
savepdb ligand $i-native-lig_leap.pdb
saveamberparm ligand $i-native-lig.prm $i-native-lig.rst
charge ligand

quit
EOF
		#sed 's/GEN/'$i'/g' ../../make_nc-generic_input.leap > make_"$i"_input.leap
		
		# Removing unwanted counterions (lower-most ions removed, by reversing line order, removing first instance of ion(s), restoring line order)
		cp $i-$t removed_ions
		
		if [[ "$full_name_from_value" =~ ^(GLU|ASP)$ ]] && ! [[ "$full_name_to_value" =~ ^(GLU|ASP|ARG|LYS|HIP)$ ]]; then
		tac $i-$t | awk '/NA/ && !p {p++;next}1' | tac > removed_ions
		elif ! [[ "$full_name_from_value" =~ ^(GLU|ASP|ARG|LYS|HIP)$ ]] && [[ "$full_name_to_value" =~ ^(ARG|LYS|HIP)$ ]]; then
		tac $i-$t | awk '/NA/ && !p {p++;next}1' | tac > removed_ions
		elif [[ "$full_name_from_value" =~ ^(ARG|LYS|HIP)$ ]] && ! [[ "$full_name_to_value" =~ ^(GLU|ASP|ARG|LYS|HIP)$ ]]; then
		tac $i-$t | awk '/CL/ && !p {p++;next}1' | tac > removed_ions
		elif ! [[ "$full_name_from_value" =~ ^(GLU|ASP|ARG|LYS|HIP)$ ]] && [[ "$full_name_to_value" =~ ^(GLU|ASP)$ ]]; then
		tac $i-$t | awk '/CL/ && !p {p++;next}1' | tac > removed_ions
		elif [[ "$full_name_from_value" =~ ^(ARG|LYS|HIP)$ ]] && [[ "$full_name_to_value" =~ ^(GLU|ASP)$ ]]; then
		tac $i-$t | awk '/CL/ && !p {p++;next}1' | awk '/CL/ && !p {p++;next}1' | tac > removed_ions
		elif [[ "$full_name_from_value" =~ ^(GLU|ASP)$ ]] && [[ "$full_name_to_value" =~ ^(ARG|LYS|HIP)$ ]]; then
		tac $i-$t | awk '/NA/ && !p {p++;next}1' | awk '/NA/ && !p {p++;next}1' | tac > removed_ions
		fi

		mv removed_ions $i-$t

		# IDENTIFY IONS TO MASK	
		cat $i-$t $t | grep -E 'NA|CL'| awk '{print $5}' | sort | uniq -u | perl -0pe 's/\n/,/g' > ions_to_mask.$t
	done
	tleap -f make_"$i".leap
	rm *ame* *numb* temp*
	cd ..
done
cd ..
