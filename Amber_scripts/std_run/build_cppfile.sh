echo "parm protein.prm" > trajin.in
ls -tr md*mdcrd|awk '{print "trajin "$0}' >> trajin.in

echo "center @CA  mass origin" >> trajin.in
echo "image origin center familiar" >> trajin.in
echo "strip :WAT,:Na+,:Cl-" >> trajin.in
echo "strip .H*" >> trajin.in
echo "strip .?H*" >> trajin.in

echo "rms first @CA,C,N out rmsd_bb.1.dat" >> trajin.in
echo "atomicfluct out rmsf_bb.1.dat @CA,C,N byres" >> trajin.in
cpptraj -i trajin.in

