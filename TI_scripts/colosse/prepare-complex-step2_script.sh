#!/bin/bash
wd=`pwd`
pre0='x-native-comp'
pre1='name-x-comp'
mask_native=':res'
mask_noncanon=':res'

for X in 1 2 3 4 5 6 7 8 9

do

cat << EOF > group_min_l${X}
-O -i mdin_min_v0_l${X} -o ${pre0}_min_v0_l${X}.out -p ${pre0}.prm -c ${pre0}.rst -r ${pre0}_min_v0_l${X}.rst
-O -i mdin_min_v1_l${X} -o ${pre1}_min_v1_l${X}.out -p ${pre1}.prm -c ${pre1}.rst -r ${pre1}_min_v1_l${X}.rst
EOF
cat << EOF > group_equi_l${X}
-O -i mdin_equi_v0_l${X} -o ${pre0}_equi_v0_l${X}.out -p ${pre0}.prm -c ${pre0}_min_v0_l${X}.rst -r ${pre0}_equi_v0_l${X}.rst
-O -i mdin_equi_v1_l${X} -o ${pre1}_equi_v1_l${X}.out -p ${pre1}.prm -c ${pre1}_min_v1_l${X}.rst -r ${pre1}_equi_v1_l${X}.rst
EOF
cat << EOF > group_prod_l${X}
-O -i mdin_prod_v0_l${X} -o ${pre0}_prod_v0_l${X}.out -p ${pre0}.prm -c ${pre0}_equi_v0_l${X}.rst -r ${pre0}_prod_v0_l${X}.rst -x ${pre0}_prod_v0_l${X}.crd
-O -i mdin_prod_v1_l${X} -o ${pre1}_prod_v1_l${X}.out -p ${pre1}.prm -c ${pre1}_equi_v1_l${X}.rst -r ${pre1}_prod_v1_l${X}.rst -x ${pre1}_prod_v1_l${X}.crd
EOF

cat << EOF > mdin_min_v0_l${X}
density minlibration
 &cntrl
  imin = 1, ntmin = 2,	ntx = 1,
  maxcyc=10000,
  ntpr = 100,
  ntf = 1,      ntc = 1,
  ntb = 1,	cut = 9.0,
  icfe=1,	clambda = 0.${X},
EOF

cp mdin_min_v0_l${X} mdin_min_v1_l${X}

cat << EOF >> mdin_min_v0_l${X}
  ifsc=1,
  crgmask='${mask_native}',
  scmask='${mask_native}',
 &end
EOF
cat << EOF >> mdin_min_v1_l${X}
  ifsc=1,
  crgmask='${mask_noncanon}',
  scmask='${mask_noncanon}',
 &end
EOF

cat << EOF > mdin_equi_v0_l${X}
density equilibration
 &cntrl
  imin = 0,	ntx = 1,	irest = 0,
  ntpr = 2500,	ntwr = 10000,	ntwx = 0,
  ntf = 1,      ntc = 1,
  ntb = 2,	cut = 9.0,
  nstlim = 50000,	dt = 0.001,
  temp0 = 300.0,	ntt = 3,	gamma_ln = 5,
  ntp = 1,	pres0 = 1.0,	taup = 0.2,
  icfe=1,	clambda = 0.${X},
EOF

cp mdin_equi_v0_l${X} mdin_equi_v1_l${X}

cat << EOF >> mdin_equi_v0_l${X}
  ifsc=1,
  crgmask='${mask_native}',
  scmask='${mask_native}',
 &end
EOF
cat << EOF >> mdin_equi_v1_l${X}
  ifsc=1,
  crgmask='${mask_noncanon}',
  scmask='${mask_noncanon}',
 &end
EOF

cat << EOF > mdin_prod_v0_l${X}
NPT production
 &cntrl
  imin = 0,	ntx = 5,	irest = 1,
  ntpr = 10000,	ntwr = 100000,	ntwx = 10000,
  ntf = 1,	ntc = 1,
  ntb = 2,	cut = 9.0,
  nstlim = 450000,	dt = 0.001,
  temp0 = 300.0,	ntt = 3,	gamma_ln = 2,
  ntp = 1,	pres0 = 1.0,	taup = 2.0,
  icfe=1,       clambda = 0.${X},
EOF

cp mdin_prod_v0_l${X} mdin_prod_v1_l${X}

cat << EOF >> mdin_prod_v0_l${X}
  ifsc=1,
  crgmask='${mask_native}',
  scmask='${mask_native}',
 &end
EOF
cat << EOF >> mdin_prod_v1_l${X}
  ifsc=1,
  crgmask='${mask_noncanon}',
  scmask='${mask_noncanon}',
 &end
EOF

cat << EOF > run.pbs.${X}

#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00

# consider changing the name below to something more explanatory
#PBS -N ${X}_complex-run

# PMK RAPI
#PBS -A cvj-675-aa # Identifiant Rap; ID

# change to the submission directory
cd \$PBS_O_WORKDIR

module purge
module load intel intelmpi fftw/3.3.0-intel-impi amber/10.0.30

# Run the parallel version of sander, using all 8 cores in the node
mpirun -np 2  sander.MPI -O -ng 2 -groupfile group_min_l${X}
mpirun -np 8  sander.MPI -O -ng 2 -groupfile group_equi_l${X}
mpirun -np 8  sander.MPI -O -ng 2 -groupfile group_prod_l${X}
EOF

msub run.pbs.${X}

done
