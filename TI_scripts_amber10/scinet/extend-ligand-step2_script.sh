#!/bin/bash
wd=`pwd`
pre0='x-native-lig'
pre1='name-x-lig'
mask_native=':res'
mask_noncanon=':res'

N=$1

for X in 1 2 3 4 5 6 7 8 9

do
cat << EOF > group_prod${N}_l${X}
-O -i mdin_prod${N}_v0_l${X} -o ${pre0}_prod${N}_v0_l${X}.out -p ${pre0}.prm -c ${pre0}_prod$(($N-1))_v0_l${X}.rst -r ${pre0}_prod${N}_v0_l${X}.rst -x ${pre0}_prod${N}_v0_l${X}.crd -e ti_${N}.en
-O -i mdin_prod${N}_v1_l${X} -o ${pre1}_prod${N}_v1_l${X}.out -p ${pre1}.prm -c ${pre1}_prod$(($N-1))_v1_l${X}.rst -r ${pre1}_prod${N}_v1_l${X}.rst -x ${pre1}_prod${N}_v1_l${X}.crd -e ti_${N}.en
EOF


cat << EOF > mdin_prod${N}_v0_l${X}
NPT production
 &cntrl
  imin = 0,	ntx = 5,	irest = 1,
  ntpr = 1000,	ntwr = 100000,	ntwx = 10000, ntwe=1000,
  ntf = 1,	ntc = 1,
  ntb = 2,	cut = 9.0,
  nstlim = 750000,	dt = 0.001,ioutfm=1,
  temp0 = 300.0,	ntt = 3,	gamma_ln = 2,
  ntp = 1,	pres0 = 1.0,	taup = 2.0,
  icfe=1,       clambda = 0.${X},
EOF

cp mdin_prod${N}_v0_l${X} mdin_prod${N}_v1_l${X}

cat << EOF >> mdin_prod${N}_v0_l${X}
  ifsc=1,
  crgmask='${mask_native}',
  scmask='${mask_native}',
 &end
EOF
cat << EOF >> mdin_prod${N}_v1_l${X}
  ifsc=1,
  crgmask='${mask_noncanon}',
  scmask='${mask_noncanon}',
 &end
EOF

cat << EOF > run.pbs.${X}

#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00

# consider changing the name below to something more explanatory
#PBS -N ${X}_ligand-run

# PMK RAPI
#PBS -A cvj-675-aa

# change to the submission directory
cd \$PBS_O_WORKDIR

module purge
module load intel/15.0.2 intelmpi/5.0.3.048 gcc/4.8.1  cuda/6.0 amber/14.0


# Run the parallel version of sander, using all 8 cores in the node
mpirun -np 8  sander.MPI -O -ng 2 -groupfile group_prod${N}_l${X}
EOF

qsub run.pbs.${X}

done
