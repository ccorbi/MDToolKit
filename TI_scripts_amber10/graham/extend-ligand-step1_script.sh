#!/bin/bash
wd=`pwd`
pre0='x-native-lig'
pre1='x-native-lig'
mask_native=':res'

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
  nstlim = 1000000,	dt = 0.001,ioutfm=1,
  temp0 = 300.0,	ntt = 3,	gamma_ln = 2,
  ntp = 1,	pres0 = 1.0,	taup = 2.0,
  icfe=1,       clambda = 0.${X},
EOF

cp mdin_prod${N}_v0_l${X} mdin_prod${N}_v1_l${X}

cat << EOF >> mdin_prod${N}_v0_l${X}
 &end
EOF
cat << EOF >> mdin_prod${N}_v1_l${X}
  crgmask='${mask_native}',
 &end
EOF

cat << EOF > run.pbs.${X}
#!/bin/sh
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4096M
#SBATCH --time=3-00:00

# consider changing the name below to something more explanatory
#SBATCH --job-name=${X}_${mask_native}_complex-run


# change to the submission directory
cd \$SLURM_SUBMIT_DIR

# module purge
module load  fftw gcc


# Run the parallel version of sander, using all 8 cores in the node
mpirun -np 8  sander.MPI -O -ng 2 -groupfile group_prod${N}_l${X}
EOF

sbatch -A rrg-pmkim run.pbs.${X}

done
