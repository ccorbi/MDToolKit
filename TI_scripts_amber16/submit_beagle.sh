#!/bin/sh
#
windows=$(seq 0.00 0.100 1.00)
amber=/usr/local/amber16
mdrun=${amber}/bin/pmemd.cuda
N='1'

for status in complex ligand; do
  cd $status
  echo $status
  for step in decharge vdw_bonded recharge;do
    cd $step
    echo $step
    for w in $windows; do
      cd $w
      echo $W
      pwd
# adapt below for your job scheduler
cat << EOF > run.gpusub
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q gpus
#$ -N job_$status_$step_$w
#$ -e sge.err
#$ -o sge.out
#$ -l gpu=1,h_rt=12:00:00

# export all environment variables to SGE
#$ -V

export AMBER_PREFIX=/usr/local/amber16
export AMBERHOME=/usr/local/amber16
export PATH="\${AMBER_PREFIX}/bin:\${PATH}"
# Add location of Amber Python modules to default Python search path
if [ -z "\$PYTHONPATH" ]; then
    export PYTHONPATH="\${AMBER_PREFIX}/lib/python2.7/site-packages"
else
    export PYTHONPATH="\${AMBER_PREFIX}/lib/python2.7/site-packages:\${PYTHONPATH}
"
fi
if [ -z "\${LD_LIBRARY_PATH}" ]; then
   export LD_LIBRARY_PATH="\${AMBER_PREFIX}/lib"
else
   export LD_LIBRARY_PATH="\${AMBER_PREFIX}/lib:\${LD_LIBRARY_PATH}"
fi
export CUDA_HOME="/usr/local/cuda"
export LD_LIBRARY_PATH="/usr/local/cuda/lib:\${LD_LIBRARY_PATH}"

device=\`setGPU\`
export CUDA_VISIBLE_DEVICES="\${device}"
echo \$CUDA_VISIBLE_DEVICES


# Run the parallel version of sander, using all 8 cores in the node
#mpirun -np 12 $mdrun -O -i min.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 -o min.out -inf min.info -e mi
#n.en -r min.rst7 -l min.log
#mpirun -lsf $mdrun -i heat.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 \
  #-O -o heat.out -inf heat.info -e heat.en -r heat.rst7 -x heat.nc -l heat.log
$mdrun -O -i min.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 -o min.out -inf min.info -e min.en -r min.rst7 -l min.log

$mdrun -O -i heat.in -c min.rst7 -ref min.rst7 -p ti.parm7 -o heat.out -inf heat.info -e heat.en -r heat.rst7 -x heat.nc -l heat.log

$mdrun -O -i prod.in -c heat.rst7 -p ti.parm7 -o ti00${N}.out -inf ti00${N}.info -e ti00${N}.en -r ti00${N}.inpcrd -x ti00${N}.nc -l ti00${N}.log
EOF

# adapt above for your job scheduler
qsub run.gpusub

cd ..
done
  cd ..
  done
    cd ..
    done
