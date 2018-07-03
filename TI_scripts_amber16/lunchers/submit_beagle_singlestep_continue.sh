#!/bin/sh
#
# adapt below for each User or enviorment
amber=/usr/local/amber16
mdrun=${amber}/bin/pmemd.cuda

# by default 10 lambdas
# Adjust the number of lambdas if it is need it
windows=$(seq 0.00 0.100 1.0)


for status in complex ligand; do
  cd $status
  echo $status

    for w in $windows; do
      cd $w
      echo $W
      pwd
# adapt below for your job scheduler
cat << EOF > con.gpusub
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q gpus
#$ -N job_$status_$w
#$ -e sge.err
#$ -o sge.out
#$ -l gpu=1

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

export TEMPDIR='/tmp/'



for N in {2..5}
do
$mdrun -O -i prod.in -o ti.\$N.out -p ti.parm7 -c  ti.\$((\$N-1)).rst7 -x ti.\$N.nc  -r ti.\$N.rst7 -e ti00\$N.en
done


EOF

# adapt above for your job scheduler
qsub con.gpusub

cd ..
done
    cd ..
    done
