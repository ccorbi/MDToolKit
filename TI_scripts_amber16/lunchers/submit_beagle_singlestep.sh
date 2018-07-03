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
cat << EOF > run.gpusub
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

$mdrun -O -i 1-min.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 -o min.1.out -inf min.1.info -e min.1.en -r min.1.rst7 -l min.log
$mdrun -O -i 2-min.in -c min.1.rst7  -p ti.parm7 -o min.2.out -inf min.2.info -e min.2.en -r min.2.rst7 -l min.log

$mdrun  -O -i 3-equil.in -o equil.3.out -p ti.parm7 -c min.2.rst7 -r equil.3.rst7 -ref min.2.rst7
$mdrun  -O -i 4-equil.in -o equil.4.out -p ti.parm7 -c  equil.3.rst7 -r equil.4.rst7 -ref equil.3.rst7

$mdrun  -O -i  prod.in -o ti.1.out -p  ti.parm7 -c  equil.4.rst7 -x ti.1.nc  -r ti.1.rst7 -ref equil.4.rst7 -e ti001.en


for N in {2..5}
do
$mdrun -O -i prod.in -o ti.\$N.out -p ti.parm7 -c  ti.\$((\$N-1)).rst7 -x ti.\$N.nc  -r ti.\$N.rst7 -e ti00\$N.en
done

EOF

# adapt above for your job scheduler
qsub run.gpusub

cd ..
done
    cd ..
    done
