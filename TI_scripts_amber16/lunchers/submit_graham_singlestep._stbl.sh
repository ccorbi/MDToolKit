#!/bin/sh


# adapt below for each User or enviorment
amber=/home/ccorbi/amber16
mdrun=${amber}/bin/pmemd.cuda

# by default 10 lambdas
# Adjust the number of lambdas if it is need it
windows=$(seq 0.00 0.100 1.0)

for status in unfolded folded; do
  cd $status
  echo $status
    for w in $windows; do
      cd $w
      echo $W
      pwd


# adapt below for your job scheduler
export LD_LIBRARY_PATH=$AMBERHOME/lib:$LD_LIBRARY_PATH

cat << EOF > run.gpusub
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1              # request GPU "generic resource"
#SBATCH --mem=8g               # memory per node
#SBATCH --time=0-5:00            # time (DD-HH:MM)
#SBATCH --output=%N-%j.out        # %N for node name, %j for jobID
# nvidia-smi

# change to the submission directory
cd \$SLURM_SUBMIT_DIR

# module purge
module load gcc/4.8.5 openmpi fftw

export PATH=/usr/local/cuda-8.0/bin:\$PATH
export LPATH=/usr/lib64/nvidia-current:\$LPATH
export LPATH=/usr/lib64/nvidia:\$LPATH
export LIBRARY_PATH=/usr/lib64/nvidia-current:\$LIBRARY_PATH
export LIBRARY_PATH=/usr/lib64/nvidia:\$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64/nvidia-current:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64/nvidia:/usr/local/cuda-8.0/lib64:/usr/local/cuda-8.0/lib:\$LD_LIBRARY_PATH
export CUDA_HOME=/usr/local/cuda-8.0

# load amber
source ${amber}/amber.sh


$mdrun -O -i 1-min.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 -o min.1.out -inf min.1.info -e min.1.en -r min.1.rst7 -l min.log
$mdrun -O -i 2-min.in -c min.1.rst7  -p ti.parm7 -o min.2.out -inf min.2.info -e min.2.en -r min.2.rst7 -l min.log
$mdrun  -O -i 3-equil.in -o equil.3.out -p ti.parm7 -c min.2.rst7 -r equil.3.rst7 -ref min.2.rst7
$mdrun  -O -i 4-equil.in -o equil.4.out -p ti.parm7 -c  equil.3.rst7 -r equil.4.rst7 -ref equil.3.rst7


$mdrun  -O -i  prod.in -o ti.1.out -p  ti.parm7 -c  equil.4.rst7 -x ti.1.nc  -r ti.1.rst7 -ref equil.4.rst7 -e ti001.en


for N in {2..3}
do
$mdrun -O -i prod.in -o ti.\$N.out -p ti.parm7 -c  ti.\$((\$N-1)).rst7 -x ti.\$N.nc  -r ti.\$N.rst7 -e ti00\$N.en
done
EOF

# adapt above for your job scheduler & User
sbatch -A def-pmkim run.gpusub

  cd ..
  done
cd ..
done
