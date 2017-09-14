#!/bin/sh
#

windows=$(seq 0.00 0.100 1.0)
amber=/home/ccorbi/amber16
mdrun=${amber}/bin/pmemd.cuda
N='1'
#cd complex
for status in complex ligand; do
  cd $status
  for w in $windows; do
    cd $w

  # adapt below for your job scheduler
export LD_LIBRARY_PATH=$AMBERHOME/lib:$LD_LIBRARY_PATH

cat << EOF > run.gpusub
#!/bin/bash
#SBATCH --account=def-mikeuoft
#SBATCH --gres=gpu:1              # request GPU "generic resource"
#SBATCH --mem=4000M               # memory per node
#SBATCH --time=0-10:00            # time (DD-HH:MM)
#SBATCH --output=%N-%j.out        # %N for node name, %j for jobID
# nvidia-smi

# change to the submission directory
cd \$SLURM_SUBMIT_DIR

# module purge
module load gcc/4.8.5 openmpi fftw

export PATH=/usr/local/cuda-7.5/bin:\$PATH
export LPATH=/usr/lib64/nvidia-current:\$LPATH
export LPATH=/usr/lib64/nvidia:\$LPATH
export LIBRARY_PATH=/usr/lib64/nvidia-current:\$LIBRARY_PATH
export LIBRARY_PATH=/usr/lib64/nvidia:\$LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64/nvidia-current:/usr/local/cuda/lib64:/usr/local/cuda/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64/nvidia:/usr/local/cuda/lib64:/usr/local/cuda/lib:\$LD_LIBRARY_PATH
export CUDA_HOME=/usr/local/cuda-7.5
export LD_LIBRARY_PATH="/usr/local/cuda-7.5/lib:\${LD_LIBRARY_PATH}"
# load amber
source ${amber}/amber.sh

# Run the parallel version of sander, using all 8 cores in the node
$mdrun -O -i min.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 -o min.out -inf min.info -e min.en -r min.rst7 -l min.log

$mdrun -O -i heat.in -c min.rst7 -ref ti.rst7 -p ti.parm7 -o heat.out -inf heat.info -e heat.en -r heat.rst7 -x heat.nc -l heat.log

$mdrun -O -i ti.in -c heat.rst7 -p ti.parm7 -o ti00${N}.out -inf ti00${N}.info -e ti00${N}.en -r ti00${N}.inpcrd -x ti00${N}.nc -l ti00${N}.log
EOF

# adapt above for your job scheduler
sbatch -A def-pmkim run.gpusub

  cd ..
  done
cd ..
done
