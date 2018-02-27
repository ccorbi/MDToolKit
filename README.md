# MD ToolKit

Collection of scripts to lunch and use analysis of amber Molecular Dynamics.

*NOTE: All scripts are currently in early development and are not guaranteed to function properly. Some knowledge of shell scripting, Linux systems and, Amber package will be of use. 

# Dependencies

  - [Amber](http://ambermd.org/) 
  - [VMD](http://www.ks.uiuc.edu/Research/vmd/)
  - Python. I personally recommend [Anaconda](https://anaconda.org/anaconda/python) distribution and package managment. 

# Python Requirements

  - Numpy
  - Pandas
  - seaborn
  - scipy

# Installation

To use these scripts, simply copy them into a desired directory and add that directory to your PATH. It is recommended to use these scripts on a Linux system.


---

# Thermodynamic Integrations Sctipts

Thermodynamic integration (TI) is a method used to compare the difference in free energy between two given states (e.g., A and B) whose potential energies U A {\displaystyle U_{A}} U_{A} and U B {\displaystyle U_{B}} U_{B} have different dependences on the spatial coordinates. Because the free energy of a system is not simply a function of the phase space coordinates of the system, but is instead a function of the Boltzmann-weighted integral over phase space (i.e. partition function), the free energy difference between two states cannot be calculated directly. In thermodynamic integration, the free energy difference is calculated by defining a thermodynamic path between the states and integrating over ensemble-averaged enthalpy changes along the path. Such paths can either be real chemical processes or alchemical processes. An example alchemical process is the Kirkwood's coupling parameter method.

## Usage

* Scripts Amber 10

  **DEPRICATED**. Scripts to setup, lunch and analyse TI  using Amber 12 or earlier.

* Scripts Amber 16

Scripts to setup, lunch and analyze TI using Amber 16 or later. There are 2 different protocols, single step and multiple steps. Multiple step protocols run independently for wildtype discharge, vdw+bonded, and mutation recharge.  Three-step is more resource consuming, however, it trend to provide more accurate dG for mutations where a charge is involved. Single topology approach is used for both protocols. These scripts intend to provide a generic solution, under certain circumstances custom tuning must need it. 

setup_ti_single script will write all inputs for amber run a one-step TI, and setup_ti_multiple will do it for three-steps. The only requirements in both scripts are to provide a PDB file with the target, a separate PDB with the ligand, and the mutation. Amber renumber files assigning 1 to the first amino acid in the PDB, please use this numbering to indicate the script the mutation. 

```bash

setup_ti_single.py --target targt.pdb --ligand ligre_ext.pdb --mutation PRO12ALA

setup_ti_multi.py --target targ_init.pdb --ligand lig_init.pdb --mutation ASP1THR

```

Options are common in both scripts. 

```bash
usage: setup_ti_single.py [-h] [--ligand LIGAND] [--chain CHAIN]
                          [--mutation MUTATION] [--target TARGET]
                          [--increment INCREMENT]
                          [--scmask-ignore SCMASK_IGNORE]
                          [--decouple-mask DECOUPLE_MASK] [--psteps PSTEPS]
                          [--ouput-folder OUTPUT_FOLDER]

optional arguments:
  -h, --help            show this help message and exit
  --ligand LIGAND
  --chain CHAIN
  --mutation MUTATION
  --target TARGET
  --increment INCREMENT
  --scmask-ignore SCMASK_IGNORE
  --decouple-mask DECOUPLE_MASK
  --psteps PSTEPS
  --ouput-folder OUTPUT_FOLDER


```

* ANALYSIS


Scripts to preprocess Amber TI output. Both scripts will return plots, lambdas averages, and final deltaG. Recommend to skip the first 500 steps (1 ns if integration step is set 2fs) to improve accuracy, however, this may be adjusted for noise trajectories.  

```bash

analysis_multistep.py --skip 500

analysis_onestep.py --skip 500
```

Options are common in both scripts. 


```bash
usage: analysis_onestep.py [-h] [--no-extrapolation]
                           [--ignore [IGNORE [IGNORE ...]]] [--skip SKIP]
                           [--limit LIMIT] [--no-plots] [--no-rms]

optional arguments:
  -h, --help            show this help message and exit
  --no-extrapolation    by default if lambda 0 or 1 is missing the value is
                        extrapolated. this argument disables this behaibour
  --ignore [IGNORE [IGNORE ...]]
                        Ignore a lambda or lambas, format -> step-lambda step-lambda2 e.g. complex-0.100 
  --skip SKIP           skip # of step to calculate the average on each
                        lambda, 500 or 1ns recommended value, default = 0
  --limit LIMIT         limit on the # of step to calculate the average dvdl
                        for each lambda, default until the end
  --no-plots            do not generate plots, only the ddg
  --no-rms              do not generate RMSD and RMSF data

```


* LUNCHERS

- for beagle: submit_beagle.sh
- for graham submit_graham.sh



# VMD Scripts

Common md analysis, rmsd, rmsf, SASA, etc.. 

# Amber Scripts

Generic Run
