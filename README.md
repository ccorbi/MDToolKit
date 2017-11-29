# COMMON Tools for MD Analysis

## Thermodynamic Integrations

Thermodynamic integration (TI) is a method used to compare the difference in free energy between two given states (e.g., A and B) whose potential energies U A {\displaystyle U_{A}} U_{A} and U B {\displaystyle U_{B}} U_{B} have different dependences on the spatial coordinates. Because the free energy of a system is not simply a function of the phase space coordinates of the system, but is instead a function of the Boltzmann-weighted integral over phase space (i.e. partition function), the free energy difference between two states cannot be calculated directly. In thermodynamic integration, the free energy difference is calculated by defining a thermodynamic path between the states and integrating over ensemble-averaged enthalpy changes along the path. Such paths can either be real chemical processes or alchemical processes. An example alchemical process is the Kirkwood's coupling parameter method.


### Scripts Amber 10

**DEPRICATED**. Scripts to setup, lunch and analyse TI. 

### Scripts Amber 16

Scripts to setup , lunch and analyse TI.

* SETUP ONE STEP

```bash

setup_ti_single.py --scmask-ignore CA,C,O,N,HA,H1,H2,H3 --target targt.pdb --ligand ligre_ext.pdb --mutation PRO12ALA

```


* SETUP MULTIPLE STEP

```bash

setup_ti_single.py --scmask-ignore CA,C,O,N,HA,H1,H2,H3 --target targt.pdb --ligand ligre_ext.pdb --mutation PRO12ALA

```



* SETUP MULTIPLE STEP

```bash

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

* LUNCH MULTIPLE STEP

- for beagle: submit_beagle.sh
- for graham submit_graham.sh

* ANALYSIS

```bash

analysis_multistep.py --skip 500

analysis_onestep.py --skip 500
```

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
