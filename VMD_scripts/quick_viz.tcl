#LOAD TRAJ

set tra [glob *mdcrd]
set trajs [lsort -ascii  $tra]
set pram [glob *prm]

set trajs [glob *md*mdcrd]
set pram [glob *prm]

mol addrep 0
mol new $pram type {parm7} first 0 last -1 step 1 waitfor all

foreach trj $trajs {
	mol addfile $trj type crdbox first 0 last -1 step 500 waitfor all 0
}

