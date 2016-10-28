## remove.tcl
## DESCRIPTION:
##    Simple script to remove chains from a pdb
##  PDB file or just remove Hydrogens
##
##
## EXAMPLE USAGE:
##         vmd -dispdev text -e  remove_chain.tcl -args  pdbfile chain
##
##
##   AUTHORS:
##  C corbiverge
##
##
################################################################



set pdbfile [lindex $argv 0 ]
set target_chain [lindex $argv 1 ]

mol new $pdbfile type pdb
set clean [atomselect 0 "chain $target_chain"]
#set clean [atomselect 0 "all and not (resid XXXX and sidechain)"]
#set clean [atomselect 0 "noh  and   (not resid  XXX  or backbone)"]
set fbasename [file rootname [file tail $pdbfile]]
$clean writepdb $fbasename-$target_chain.pdb
quit
