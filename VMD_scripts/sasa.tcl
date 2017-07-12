## sasa.tcl
## DESCRIPTION:
##    This script is quick and easy to provide procedure
## for computing the Solvent Accessible Surface Area (SASA)
## of Protein and allows Users to select regions of protein.
##
## EXAMPLE USAGE:
##         vmd -dispdev text -e  sasa.tcl -args protein >out.log
##         Selection: chain A, resid 4 to 5 and protein
##
##   AUTHORS:
##
##
##
################################################################
#

#LOAD AMBER TRAJ
# load all files
set tra [glob *mdcrd]
set trajs [lsort -ascii  $tra]
set pram [glob *prm]

mol addrep 0
mol new $pram type {parm7} first 0 last -1 step 1 waitfor all

foreach trj $trajs {
	mol addfile $trj type crdbox first 0 last -1 step 1 waitfor all 0
}

puts -nonewline "\n \t \t Selection: "
#gets stdin selmode
set selmode [lindex $argv 0 ]
## selection
set sel [atomselect top "$selmode"]
set protein [atomselect top "protein"]
set n [molinfo top get numframes]
set output [open "SASA_$selmode.dat" w]
# sasa calculation loop
for {set i 0} {$i < $n} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $protein -restrict $sel]
	puts "\t \t progress: $i/$n"
	puts $output "$sasa"
}
puts "\t \t progress: $n/$n"
puts "Done."
puts "output file: SASA_$selmode.dat"
close $output
#
quit
