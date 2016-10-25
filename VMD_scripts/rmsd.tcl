###############################################################
## rmsd.tcl                                                    #
## DESCRIPTION:                                                #
##    This script is quick and easy to provide procedure       #
## for computing the root mean sqrt distance (RMSA)    #
## of Protein and allows Users to select regions of protein.   #
##                                                             #
## EXAMPLE USAGE:                                              #
##         vmd -dispdev text -e  rmsd.tcl >out.log             #
##         Selection: chain A and resid 1                      #
##                                                             #
##   AUTHORS:                                                  #
##  C corbiverge                                              #
##                                                             #
##                                                              #
################################################################
#


#LOAD TRAJ



set tra [glob *mdcrd]
set trajs [lsort -ascii  $tra]
set pram [glob *prm]

mol addrep 0
mol new $pram type {parm7} first 0 last -1 step 1 waitfor all

foreach trj $trajs {
    mol addfile $trj type crdbox first 0 last -1 step 1 waitfor all 0
}


#LOAD REFERENCE

set origen [glob amber*pdb]
mol addrep 1
mol new $origen type pdb

#align to the core
#set ref [atomselect 1 "(backbone and resid 4 to 7) or  (backbone and  resid 25 to 30) or  (backbone and resid 36 to 40) or (backbone and resid 47 to 50) " frame 0]
#set sel [atomselect 0 "(backbone and resid 4 to 7) or  (backbone and resid 25 to 30) or (backbone and resid 36 to 40) or (backbone and resid 47 to 50) "]
set ref [atomselect 1 "backbone" frame 0]
set sel [atomselect 0 "backbone"]

#TOTAL PROTEIN
set ref_total [atomselect 1  "backbone"]
set trg_total [atomselect 0  "backbone"]

#BINDING INTERFACE
#set ref_area8 [atomselect 1  "protein and resid 8 and noh"]
#set trg_area8 [atomselect 0  "protein and resid 8 and noh"]



set output_core [open "rmsd_resid_core.md1.dat" w]
puts $output_core "frame total"

set output_binding [open "rmsd_resid_binding.md1.dat" w]
puts $output_binding "#frame   #total  #8  #10  #13  #17  #33  #35  #36  #49  #51  #53  #54  "
#SYSTEM
set all [atomselect 0 all]
#FRAMES
set n [molinfo 0 get numframes]

  for { set i 0 } { $i < $n } { incr i } {
    $sel frame $i
    $all frame $i
    $trg_total frame $i

    set matrix [measure fit $sel $ref]
    #$all move [measure fit $sel $ref]
    $all move $matrix
    set rmsd [measure rmsd $trg_total $ref_total]


                # print the RMSD
                puts $output_core "$i $rmsd"
#       puts $output_core "$i $rmsd $rmsd18 $rmsd20 $rmsd26 $rmsd28 $rmsd39 $rmsd50 $rmsd55"
  }


close $output_binding
close $output_core

quit