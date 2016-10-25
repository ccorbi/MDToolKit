###############################################################
## distances.tcl                                                    #
## DESCRIPTION:                                                #
##    This script is quick and easy to provide procedure       #
## for computing the  distance between sidechains of two AA    #
## of Protein                  #
##                                                             #
## EXAMPLE USAGE:                                              #
##         vmd -dispdev text -e  distances.tcl >out.log                                #
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

#set ref_rt [atomselect 1 "backbone and (resid 4 to 7 or resid 25 to 30 or resid 36 to 40 or resid 47 to 50) " frame 0]
#set sel_rt [atomselect 0 "backbone and (resid 4 to 7 or resid 25 to 30 or resid 36 to 40 or resid 47 to 50) "]

set ref [atomselect 1 "backbone " frame 0]
set sel [atomselect 0 "backbone "]

#HELIX - Loop defintion, bb backbone, and full side chain sd

set refbb_hx [atomselect 1  "backbone and resid 1 to 12 and noh"]

set trgbb_hx [atomselect 0  "backbone and resid 1 to 12 and noh"]

set refsd_hx [atomselect 1  "resid 1 to 12 and noh"]

set trgsd_hx [atomselect 0  "resid 1 to 12 and noh"]

#n-SRC defintion, bb backbone, and full side chain sd

set refbb_bh [atomselect 1  "backbone and resid 13 to 34 and noh"]

set trgbb_bh [atomselect 0  "backbone and resid 13 to 34 and noh"]

set refsd_bh [atomselect 1  "resid 1 to 34 and noh"]

set trgsd_bh [atomselect 0  "resid 1 to 34 and noh"]

set refbb_domain [atomselect 1  "resid 35 to 149 and noh"]
set trgbb_domain [atomselect 0  "resid 35 to 149 and noh"]


                    #PHE28 -TH13
                    #PHE28 -PHE10
                    #PHE28 - PHE26
                    #PHE26 - PHE10
                    #PHE26 - PHE14


#Domain
#PHE553
set res28 [atomselect 0  "resid 28 and not backbone and noh"]
#PHE551
set res26 [atomselect 0  "resid 26 and not backbone and noh"]
#VAL557
set res32 [atomselect 0  "resid 32 and not backbone and noh"]


#Helix interesting
#MET531
set res6 [atomselect 0  "resid 6 and not backbone and noh"]
#VAL534
set res9 [atomselect 0  "resid 9 and not backbone and noh"]
#PHE535
set res10 [atomselect 0  "resid 10 and not backbone and noh"]
#THR538
set res13 [atomselect 0  "resid 13 and not backbone and noh"]
#PHE539
set res14 [atomselect 0  "resid 14 and not backbone and noh"]

#set res26 [atomselect 0  "resid 51 and not backbone and noh"]

#outputfiles
set output [open "rmsd.dat" w]
set output_dist [open "distances.dat" w]

puts $output_dist "frame\t28-10\t28-13\t28-26\t26-10\t26-14\t32-6\t32-9"

puts $output "frame\trmsd_glb\trmsd_bbhx\trmsd_sdhx\trmsd_bbbh\trmsd_sdbh\trmsd_bb_domain"
#puts $output "frame\tglobal\tbb_hx\tbb_bh"


#SYSTEM
set all [atomselect 0 all]
#FRAMES
set n [molinfo 0 get numframes]

  for { set i 0 } { $i < $n } { incr i } {
    $sel frame $i
    $all frame $i
    $trgbb_hx frame $i
    $trgsd_hx frame $i
    $trgbb_bh frame $i
    $trgsd_bh frame $i
    $res26 frame $i
    $res28 frame $i
    $res10 frame $i
    $res13 frame $i
    $res14 frame $i
    $res32 frame $i
    $res6 frame $i
    $res9 frame $i

    $trgbb_domain frame $i
    $refbb_domain frame $i

    set matrix [measure fit $sel $ref]
    #$all move [measure fit $sel $ref]
    $all move $matrix
    set rmsd_glb [measure rmsd $sel $ref]

    set matrix [measure fit $trgbb_hx $refbb_hx]
    #$all move [measure fit $sel $ref]
    $all move $matrix
    set rmsd_bbhx [measure rmsd $trgbb_hx $refbb_hx]
    set rmsd_sdhx [measure rmsd $trgsd_hx $refsd_hx]

    set matrix [measure fit $trgbb_bh $refbb_bh]
    #$all move [measure fit $sel $ref]
    $all move $matrix

    set rmsd_bbbh [measure rmsd $trgbb_bh $refbb_bh]
    set rmsd_sdbh [measure rmsd $trgsd_bh $refsd_bh]


    set matrix [measure fit $trgbb_domain $refbb_domain]
    #$all move [measure fit $sel $ref]
    $all move $matrix

    set rmsd_bb_domain [measure rmsd $trgbb_domain $refbb_domain]
    #set rmsd_sdbh [measure rmsd $trgsd_bh $refsd_bh]
    # print the RMSD
    puts $output "$i\t$rmsd_glb\t$rmsd_bbhx\t$rmsd_sdhx\t$rmsd_bbbh\t$rmsd_sdbh\t$rmsd_bb_domain"

    set coord1 [measure center $res28 weight mass]
    set coord2 [measure center $res10 weight mass]
    set d28_10 [vecdist $coord1 $coord2]


    set coord2 [measure center $res13 weight mass]
    set d28_13 [vecdist $coord1 $coord2]


    set coord2 [measure center $res26 weight mass]
    set d28_26 [vecdist $coord1 $coord2]

    set coord1 [measure center $res10 weight mass]

    set d26_10 [vecdist $coord1 $coord2]
    set coord1 [measure center $res14 weight mass]

    set d26_14 [vecdist $coord1 $coord2]
    # print distances

    set coord1 [measure center $res32 weight mass]
    set coord2 [measure center $res6 weight mass]
    set d32_6 [vecdist $coord1 $coord2]

    set coord2 [measure center $res9 weight mass]
    set d32_9 [vecdist $coord1 $coord2]

    puts $output_dist "$i\t$d28_10\t$d28_13\t$d28_26\t$d26_10\t$d26_14\t$d32_6\t$d32_9"
  }


close $output
close $output_dist

                    #PHE28 -TH13
                    #PHE28 -PHE10
                    #PHE28 - PHE26
                    #PHE26 - PHE10
                    #PHE26 - PHE14
                    #vmd -dispdev text -e rmsd.1.areas.tcl
quit
