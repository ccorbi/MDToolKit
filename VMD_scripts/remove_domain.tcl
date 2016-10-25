###############################################################
## remove.tcl                                                    #
## DESCRIPTION:                                                #
##    This script is quick and easy to remove parts from a   #
##  PDB file or just remove Hydrogens                          #
##   #
##                                                             #
## EXAMPLE USAGE:                                              #
##         vmd -dispdev text -e  remove.tcl -args             #
##                                                             #
##                                                             #
##   AUTHORS:                                                  #
##  C corbiverge                                              #
##                                                             #
##                                                              #
################################################################
#



foreach i $argv {
mol new $i type pdb
set clean [atomselect 0 "chain XXXX"]
#set clean [atomselect 0 "all and not (resid XXXX and sidechain)"]
#set clean [atomselect 0 "noh  and   (not resid  XXX  or backbone)"]
$clean writepdb peptide_ncaa.pdb
}
quit
