# usage:
# # # vmd -e check_aa.tcl  -f loop.pdb -args position
# # #
# # #
# # #
# # #
# #
# # # Parse command line
#
set position [lindex $argv 0 ]
#
#
#
#
# #NewCartoon by Chain
mol delrep 0 top
mol representation NewCartoon
mol color chain
mol selection {all}
mol material Opaque
mol addrep top
#

# #Licorice Mutation point
mol color Element
mol representation licorice
set seltext "{resid $position}"
#set sel [atomselect 0 $seltext]
mol selection $seltext
mol addrep top

# #Target protein surface
#
