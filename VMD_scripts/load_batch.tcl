###############################################################
## load_batches.tcl                                            #
## DESCRIPTION:                                                #
##    This script is quick and easy to visulize and compare    #
## a dataset of strucutres                                      #
##                #
##                                                             #
## EXAMPLE USAGE:                                              #
## vmd -dispdev  -e  load_batches.tcl -args pdb1.pdb pdb2.pdb                             #
##                        #
##                                                             #
##   AUTHORS:                                                  #
##  C corbiverge                                              #
##                                                             #
##                                                              #
################################################################
#




foreach i $argv {
  mol new $i type pdb
#NewCartoon by Chain
mol delrep 0 top
mol representation NewCartoon
mol color chain
mol selection {all}
mol material Opaque
mol addrep top

#Licorice Mutation point
mol color ColorID 1
mol representation licorice
mol selection {resid 185 and noh}
mol addrep top

#Target protein surface
mol color ResType
mol representation MSMS
mol selection {chain A }
mol addrep top

#Peptide Surface
#mol color ColorID 1
#mol representation MSMS
#mol selection {chain G}
#mol addrep top
mol rename top $i
}
