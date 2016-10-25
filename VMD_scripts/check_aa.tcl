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
mol color Element 
mol representation licorice
mol selection {resid 176 }
mol addrep top
}
#Target protein surface
