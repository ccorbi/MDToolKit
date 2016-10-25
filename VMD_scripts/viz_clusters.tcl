#LOAD Clusters



set clusters [glob *pdb]


foreach trj $clusters {
	mol new $trj type {pdb} first 0 last -1 step 1 waitfor all
  mol delrep 0 top
  mol representation NewCartoon
  mol color chain
  mol selection {all}
  mol material Opaque
  mol addrep top
}
