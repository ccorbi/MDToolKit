foreach i $argv {
mol new $i type pdb
set clean [atomselect 0 "noh"]
$clean writepdb template.pdb
}
quit
