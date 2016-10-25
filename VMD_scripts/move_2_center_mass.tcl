foreach i $argv {
mol new $i type pdb
set sel [atomselect top "protein"]
set com [measure center $sel weight mass]
$sel moveby [vecscale -1.0 $com]
$sel writepdb reset_$i
}
quit

