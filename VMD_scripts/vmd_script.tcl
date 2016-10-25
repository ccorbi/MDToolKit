
proc fitfram { molref selref moltrg seltrg } {
  set ref [atomselect $molref $selref]
  set sel [atomselect $moltrg $seltrg]
  set all [atomselect $moltrg all]
  set n [molinfo $moltrg get numframes]
   
  for { set i 0 } { $i < $n } { incr i } {
    $sel frame $i
    $all frame $i
    $all move [measure fit $sel $ref]
  }
  return
}

#vmd  -dispdev text -eofexit < rot.tcl > output.log


#set trajs [glob *MD*mdcrd]
#set pram [glob *MD*prmtop]
#mol addrep 0
#mol new $pram type {parm7} first 0 last -1 step 1 waitfor all 
#mol addfile $trajs type crdbox first 0 last -1 step 1 waitfor all 0
#set origen [glob dry*pdb]
#mol addrep 1
#mol new $origen type pdb

#fitfram 1 "protein noh" 0 "protein noh"
#set a get numframes
#echo $a
#animate write crdbox {aligned.mdcrd} beg 0 end -1 0


set trajs [glob *mdcrd]
set pram [glob *prm]
mol addrep 0
mol new $pram type {parm7} first 0 last -1 step 1 waitfor all
foreach trj $trajs {
	mol addfile $trj type crdbox first 0 last -1 step 1 waitfor all 0
}
set sel [atomselect top "resid 10 36 49 51 54"]
set protein [atomselect top "protein"]
set n [molinfo top get numframes]
set output [open "SASA_pocket.dat" w]
# sasa calculation loop
for {set i 0} {$i < $n} {incr i} {
        molinfo top set frame $i
        set sasa [measure sasa 1.4 $protein -restrict $sel]
        puts "\t \t progress: $i/$n"
        puts $output "$sasa"
}
puts "\t \t progress: $n/$n"
puts "Done."
puts "output file: SASA_pocket.dat"
close $output

