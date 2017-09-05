#!/usr/bin/perl -w

open FILE, $ARGV[0];
if ($#ARGV>0) {
	$quiet = $ARGV[1];
} else {
	$quiet = 0;
}
$i=0;

while ($line=<FILE>) {
	$line =~ /\s*([0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)/;
	($lambda[$i],$dvdl[$i],$err[$i])=($1,$2,$3);
# uncomment for attractive
#	if ( $dvdl[$i]<0 ) { $dvdl[$i]=0 ;}
# or repulsive part only
#	if ( $dvdl[$i]>0 ) { $dvdl[$i]=0 ;}
	$i++;
}

# interpolate to l=0
if ($lambda[0] != 0) {
	$ynew = ($lambda[0]*$dvdl[1] - $lambda[1]*$dvdl[0]) / ($lambda[0]-$lambda[1]);
	unshift @lambda,0;
	unshift @dvdl,$ynew;
	unshift @err,$err[0];
#	print "Added new value at $lambda[0] of $dvdl[0] +- $err[0]\n";
	$i++;
}
#and to l=1
if ($lambda[$i-1] != 1.0) {
	$ynew = ( $dvdl[$i-1]*($lambda[$i-2] - 1 ) + $dvdl[$i-2]*(1-$lambda[$i-1]) ) / ($lambda[$i-2]-$lambda[$i-1]);
	push @lambda,1;
	push @dvdl,$ynew;
	push @err, $err[$i-1];
	$i++;
#	print "Added new value at $lambda[$i-1] of $dvdl[$i-1] +- $err[$i-1]\n";
}
for ($j=0 ; $j<$i ; $j++) {
	if ($i == 1) {
		$width[$j] = 1;
	} elsif ($j==0) {
		$width[$j] = 0.5 * ( $lambda[$j] + $lambda[$j+1] );
	} elsif ($j == $i-1) {
		$width[$j] = 1 - 0.5 * ( $lambda[$j] + $lambda[$j-1] );
	} else {
		$width[$j]= 0.5 * ( $lambda[$j+1] - $lambda[$j-1] );
	}
}

$tot_dvdl=0;
$tot_err=0;
for ($j=0 ; $j<$i ; $j++) {
	$tot_dvdl+=$width[$j]*$dvdl[$j];
#	Simple additive error propagation
	$tot_err+=$width[$j]*$err[$j];
#	Gaussian Error Propagation
#	$tot_err+=($width[$j]*$err[$j])**2;
	if ($quiet==0) {
		printf "l: %4.2f (w=%4.3f)\t cont: %4.4f+-%4.4f", $lambda[$j], $width[$j], $width[$j]*$dvdl[$j], $width[$j]*$err[$j];
#		print "l: $lambda[$j] (w=$width[$j])\t cont:",$width[$j]*$dvdl[$j],"+-",$width[$j]*$err[$j];
		print "\n";
#		printf "%4.4f  ", $width[$j]*$dvdl[$j];
	}
}
# For Gaussian error proagation
# $tot_err = sqrt($tot_err);
if ($quiet==0) {
	print "Total DV/DL:", $tot_dvdl, "+-",$tot_err,"\n";
} else {
#	print $tot_dvdl,"\t",$tot_err,"\n";
	printf ("%4.2f +- %2.2f\n", $tot_dvdl, $tot_err);
}

