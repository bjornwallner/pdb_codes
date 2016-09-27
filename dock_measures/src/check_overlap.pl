#!/usr/bin/perl -w
$n=0;
while(<>) {
    @temp=split(/\s+/);
    if(/MODEL/) {
	$model{"$temp[1] $temp[2]"}=1
    }
    if(/NATIVE/) {
	$native{"$temp[1] $temp[2]"}=1;
	$n++;
	    
    }
    
}
$c=0;
foreach $key(keys(%native)) {

    if(defined($model{$key})) {
	$c++;
	print "PERL: Overlap $key\n";
    }
}
$frac=$c/$n;
print "$c $n $frac\n";
