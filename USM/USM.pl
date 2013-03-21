#!/usr/bin/perl -w

$pdb1=$ARGV[0];
$pdb2=$ARGV[1];

$tmpfile="/tmp/KCS.temp.$$.";
$tmpfile1=$tmpfile."1";
$tmpfile2=$tmpfile."2";
$tmpfile12=$tmpfile."12";
$tmpfile21=$tmpfile."21";
#generate_contact_maps
#this file produce $tmpfile1 $tmpfile2 $tmpfile3 ($tmpfile1+$tmpfile2)
`/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/USM/contact_map $pdb1 $pdb2 $tmpfile 2 6.5`;

`cat $tmpfile1 $tmpfile2 > $tmpfile12`;
`cat $tmpfile2 $tmpfile1 > $tmpfile21`;

$str=`gzip -c $tmpfile1`;
$bytes1=length($str);
$str=`gzip -c $tmpfile2`;
$bytes2=length($str);

$str=`gzip -c $tmpfile12`;
$bytes12=length($str);
$str=`gzip -c $tmpfile21`;
$bytes21=length($str);


#print "compute_KCS_complexity($bytes1,$bytes2,$bytes12,$bytes21)\n";
$dist=compute_KCS_complexity($bytes1,$bytes2,$bytes12,$bytes21);
print "DISTANCE: $dist\n";
`rm $tmpfile1 $tmpfile2 $tmpfile12 $tmpfile21`;


sub compute_KCS_complexity
{
    my ($x,$y,$xy,$yx)=@_;
    my $maxNum=1;
    my $maxDem=1;
    if( ($yx-$y)>($xy-$x) )
    {
	$maxNum = $yx-$y;
    }
    else
    {
	$maxNum = $xy-$x;
    }
    
    if($x>$y)
    {
	$maxDen = $x;
    }
    else
    {
	$maxDen = $y;
    }
    return ($maxNum/$maxDen);
    
}
