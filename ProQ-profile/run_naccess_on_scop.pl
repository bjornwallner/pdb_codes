#!/usr/bin/perl -w 

$infile=$ARGV[0];
$scopdir="/afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.63/pdb/";

open(IN,$infile);
chdir('/scratch/bjornw/');
`rm pdc.*` if(-e "pdc.rsa");
while(<IN>)
{
    chomp;
    $name=$_;
    $pdbfile=$scopdir.$name.".pdb";
    next if(!-e $pdbfile);
    $outfile=$name.".rsa";
    @temp=split(/\//,$outfile);
    $outfile=$temp[1];
    next if(-e $outfile);
    $exec_string="/afs/pdc.kth.se/home/b/bjornw/modules/naccess $pdbfile";
    print $exec_string."\n";
    $out=`$exec_string`;
    print $out."\n";
    `mv pdc.rsa $outfile`;
    `rm pdc.*`;
   # exit;
    


}
