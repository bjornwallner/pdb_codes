#!/usr/bin/perl -w
use Benchmark;

my $scopdir="/afs/pdc.kth.se/home/a/arnee/structpredict/SCOP/scop-1.63/pdb/";
my $model_path="/afs/pdc.kth.se/home/b/bjornw/modelling/ProQhires/models.fixed/prof-prof/";
my $tmpdir="/scratch/lgscore/";
mkdir($tmpdir) if(!-e $tmpdir);
`copy_to_scratch`;




my $lgscore="./lgscore";

my $file="selected_models";
#my @L=(4,6,8,10);
#
#my @factor=(5,2.24);
#
#my @minsim=(3,4,9,16,25,36,49,64,81,100,121,144,169,196,225,256,289);


my @minsim=(36,49,64,81,100,121,144,169,196,225,256,289,324,361,400);

my @L=(4,8,6,10,12);
my @factor=(1,1.5,2,2.23);
#my $L=$ARGV[0];
#my $minsim=$ARGV[1];
#my $factor=$ARGV[2];


#@L=(4);
#@minsim=(3);
#@factor=(5);

chdir($tmpdir);
foreach my $L(@L)
#for(my $L=4;$L<=4;$L++)
{
#    my $L=4;
    foreach my $minsim(@minsim)
    {
	foreach my $factor(@factor)
	{
	    #print "$L $minsim $factor\n";
	   # next;
	    my $outfile="/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/benchmark4/L_$L-minsim_$minsim-factor_$factor.dat";
	    my $local_outfile="L_$L-minsim_$minsim-factor_$factor.dat";
	    #print $outfile."\n";
	    #next;
	    #`rm $outfile`;
	    #next;
	    next if(-e $outfile);
	    `touch $outfile`;

	    #exit;
	    my $start_time = new Benchmark;
	    my $out=`$lgscore -list set1 -L $L -minsim $minsim -factor $factor`; #|egrep 'Ssum|SCORE'`;
	    my $end_time = new Benchmark;
	    my $difference = timediff($end_time, $start_time);
	    my $timstr=timestr($difference);
		
	    open(OUT,">$local_outfile");
	    print OUT "$out";
	    print OUT "TIME: $timstr\n";
	    close(OUT);
	    `mv $local_outfile $outfile`;
	}
    }
}
