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


#my @minsim=(36,49,64,81,100,121,144,169,196,225,256,289,324,361,400,441,484,529,576,625);

#my @L=(4,8,6,10,12,14,16,18,20,22,24,26);
#my @factor=(0.25,0.5,0.75,1,1.5,2,2.23);
#my $L=$ARGV[0];
#my $minsim=$ARGV[1];
#my $factor=$ARGV[2];

my @minsim=(121);
my @factor=(0.500000);
my @L=(4);
my @step=(2); #,3,4,5,6);


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
	    foreach my $step(@step)
	    {
		#print "$L $minsim $factor\n";
	   # next;

		my $outfile="/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/benchmark3/L_$L-minsim_$minsim-factor_$factor.dat";
		my $local_outfile="L_$L-minsim_$minsim-factor_$factor.dat";
		if($step>1)
		{
		    $outfile="/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/benchmark3/L_$L-minsim_$minsim-factor_$factor-step$step.dat";
		    $local_outfile="L_$L-minsim_$minsim-factor_$factor-step$step.dat";
		}

		#print $outfile."\n";
		#next;
		#`rm $outfile`;
		#next;
		next if(-e $outfile);
		#print $outfile."\n";
		#next;
		#print "$L $minsim $factor\n";
		#next;
		`touch $outfile`;
		
		open(OUT,">$local_outfile");
		open(FILE,$file);
		my $counter=0;
		while(<FILE>)
		{
		    $counter++;
		    #	print $counter."\n";
		    chomp;
		    my ($model,$mx)=split(/\s+/);
		    my @temp=split(/\_on\_/,$model);
		    my $target_id=$temp[0];
		    my $subdir=substr($model,8,1);
		    #	my $correct_file="$scopdir$subdir/$target_id.pdb";
		    my $correct_file="$target_id.pdb";
		    #print $correct_file."\n";
		    #next;
		    #my $modelfile=$model_path.$model;
		    my $modelfile=$model;
#		print $modelfile."\n";
#		next;
		    #		my $out=`$lgscore $modelfile $correct_file -L $L -minsim $minsim -factor $factor|egrep 'Ssum|SCORE'`;
		    #print "$lgscore $modelfile $correct_file\n";
		    #exit;
		    my $start_time = new Benchmark;
		    #	print "$modelfile $correct_file\n";
		    #	next;
		    my $out=`$lgscore -L $L -minsim $minsim -factor $factor $modelfile $correct_file`; #|egrep 'Ssum|SCORE'`;
		    my $end_time = new Benchmark;
		    $out=~s/\n/ /g;
		    #	print $out;
		    #	exit;
		    #  my $out=`./lgscore_old $modelfile $correct_file -L $L -minsim $minsim -factor $factor|grep Ssum`;
		    
		    #my @line=grep(@out,/Ssum/);
		    #print @line,"\n";
		    chomp($out);
		    my $difference = timediff($end_time, $start_time);
		    my $timstr=timestr($difference);
		    
		    print OUT "$out (MX: $mx) $timstr\n";
		}
		#    exit;
		close(FILE);
		close(OUT);
		`cp $local_outfile $outfile`;
	    }
	}
    }
}
