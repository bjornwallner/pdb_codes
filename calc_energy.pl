#!/usr/bin/perl -w 

$prosa="/afs/pdc.kth.se/home/b/bjornw/bjorn/modules/ProsaII/bin/prosaII -f";
#$working_dir=`pwd`;
#chomp($working_dir);
#$working_dir.="/";
#$MODEL_DIR="/afs/pdc.kth.se/home/b/bjornw/bjorn/modelling/MODELLER_FILES/models/";
#$MODEL_DIR="";
#$outfile="z-scores";
#`touch $outfile`;
#$done="done";
#$host=`hostname`;
#chomp($host);
$model=$ARGV[0];
if(!-e "$model.slp" || -s "$model.slp"<250)
{
    $cmd_file=$model.".cmd";
    #print $cmd_file."\n";
    open(TMP,">$cmd_file");
    #print TMP "pdb_dir = $MODEL_DIR\n";
    print TMP "read pdb $model obj\n";
    print TMP "init zscore pII3.0.long.ply\ncombine type sdev\n";
    print TMP "zscore obj $model\n";
    close(TMP);
    print "Calculating energy for $model ...\n";
    
    `$prosa $cmd_file`;
    
    #print "grep ^[0-9] $model.slp";
    #$scores=`grep ^[0-9] $model.slp`;
    #chomp($scores);
    #$score=substr($scores,10,190);
    #print $scores."\n";
    #print $score."\n";
    #open(OUT,">>$outfile");
    #printf OUT ("%-40s %-190s\n",$model,$score);
    #close(OUT);
    #`rm $model.sor $model.cmd`;
}


sub check_done
{
    my $words=shift;
    #print $file."\n";
    my $match=1;
    open(DONE,"$done") || die "Cannot open $done.\n";
    while(<DONE>)
    {
	#chomp;
	#$line=$_;
	if(/$words/)
	{
	    #print "YES";
	    $match=0;
	    last;
	}
    }
    close(DONE);
    #print $match."\n";
    $match;
}
sub module 
{
  eval `$ENV{MODULESHOME}/bin/modulecmd perl @_`;
}
sub source
{
    print @_;
    print "\n";
    eval `. @_`;
}








