#!/usr/bin/perl -w

$infile=$ARGV[0];


@hcut=(3.0,3.2,3.4,3.6,3.8,4.0);
@hangle=(1.0,1.2,1.4,1.6,1.8,2.0);

@hcut=(3.6);
@hangle=(1.2);

foreach $hcut(@hcut)
{
    foreach $hangle(@hangle)
    {
	my @stride_h=();
	my @dssp_h=();
	my @both_h=();
	
	my @stride_e=();
	my @dssp_e=();
	my @both_e=();
	
	my @stride_c=();
	my @dssp_c=();
	my @both_c=();
	
	my @stride=();
	my @dssp=();
	my @both=();
	open(IN,"$infile");
	while(<IN>)
	{

	    chomp;
	    $pdb=$_;
	    print $pdb."\n";  
	    $outfile="out/cuts-$hcut-$hangle";
	    ($my_seq,$my_ss)=`./cont_surf $pdb $hcut $hangle|egrep '^SEQ|^SS'|awk '{print \$2}'`;
	    ($stride_seq,$stride_ss)=`/afs/pdc.kth.se/home/a/arnee/MODULES/scientific/stride/linux/bin/stride $pdb|stride_parser.pl`;
	    print "/afs/pdc.kth.se/home/a/arnee/MODULES/scientific/stride/linux/bin/stride $pdb|stride_parser.pl\n";
	    ($dssp_seq,$dssp_ss)=`/afs/pdc.kth.se/home/a/arnee/MODULES/scientific/dssp/linux/bin/dsspcmbi $pdb|dssp_parser.pl`;
	    
	    chomp($my_seq);
	    chomp($my_ss);
	    chomp($stride_seq);
	    chomp($stride_ss);
	    chomp($dssp_seq);
	    chomp($dssp_ss);
	    
	    $len=length($my_ss);
	    
	    if($len == length($stride_ss) && $len == length($dssp_ss))
	    {
		##print "STRIDE: ";
		($frac_stride1,$frac_stride2,$frac_stride3,$frac_stride4)=`compare_str.pl $my_ss $stride_ss`;
		#print $out."\n";
		##print "DSSP: ";
		($frac_dssp1,$frac_dssp2,$frac_dssp3,$frac_dssp4)=`compare_str.pl $my_ss $dssp_ss`;
		#print $out."\n";
		print "compare_str.pl $my_ss $dssp_ss long\n";
		##print "STRIDE_DSSP: ";
		($frac_both1,$frac_both2,$frac_both3,$frac_both4)=`compare_str.pl $stride_ss $dssp_ss`;
		#print "compare_str.pl $stride_ss $dssp_ss\n";
		#print $out;
		#print "\n";
		if($frac_stride1 ne '' && $frac_dssp1 ne '' && $frac_both1 ne '')
		{
		    chomp($frac_stride1);
		    chomp($frac_stride2);
		    chomp($frac_stride3);
		    chomp($frac_stride4);
		    chomp($frac_dssp1);
		    chomp($frac_dssp2);
		    chomp($frac_dssp3);
		    chomp($frac_dssp4);
		    chomp($frac_both1);
		    chomp($frac_both2);
		    chomp($frac_both3);
		    chomp($frac_both4);
		    
		    #print $frac_both1."\n\n";
		    push(@stride_h,$frac_stride1);
		    push(@dssp_h,$frac_dssp1);
		    push(@both_h,$frac_both1);
		    
		    push(@stride_e,$frac_stride2);
		    push(@dssp_e,$frac_dssp2);
		    push(@both_e,$frac_both2);
		    
		    push(@stride_c,$frac_stride3);
		    push(@dssp_c,$frac_dssp3);
		    push(@both_c,$frac_both3);
		    
		    push(@stride,$frac_stride4);
		    push(@dssp,$frac_dssp4);
		    push(@both,$frac_both4);
		}
	    }
	    my $stride_h=mean(@stride_h);
	    my $dssp_h=mean(@dssp_h);
	    my $both_h=mean(@both_h);
	    
	    my $stride_e=mean(@stride_e);
	    my $dssp_e=mean(@dssp_e);
	    my $both_e=mean(@both_e);
	    
	    my $stride_c=mean(@stride_c);
	    my $dssp_c=mean(@dssp_c);
	    my $both_c=mean(@both_c);
	    
	    
	    my $stride=mean(@stride);
	    my $dssp=mean(@dssp);
	    my $both=mean(@both);
	    
	    #open(OUT,">$outfile");
	    printf ("$outfile %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n",$stride_h,$dssp_h,$both_h,$stride_e,$dssp_e,$both_e,$stride_c,$dssp_c,$both_c,$stride,$dssp,$both);
	}
    }
}


sub mean
{
    my @vec=@_;
    my $mean=0;
    my $len=scalar @vec;
    if($len != 0)
    {
	foreach my $val(@vec)
	{
	    $mean+=$val;
	}
	$mean=$mean/$len;
    }
    else
    {
	return 0;
    }
    return $mean;


}
