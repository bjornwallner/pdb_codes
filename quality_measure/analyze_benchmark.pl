#!/usr/bin/perl -w 




require '/afs/pdc.kth.se/home/b/bjornw/source/perl/bjornlib.pl';

my %Ssum=();
my %MeanS=();
my %mx=();
my %cusr=();

# Read all files
foreach my $file(glob("benchmark/*.dat"))
{
   # print $file."\n";
    next if(-s $file<2000);
    if($file=~/L\_(\d+)\-minsim\_(\d+)\-factor\_([\d\.]+)\./)
    {
	my $L=$1;
	my $minsim=$2;
	my $factor=$3;
#	print "$L $minsim $factor\n";
	#next if($factor!= 5);
	#next if($minsim!=2 && $minsim!=3 && $minsim!=1); # && $minsim!=2 && $minsim!=3);
	#next if($L!=4);
	#
	#print "$L $minsim $factor\n";
	#exit;
	open(FILE,"$file");
	while(<FILE>)
	{
	    m/^Ssum:\s+([\d\.]+).+of\s+([\d\.]+)\s\(MX:\s([\d\.]+)\).+\=\s+([\d\.]+)\sCPU/;
	    #print;
	    chomp;
	    my $Ssum=$1;
	    my $MeanS=$2;
	    my $mx=$3;
	    my $cusr=$4;
	    #	print "$Ssum $MeanS $mx $cusr\n";
	    push(@{$Ssum{$L}{$minsim}{$factor}},$Ssum);
	    #push(@{$MeanS{$L}{$minsim}{$factor}},$MeanS);
	    push(@{$mx{$L}{$minsim}{$factor}},$mx);
	    push(@{$cusr{$L}{$minsim}{$factor}},$cusr);
	}
	close(FILE);
    }
}
#exit;

# Statistics


my %correlation=(); #to L=4,factor=5,minsim=3
my %times=(); 
foreach my $L(keys(%Ssum))
{
    foreach my $minsim(keys(%{$Ssum{$L}}))
    {
	foreach my $factor(keys(%{$Ssum{$L}{$minsim}}))
	{
	    #print "$L $minsim $factor\n-------------------\n";
	    #for(my $i=0;$i<10;$i++)
	    #{
	    #	print "${$Ssum{$L}{$minsim}{$factor}}[$i] ${$Ssum{4}{3}{5}}[$i]\n";
	    #	
	    #}
	    next if(scalar @{$Ssum{$L}{$minsim}{$factor}}!=500);
	    #print "\n";
	    
	    $correlation{$L}{$minsim}{$factor}=corrcoef([@{$Ssum{$L}{$minsim}{$factor}}],[@{$Ssum{4}{3}{5}}]);
	    $times{$L}{$minsim}{$factor}=sum(@{$cusr{$L}{$minsim}{$factor}});
	    
	    printf("%10.7f %8.2f (L=%2d minsim=%3d factor=%5.2f)\n",$correlation{$L}{$minsim}{$factor}, $times{$L}{$minsim}{$factor},$L,$minsim, $factor);
	  #  print "TIME       :  $times{$L}{$minsim}{$factor} ($L $minsim $factor)\n";
	
	}
			   
	
    }




}




