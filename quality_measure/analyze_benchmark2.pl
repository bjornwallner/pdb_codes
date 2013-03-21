#!/usr/bin/perl -w 




require '/afs/pdc.kth.se/home/b/bjornw/source/perl/bjornlib.pl';

my %Ssum=();
my %LGscore=();
my %MeanS=();
my %mx=();
my %cusr=();

# Read all files
foreach my $file(glob("benchmark3/*.dat"))
{
   # print $file."\n";
    next if(-s $file<2000);
    if($file=~/L\_(\d+)\-minsim\_(\d+)\-factor\_([\d\.]+)[\.\-]/)
    {
	my $L=$1;
	my $minsim=$2;
	my $factor=$3;
	my $step="1";
	if($file=~/step(\d+)/)
	{
	    $step=$1;
	}

	#print "$L $minsim $factor $step\n";
	#if(not(($factor == 5 && $L==4 && $minsim==3)))
	#{
	#    next if($factor!= 2.23);
	#    #next if($minsim!=2 && $minsim!=3 && $minsim!=1); # && $minsim!=2 && $minsim!=3);
	next if($L!=4); # && $L!=8);
	#}

	#
	#print "$L $minsim $factor\n";
	#exit;
	#$outfile=$file."a";
#	open(OUT,">$outfile");
	#open(FILE,"gunzip -c $file|");
	open(FILE,"$file"); #gunzip -c $file|");
	while(<FILE>)
	{

	    
	    if(/^Ssum:\s+([\d\.]+).+of\s+[\d\.]+.+\s+([\d\.]+)\s+\(MX:\s([\d\.]+)\).+\=\s+([\d\.]+)\sCPU/ ||
	       /^SCORE:\s+([\d\.]+)\s[\d\.]+\s([\d\.]+).+\(MX:\s([\d\.]+)\).+\=\s+([\d\.]+)\sCPU/)
	    {
	#	print;
		chomp;
		my $Ssum=$1;
		my $LGscore=$2;
		my $mx=$3;
		my $cusr=$4;
		#print ;
		#print "\n";
		#next if($mx<0.8);
		
#	    print OUT "$Ssum $LGscore $mx $cusr\n";
		#print "$Ssum $LGscore $mx $cusr\n";
		push(@{$Ssum{$L}{$minsim}{$factor}{$step}},$Ssum);
		push(@{$LGscore{$L}{$minsim}{$factor}{$step}},$LGscore);
		push(@{$mx{$L}{$minsim}{$factor}{$step}},$mx);
		push(@{$cusr{$L}{$minsim}{$factor}{$step}},$cusr);
	    }
	    elsif(/^SCORE:\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)/)
	    {
		#print;
		chomp;
		my $Ssum=$1;
		my $LGscore=$3;
		#print ;
		#print "\n";
		#next if($mx<0.8);
		
#	    print OUT "$Ssum $LGscore $mx $cusr\n";
		#print "$Ssum $LGscore\n";
		push(@{$Ssum{$L}{$minsim}{$factor}{$step}},$Ssum);
		push(@{$LGscore{$L}{$minsim}{$factor}{$step}},$LGscore);
	    }
	    elsif(/User Time\:\s+([\d\.]+)/)
	    {
		#print;
		#print $1."\n";
		push(@{$cusr{$L}{$minsim}{$factor}{$step}},$1/100);
	    }
	}
	close(FILE);
#	close(OUT);
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
	    foreach my $step(keys(%{$Ssum{$L}{$minsim}{$factor}}))
	    {
	    #print "$L $minsim $factor\n-------------------\n";
	    #for(my $i=0;$i<10;$i++)
	    #{
	    #	print "${$Ssum{$L}{$minsim}{$factor}{$step}}[$i] ${$Ssum{4}{3}{5}}[$i]\n";
	    #	
	    #}
		if(not(defined(@{$Ssum{4}{3}{5}{1}})))
		{
		    print "4,3,5 not defined\n";
		    exit;
		}
		#print scalar @{$LGscore{$L}{$minsim}{$factor}{$step}};
		#print "\n";
		next if(scalar @{$LGscore{$L}{$minsim}{$factor}{$step}}!=500); #scalar @{$LGscore{4}{3}{5}});
		#  print "$L $minsim $factor\n-------------------\n";
		#  print "\n";
		#   next;
		$correlation{$L}{$minsim}{$factor}{$step}{Ssum}=corrcoef([@{$Ssum{$L}{$minsim}{$factor}{$step}}],[@{$Ssum{4}{3}{5}{1}}]);
		$correlation{$L}{$minsim}{$factor}{$step}{LGscore}=corrcoef([@{$LGscore{$L}{$minsim}{$factor}{$step}}],[@{$LGscore{4}{3}{5}{1}}]);
		$times{$L}{$minsim}{$factor}{$step}=sum(@{$cusr{$L}{$minsim}{$factor}{$step}});
		
		printf("%10.7f %10.7f %8.2f (L=%2d minsim=%3d factor=%s step=%2d)\n",$correlation{$L}{$minsim}{$factor}{$step}{Ssum}, $correlation{$L}{$minsim}{$factor}{$step}{LGscore},$times{$L}{$minsim}{$factor}{$step},$L,$minsim, $factor,$step);
		###exit;
	    }
	  #  print "TIME       :  $times{$L}{$minsim}{$factor}{$step} ($L $minsim $factor)\n";
	
	}
			   
	
    }




}




