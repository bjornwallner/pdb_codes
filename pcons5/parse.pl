#!/usr/bin/perl -w 

@para=`grep PARA read_this4`;
@score=`grep SCORE read_this4`;
#print scalar @score," ",scalar @para,"\n";
print "<pre>";
for($i=0;$i<scalar @score;$i++)
{
    chomp($score[$i]);
    chomp($para[$i]);
    
   
    @temp=split(/\s+/,$score[$i]);
    @temp3=split(/\//,$temp[1]);
    $name=$temp3[$#temp3];
	
    @temp2=split(/\s+/,$para[$i]);
    $str=join(" ",@temp[2..7]);
    print "<pre><a href=$name.pdb>$name.pdb</a> $str ";
    print "$temp2[$#temp2-1]\n";
    #$str=join(" ",@score[1..6]);

    #exit;
}
