#!/usr/bin/perl -w

use POSIX;
my @vec=();

my $number=500;
while(<>)
{
    chomp;
    my @temp=split(/\s+/);
    
    push(@vec,$_) if(scalar @temp==2);
}
my @selected=();

for(my $i=0;$i<$number;$i++)
{
    my $rand=floor(rand()*scalar(@vec));
    push(@selected,"$vec[$rand]\n");
    splice(@vec,$rand,1);
   # print $rand."\n";
}
print @selected;
