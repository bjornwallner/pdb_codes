#!/usr/bin/perl -w


print "\{";
while(<>)
{
    print"\{";
    @temp=split(/\s+/);
    foreach $number(@temp)
    {
	print "$number,"; 
    }
    print "\},\n";
}
