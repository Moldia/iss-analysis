#!usr/bin/perl -w
use strict;

open my $out,'>',"symbollist.csv";

my @list = <TransparentMarkers/*.png>;
foreach my $file (@list) {
	my @fullpath = split('/',$file);
	my @fullname = split(/\./,$fullpath[1]);
	print $out "$fullname[0]\n";
}
 
print $out "\n";
close($out);
