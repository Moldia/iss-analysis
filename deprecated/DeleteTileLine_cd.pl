#this is used to delete information regarding specific tiles in .csv file
#the script will run in the current working directory

#!usr/bin/perl -w
use strict;

#here to specify which positions want to delete
#use .. to indicate range, eg 1..5
my @delete = (6, 11, 26, 31, 46, 51, 66, 71, 86, 91, 106, 111, 126, 131, 146, 151, 166, 171, 186, 191, 206, 211, 226, 231, 246, 251, 266, 271, 286, 291);
my $deleteline = 0;
#return the number of valuables in the array 
my $deletenumber = scalar @delete;

open my $file, '<', "CPfile.txt" or die "Couldn't open CPfile.txt: $!";
open my $out, '>', "CPfile_Selected.txt";

#read $file line by line
my $firstline = <$file>;
print $out $firstline;

while (my $line = <$file>) {
	my @Readline = split 'hyb', $line;
	
	#smart match function of perl 5.10.1, verifies if a value exists in the array or not.
	if ($Readline[0] ~~ @delete) {
		$deleteline ++;
	}
	else {
		print $out $line;
	}
}

print "\n$deletenumber tiles to delete.\n$deleteline of lines deleted.\nNew CPfile created.\n";

closedir (D);
close ($file);
close ($out);