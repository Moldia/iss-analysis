#this is used to delete information regarding specific tiles in .csv file
#the script will run in the current working directory

#!usr/bin/perl -w
use strict;

#here to specify which positions want to delete
#use .. to indicate range, eg 1..5
my @delete = (6, 11, 26, 31, 46, 51, 66, 71, 86, 91, 106, 111, 126, 131, 146, 151, 166, 171, 186, 191, 206, 211, 226, 231, 246, 251, 266, 271, 286, 291);

#here to specify the reference layer to not to delete
my $Ref = 01;

my $deleteline = 0;
#return the number of valuables in the array 
my $deletenumber = scalar @delete;
my $Refline = 0;

open my $file, '<', "CPfile.txt" or die "Couldn't open CPfile.txt: $!";
open my $out, '>', "CPfile_Selected_RefSaved.txt";

#read $file line by line
my $firstline = <$file>;
print $out $firstline;

while (my $line = <$file>) {
	my @Readline = split 'hyb', $line;
	my @Hybnum = split 'base', $Readline[1];
	
	#smart match function of perl 5.10.1, verifies if a value exists in the array or not. magic!
	if ($Readline[0] ~~ @delete && $Ref != $Hybnum[0]) {
		$deleteline ++;
	}
	elsif ($Readline[0] ~~ @delete && $Ref == $Hybnum[0]) {
		$Refline ++;
		print $out $line;
	}
	else {
		print $out $line;
	}
}

print "\n$deletenumber tiles to delete.\n$deleteline of lines deleted.\n$Refline of lines regarding reference layer saved.\nNew CPfile created.\n";

closedir (D);
close ($file);
close ($out);