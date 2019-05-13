#Modify the csv file created by Python script
#2014.8.6, Xiaoyan

#!usr/bin/perl -w
use strict;
use Image::Size;

open my $file, '<', "forward_GCAT.csv" or die "Couldn't open csv file: $!";
open my $out, '>', "CSVmodified.csv";

my $firstline = <$file>;
my @first = split'\n', $firstline;
my @Ln1 = split',', $first[0];;
$firstline = join',',(@Ln1, 'Image_PathName_Spec_blob','Image_FileName_Spec_blob','Dimension_X','Dimension_Y');
print $out "$firstline\n";


my $refimagepath;
my $refimagefile;
my $refsizex;
my $refsizey;
my $match;
my $specblob;
my @mismatch;
my @line;

while (my $fileline = <$file>) {
	my @li = split'\n', $fileline;
	@line = split',', $li[0];;
	my @hyb = split'hyb', $line[3];
	my @blobpath = split'base', $line[12];
	my @blobpath2 = split'_c2', $blobpath[2];
	$specblob = join'',($blobpath[0], 'base', $hyb[1]+1, '/base', $hyb[1]+1, '_c2', $blobpath2[1]);
	#print("$specblob\n");
	my $image = join'',($specblob, $line[13]);
	(my $sizex,my $sizey) = imgsize($image);
	if ($hyb[0] == 0) {
		$refsizex = $sizex;
		$refsizey = $sizey;
	}
	else {
		if ($sizex == $refsizex && $sizey == $refsizey) {}
		else {
			push (@mismatch,$line[0]);
		}	
	}
	my $newline = join',', (@line, $specblob, $line[13], $sizex, $sizey);
	print $out "$newline\n";
}

if (@mismatch) {
	print "Image dimension mismatch detected.\n";
	open my $out2, '>', "DimensionMismatches.txt";
	foreach my $mis (@mismatch) {
		print $out2 "$mis\t";
	}
	close $out2;
}

close ($file);
close ($out);