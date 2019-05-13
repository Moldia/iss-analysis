#rename blob_align.png and blob_prealign.png files

#!usr/bin/perl -w
use strict;

my @listAlign = <*blobs_align*.png>; 
my @listPrealign = <*blobs_prealign*.png>;

my @del = ();

#one output file for renamed files
open my $rename, '>', "FilesRenamed.txt";

print $rename "Files renamed:\n";

#-------------------------------------------------------------
my $countAlignRe;
my $countAlignFile = 0;

foreach my $fileAlign (@listAlign) {
	
	my @Alignname = split 'blobs_align', $fileAlign;
	my @Alignnr = split '.png', $Alignname[1];
	my $numberAlign = @Alignnr[0];
	
	#delete the files from base 2-4
	for (my $i =1; $i<5; $i ++) {
		if (($numberAlign-$i) % 4 == 0) {
			$countAlignRe = ($numberAlign-$i)/4+1;
			if ($cuontAlignRe ~~ @del) {
				my $index  = 0;
				$index ++
				until $del[$index] == $countAlignRe;
				$countAlignRe = $countAlignRe + $index;
				}
			my $nameAlign = join '', ('blobs_b', $i, '_', $countAlignRe, '.png');
			rename $fileAlign, $nameAlign;
			$countAlignFile ++;
			print $rename "$fileAlign to $nameAlign\n";
		}
	}
}

print "$countAlignFile of blob_align files renamed.\n";

#-------------------------------------------------------------
my $countPrealignRe;
my $countPrealignFile = 0;

foreach my $filePrealign (@listPrealign) {
	
	my @Prealignname = split 'blobs_prealign', $filePrealign;
	my @Prealignnr = split '.png', $Prealignname[1];
	my $numberPrealign = @Prealignnr[0];
	
	#delete the files from base 2-4
	for (my $i =1; $i<5; $i ++) {
		if (($numberPrealign-$i) % 4 == 0) {
			$countPrealignRe = ($numberPrealign-$i)/4+1;
			if ($countPrealignRe ~~ @del) {
				my $index  = 0;
				$index ++
				until $del[$index] == $countPrealignRe;
				$countPrealignRe = $countPrealignRe + $index;
				}
			my $namePrealign = join '', ('blobs_prealign_b', $i, '_', $countPrealignRe, '.png');
			rename $filePrealign, $namePrealign;
			$countPrealignFile ++;
			print $rename "$filePrealign to $namePrealign\n";
		}
	}
}

print "$countPrealignFile of blob_prealign files renamed.\n";

#-------------------------------------------------------------
close ($rename);