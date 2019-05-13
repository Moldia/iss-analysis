#this is used to delete and rename Nuclei.png and Nuclei.png
#the files will be renamed based on the original tile number resulted from re-tiling
#the script will run in the current working directory

#!usr/bin/perl -w
use strict;

my @listTile = <Tile_*.png>; 
my @listNuclei = <Nuclei_*.png>; 

#write the original tile numbers that should not be included here
my @del = ();

open my $rename, '>', "FilesRenamed_OriginalNum.txt";

print $rename "Files renamed according to the original tile numbers:\n";

#-------------------------------------------------------------
my $countTileRe = 0;
my $countTileFile = 0;
my %hashTile;
my @listTile_sorted;

foreach my $fileTile (@listTile) {

	my @Tilename = split 'Tile_', $fileTile;
	my @Tilenr = split '.png', $Tilename[1];
	my $numberTile = @Tilenr[0];
	
	$hashTile {$fileTile} = $numberTile;

}

foreach my $Tile (sort {$hashTile{$a} <=> $hashTile{$b}} keys %hashTile) {
		push (@listTile_sorted, $Tile);
    }
	
foreach my $sortedTile (@listTile_sorted) {

	$countTileRe ++;
	while ($countTileRe ~~ @del) {
		$countTileRe ++;
	}
	my $nameTile = join '', ('TileRe_', $countTileRe, '.png');
	rename $sortedTile, $nameTile;
	$countTileFile ++;
	print $rename "$sortedTile to $nameTile\n";
}

print "$countTileFile of Till_ files renamed.\n";

#-------------------------------------------------------------
my $countNucleiRe = 0;
my $countNucleiFile = 0;
my %hashNuclei;
my @listNuclei_sorted;

foreach my $fileNuclei (@listNuclei) {

	my @Nucleiname = split 'Nuclei_', $fileNuclei;
	my @Nucleinr = split '.png', $Nucleiname[1];
	my $numberNuclei = @Nucleinr[0];
	
	$hashNuclei {$fileNuclei} = $numberNuclei;

}

foreach my $Nuclei (sort {$hashNuclei{$a} <=> $hashNuclei{$b}} keys %hashNuclei) {
		push (@listNuclei_sorted, $Nuclei);
    }
	
foreach my $sortedNuclei (@listNuclei_sorted) {

	$countNucleiRe ++;
	while ($countNucleiRe ~~ @del) {
		$countNucleiRe ++;
	}
	my $nameNuclei = join '', ('NucleiRe_', $countNucleiRe, '.png');
	rename $sortedNuclei, $nameNuclei;
	$countNucleiFile ++;
	print $rename "$sortedNuclei to $nameNuclei\n";
}
print "$countNucleiFile of Nuclei_ files renamed.\n";

#-------------------------------------------------------------
close ($rename);
