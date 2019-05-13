#delete and rename Tile.png, Nuclei_outlines.png and Cell_outlines.png according to the original tile number
#Xiaoyan, 2014-03-04

#!usr/bin/perl -w
use strict;

my @listTile = <*Tile*.png>; 
my @listNuclei = <*Nuclei*.png>;
my @listCell = <*Cell*.png>;

#write the original tile numbers that should not be included here
my @del = ();

#two output files, one for deleted files and the other for renamed files
open my $del, '>', "FilesDeleted.txt";
open my $rename, '>', "FilesRenamed.txt";

print $del "Files deleted:\n";
print $rename "Files renamed:\n";

#-------------------------------------------------------------
my $countTile = 0;
my $countTileRe = 0;
my $countTileFile = 0;
my %hashTile;
my @listTile_sorted;

foreach my $fileTile (@listTile) {
	
	my @Tilename = split 'Tile', $fileTile;
	my @Tilenr = split '.png', $Tilename[1];
	my $numberTile = @Tilenr[0];
	
	#delete the files from base 2-4
	if (($numberTile-1) % 4 != 0) {
		unlink $fileTile;
		$countTile ++;
		#print "Delete file: $FileTile\n";
		print $del "$fileTile\n";
	}
	else {
		$hashTile {$fileTile} = $numberTile;
	}
}

print "$countTile of Tile files deleted.\n";

foreach my $Tile (sort {$hashTile{$a} <=> $hashTile{$b}} keys %hashTile) {
		push (@listTile_sorted, $Tile);
    }
	
foreach my $sortedTile (@listTile_sorted) {
	$countTileRe ++;
	while ($countTileRe ~~ @del) {
		$countTileRe ++;
	}
	my $nameTile = join '', ('Tile_', $countTileRe, '.png');
	rename $sortedTile, $nameTile;
	$countTileFile ++;
	print $rename "$sortedTile to $nameTile\n";
	#print "Renamed file $sortedTile To $nameTile\n";
}

print "$countTileFile of Tile files renamed.\n";

#-------------------------------------------------------------
my $countNuclei = 0;
my $countNucleiRe = 0;
my $countNucleiFile = 0;
my %hashNuclei;
my @listNuclei_sorted;

foreach my $fileNuclei (@listNuclei) {
	
	my @Nucleiname = split 'Nuclei_outlines', $fileNuclei;
	my @Nucleinr = split '.png', $Nucleiname[1];
	my $numberNuclei = @Nucleinr[0];
	
	#delete the files from base 2-4
	if (($numberNuclei-1) % 4 != 0) {
		unlink $fileNuclei;
		$countNuclei ++;
		#print "Delete file: $FileNuclei\n";
		print $del "$fileNuclei\n";
	}
	else {
		$hashNuclei {$fileNuclei} = $numberNuclei;
	}
}

print "$countNuclei of Nuclei files deleted.\n";

foreach my $Nuclei (sort {$hashNuclei{$a} <=> $hashNuclei{$b}} keys %hashNuclei) {
		push (@listNuclei_sorted, $Nuclei);
    }
	
foreach my $sortedNuclei (@listNuclei_sorted) {
	$countNucleiRe ++;
	while ($countNucleiRe ~~ @del) {
		$countNucleiRe ++;
	}
	my $nameNuclei = join '', ('Nuclei_', $countNucleiRe, '.png');
	rename $sortedNuclei, $nameNuclei;
	$countNucleiFile ++;
	print $rename "$sortedNuclei to $nameNuclei\n"
	#print "Renamed file $sortedNuclei To $nameNuclei\n";
}

print "$countNucleiFile of Nuclei files renamed.\n";

#-------------------------------------------------------------
my $countCell = 0;
my $countCellRe = 0;
my $countCellFile = 0;
my %hashCell;
my @listCell_sorted;

foreach my $fileCell (@listCell) {
	
	my @Cellname = split 'Cell_outlines', $fileCell;
	my @Cellnr = split '.png', $Cellname[1];
	my $numberCell = @Cellnr[0];
	
	#delete the files from base 2-4
	if (($numberCell-1) % 4 != 0) {
		unlink $fileCell;
		$countCell ++;
		#print "Delete file: $FileCell\n";
		print $del "$fileCell\n";
	}
	else {
		$hashCell {$fileCell} = $numberCell;
	}
}

print "$countCell of Cell files deleted.\n";

foreach my $Cell (sort {$hashCell{$a} <=> $hashCell{$b}} keys %hashCell) {
		push (@listCell_sorted, $Cell);
    }
	
foreach my $sortedCell (@listCell_sorted) {
	$countCellRe ++;
	while ($countCellRe ~~ @del) {
		$countCellRe ++;
	}
	my $nameCell = join '', ('Cell_', $countCellRe, '.png');
	rename $sortedCell, $nameCell;
	$countCellFile ++;
	print $rename "$sortedCell to $nameCell\n";
	#print "Renamed file $sortedCell To $nameCell\n";
}

print "$countCellFile of Cell files renamed.\n";

#-------------------------------------------------------------

close ($del);
close ($rename);