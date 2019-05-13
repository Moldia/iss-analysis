#Create csv file for CellProfiler as input file
#To run the script on windows, copy the file to the folder, from command window change directory, and type perl filename.pl
#2014.8.6, Xiaoyan

#!usr/bin/perl -w
use strict;
use Image::Size;

#create needed variables to store data
my $sizex;
my $sizey;
my @refsizeX;
my @refsizeY;
my @mismatchx;
my @mismatchy;

#variables required in CP file
#each line defined as array @god
#initial values in @god are the names of each column
my @god = ('HybStep', 'Image_PathName_T', 'Image_FileName_T', 'Image_PathName_G', 'Image_FileName_G', 'Image_PathName_C', 'Image_FileName_C', 'Image_PathName_A', 'Image_FileName_A', 'Image_PathName_Nuclei', 'Image_FileName_Nuclei', 'Image_PathName_General_blob', 'Image_FileName_General_blob', 'Image_PathName_Spec_blob', 'Image_FileName_Spec_blob');
my @godD = ('Metadata_position', 'HybStep', 'Image_FileName_General_blob', 'Size_Reference_X', 'Size_Reference_Y', 'Image_FileName_Spec', 'Size_Current_X', 'Size_Current_Y', 'X_match', 'Y_match');

my @file;
my @refc1;
my @refc2;
my $refpathname;
my $hybstep = 0;
my $count;
my $printline = 0;
#$foldnr to define base* to read
#if doing forward sequencing, initial value of $foldnr should be 1
#if doing reverse sequencing, initial value of $foldnr should be 4
my $foldnr = 1;

#read the contents of the directory using glob function
my @list = <*>; 
#grep !/^\.\.?$/ reads all contents except '.' and '..'
#sort @list into alphabetic order using binary "cmp"
@list = sort {$a cmp $b} @list;

print "Contents of the current directory:\n";
foreach my $file (@list) {
	print "$file\n";
}
print "\n";

#create a file handle $out for writing
#the output file will be created in the input directory
open my $out, '>', "CPinput.csv";
print $out "Metadata_position";
foreach (@god) {
	print $out ",$_";
}
print $out "\n";

open my $out2, '>', "CPinput_ImageDimensions.txt";
foreach (@godD) {
	print $out2 "$_\t";
}
print $out2 "\n";

#the folders containing images for each sequencing cycle are required to name as base1, base2, ...
#no space between text and number
foreach my $folder (@list) {
	my @foldername = ('base', $foldnr);
	my $foldername = join '', @foldername;
	#scalar $count to count the number of .tif files read in a folder
	$count=0;
	#$pos to represent the matadata position
	my $pos = 0;
	#only read folder named as $foldername
	
	if ($folder ~~ $foldername) {
		#assume when start   reading another folder, starts next hybridization step
		$hybstep ++;
		#if doing reverse sequencing, $foldnr --
		$foldnr ++;
		print "Reading directory: $folder\n";
		#read all tif images in the folder
		my @listb = <"$folder/*.tif">;
		@listb = sort {$a cmp $b} @listb;
		
		if (@listb) {
			#each file name together with its path is read into scalar $imagefile
			
			foreach my $imagefile (@listb) {
				$count ++;
				#separate the path from the filename
				my @imagename = split '/', $imagefile;
				#image file name ($imagename[1]) is the original file name in the folder
				#each file name contains base number and scene and channel number
				#@id contains base* and s*c* 
				#@subid contains scene number and channel number
				my @id = split '_', $imagename[1];
				my @subid = split 'c', $id[1];				
				my $scene = $subid[0];
				my $channel = $subid[1];
				$file[$channel-1] = $imagename[1];
				
				#when every file for one scene is read output the relative information to file
				#in this case, output file is written when every sixth file is read, because there are in total 6 files for one scene (six channels)
				if ($count%6 == 0) {
					$pos ++;
					
					if ($hybstep == 1) {
						push (@refc1, $file[0]);
						push (@refc2, $file[1]);
						$refpathname = join'', ($folder, "\/");
					}
					
					if ($channel == 2) {
					#read the dimensions of the image DO_c2
					($sizex, $sizey) = imgsize ($imagefile);
					}
					#contents of HybStep column must be in the format of hyb01, hyb02, etc in CP file
					my @hybstep = ('hyb0',$hybstep);
					my $hyb = join '', @hybstep;
					#the pathname requires a slash in the end
					my @pathname = ($folder,"\/");
					my $pathname = join'', @pathname;
					#change here if the channel setting is not the same
					@god = ($hyb, $pathname, $file[2], $pathname, $file[3], $pathname, $file[4], $pathname, $file[5], $refpathname, $refc1[$pos-1], $refpathname, $refc2[$pos-1], $pathname, $file[1]); 
					
					#store the position of tiles that have different size from reference image
					my $matchx;
					my $matchy;
					if ($refsizeX[$pos-1] == $sizex) {
						$matchx = 1;
					}
					else { 
						$matchx = 0;
						push (@mismatchx, $pos);
					}
					
					if ($refsizeY[$pos-1] == $sizey) {
						$matchy = 1;
					}
					else { 
						$matchy = 0;
						push (@mismatchy, $pos);
					}
					
					@godD = ($pos, $hyb, $refc2[$pos-1], $refsizeX[$pos-1], $refsizeY[$pos-1], $file[1], $sizex, $sizey, $matchx, $matchy);
					
					print $out "$pos";
					
					foreach my $god (@god) {
						print $out ",$god";
						$printline ++;
					}
					
					foreach my $godD (@godD) {
						print $out2 "$godD\t";
					}
					
					print $out "\n";
					print $out2 "\n";
				}
			}
		}
		else {
			print "Empty directory $folder!\n";
		}
	}
}

#sort the mismatch matrices in numeric order
@mismatchx = sort {$a <=> $b} @mismatchx;
@mismatchy = sort {$a <=> $b} @mismatchy;

print $out2 "\nX dimension mismatches (if any):\n";
foreach my $misx (@mismatchx) {
	print $out2 "$misx\t";
}

print $out2 "\nY dimension mismatches (if any):\n";
foreach my $misy (@mismatchy) {
	print $out2 "$misy\t";
}
if ($hybstep && $printline) {
	print "Successfully generated CPinputs.csv. \nMay need to be resaved as CSV (MS-DOS) format in excel.\n";
}
elsif ($hybstep == 0) {
	print "No folders containing base images!\n";
	print $out "No folders containing base images!\n";
}

close ($out);
close ($out2);