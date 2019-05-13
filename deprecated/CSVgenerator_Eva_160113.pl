#Create csv file for CellProfiler as input file
#To run the script on windows, copy the file to the folder, from command window change directory, and type perl filename.pl
#2014.8.6, Xiaoyan

#!usr/bin/perl -w
use strict;
#use Image::Size;

#create needed variables to store data
my $sizex;
my $sizey;
#my @refsizeX;
#my @refsizeY;
#my @mismatchx;
#my @mismatchy;

#variables required in CP file
#each line defined as array @god
#initial values in @god are the names of each column

my @file;
my @refblob;
my $refpathname;
my $count;
my $hybstep = 0;
my $printline = 0;
#$tilenr to define base* to read
#if doing forward sequencing, initial value of $tilenr should be 1
#if doing reverse sequencing, initial value of $tilenr should be 4
my $tilenr = 1;

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
my @god = ('start_X_pos', 'start_Y_pos','HybStep','Image_PathName_T', 'Image_FileName_T', 'Image_PathName_G', 'Image_FileName_G', 'Image_PathName_C', 'Image_FileName_C', 'Image_PathName_A', 'Image_FileName_A', 'Image_PathName_General_blob', 'Image_FileName_General_blob', 'Image_PathName_Spec_blob', 'Image_FileName_Spec_blob');
foreach (@god) {
	print $out ",$_";
}
print $out "\n";

#the folders containing images for each sequencing cycle are required to name as base1, base2, ...
#no space between text and number
foreach my $tilefolder (@list) {
	my @tilename = ('tile', $tilenr);
	my $tilename = join '', @tilename;
	my @filelist = <"$tilename/*">;

	$hybstep = 0;
	foreach my $folder (@filelist) {
		my @foldername = ('base', $hybstep);
		my $foldername = join '', @foldername;

		#only read folder named as $foldername
		if ($folder ne $foldername) {
			#assume when start   reading another folder, starts next hybridization step
			$hybstep ++;
			#if doing reverse sequencing, $tilenr --

			print "Reading directory: $folder\n";
			#read all tif images in the folder
			my @listb = <"$folder/*.tif">;
			@listb = sort {$a cmp $b} @listb;
			#scalar $count to count the number of .tif files read in a folder
			$count=0;

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
					my @id = split '_', $imagename[2];
					my @subid = split 'c', $id[1];				
					my $scene = $subid[0];
					my $channel = $subid[1];
					$file[$channel-1] = $imagename[2];
								

					if ($count%5 == 0) {
					
						if ($hybstep == 1) {
							push (@refblob, $file[4]);
							$refpathname = join'', ($folder, "\/");
						}
					
						
						#contents of HybStep column must be in the format of hyb01, hyb02, etc in CP file
						my @hybstep = ('hyb0',$hybstep);
						my $hyb = join '', @hybstep;
						#the pathname requires a slash in the end
						my @pathname = ($folder,"\/");
						my $pathname = join'', @pathname;
						
						print $out "$tilenr,0,0";
						#change here if the channel setting is not the same
						@god = ($hyb, $pathname, $file[0], $pathname, $file[1], $pathname, $file[2], $pathname, $file[3], $refpathname, $refblob[$tilenr-1], $pathname, $file[4]); 
					
						foreach my $god (@god) {
							print $out ",$god";
							$printline ++;
						}
										
						print $out "\n";

					}
				}
			}
			else {
				print "Empty directory $folder!\n";
			}
		}
	}
	$tilenr ++;
}

close ($out);
