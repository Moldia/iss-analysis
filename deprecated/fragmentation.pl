#!/usr/bin/perl -w 

use FileHandle;

# Subroutine prototypes
sub usage;

my ($next_arg, $infile, $outfile, $size, $out_fh );

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-s")    { $size = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outfile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ( !$infile || !$size ) { &usage(); }

# Create input & output filehandle
my $in_fh = new FileHandle("$infile") or die "\n\nCouldn't open input file ".$infile."!\n\n";
if($outfile) {
	$out_fh = new FileHandle(">$outfile") or die "\n\nCouldn't open output file ".$outfile."!\n\n";;
}
else {
	$out_fh = "/dev/stdout";
}


my ($id, $seq);
while (<$in_fh>) {
	if($_=~m/^#/ || $_=~m/^$/) { next; }
	
	chomp;
	if ($_ =~ m/^>/) {
		$id = $_;
	}	
	else {
		$seq = $_;
		my $seqLength = length($seq);
		my $fragStart = 0;
		my $fragEnd = $size;
		while ( $fragEnd <= $seqLength ) {
			print $out_fh $id."_frag".($fragStart+1)."\n";
			print $out_fh substr($seq,$fragStart,$size/2)."|".substr($seq,($fragStart+$size/2),$size/2)."\n";
#			print $out_fh substr($seq,$fragStart,$size)."\n";
			$fragStart++;
			$fragEnd++;
		}
	}
}
$in_fh->close;
$out_fh->close;

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <infile> -o <outfile> -s <size of sliding window>\n 
 
 -i Infile
 -o Outfile [optional, default stdout]
 -s Size of the sliding window \n\n";
  exit(1);
}
