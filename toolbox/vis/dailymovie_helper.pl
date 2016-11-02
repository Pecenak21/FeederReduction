#!/usr/bin/perl
use strict;
use Time::Piece;
use Time::Seconds;

# Here's our brief usage statement
my $usage = "Usage: $0 usi_dir start_time output_path\n";
# we want the usi image directory (we'll pick the date folder), the starting time, and a path to save the movie to

# capture our inputs
my $inputdir = shift(@ARGV) or die "Must supply an input directory\n$usage";
chomp($inputdir = `readlink -f $inputdir`);
my $starttime = shift @ARGV or die "Must supply a start time\n$usage";
my $outputfile = shift @ARGV or die "Must supply an output filename\n$usage";

# construct a list of the input files
my $t = Time::Piece->strptime($starttime,"%Y-%m-%d %H:%M:%S");
my $idir2 = $inputdir;
$inputdir .= '/' . $t->strftime('%Y%m%d');
my $tend = $t->new();
$tend += ONE_DAY;
$idir2    .= '/' . $tend->strftime('%Y%m%d');
my $filelist = `ls "$inputdir/"*_prev.jpg "$idir2/"*_prev.jpg`;
my @files = split(/\n/,$filelist);

# create a temp directory to work in - we're going to make a bunch of symlinks that we don't want to clutter up other parts of the disk with
my $outputdir = "/tmp/makemovie/" . `uuidgen`;
chomp $outputdir;
system("mkdir -p $outputdir");

# symlink the movie frames into a temp directory
my $file;
my $index = 1;
foreach $file (@files) {
	$file =~ /(\d{14})_prev.jpg$/;
	my $ftime = Time::Piece->strptime($1,"%Y%m%d%H%M%S");
	if($ftime < $t || $ftime > $tend) { next; }
	my $dstname = sprintf("%s/%04d.jpg",$outputdir,$index);
	$index ++;
	system("ln -s \"$file\" \"$dstname\"");
}

# try passing the output to matlab so we can do the video encoding directly from there
print $outputdir;
exit 0;

# output the video
system("x264 $outputdir/%04d.jpg --crf 18 -o \"$outputfile\"");

# and clean up the temp dir
rmdir $outputdir;
