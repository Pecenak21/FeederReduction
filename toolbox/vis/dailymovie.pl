#!/usr/bin/perl

use strict;
use Cwd 'abs_path';
use Data::Dumper;
use Time::Piece;
use Time::Seconds;

my $conf_file = abs_path("$0/../")."/../../conf/local/setup.conf";
our $setup = readConf($conf_file);

# usage:
sub print_usage {
	print STDERR "usage:\n\tdailymovie.pl usi_name [yyyy-mm-dd]\n";
	exit 1;
}

# get input arguments
my $iname = shift(@ARGV);
if($iname eq '') { print_usage(); }
my $imager = siImager($iname);
if(!defined $imager) { die "couldn't load imager \"$iname\""; }

my $day = shift(@ARGV);
if(!defined $day) {
	# default to today in the local timezone
	$day = localtime(time) - ONE_DAY;
	$day = $day->strftime('%Y-%m-%d');
}
# put day in the right format
$day = Time::Piece->strptime($day, "%Y-%m-%d");

# name the output file
my $output = siNormalizePath($setup->{MOVIE_DIR});
$iname =~ /[Uu][Ss][Ii]_?(\d+)_(\d+)/;
$output .= "/usi$1-$2";
if(! -d $output) {
	mkdir($output);
}
$output .= "/" . $day->strftime('%Y%m%d') . ".mp4";
if(-e $output) {
	die "The automatically chosen output file name\n\n$output\n\nalready exists.  Cancelling";
}
print "output to $output\n";

# adjust for UTC time on the day (15 degrees of longitude per hour; just need to get close enough that this is in the middle of the night)
$day = $day - ONE_DAY * ($imager->{longitude}/15)/24;

# run the perl script
print "Preparing frames...\n";
my $framedir = prep_frames( $imager, $day, $output);

# encode the video
if( !hasx264() ) {
	die "dailymovie requires x264 to be installed.\ntry 'sudo apt-get install x264'.\n";
}
system("x264 \"$framedir/%04d.jpg\" --crf 18 -o \"$output\"");

# and clean up the temp dir
system("rm -r '$framedir'");

print "\nDone\n\n";

sub hasx264 {
# checks for presence of x264
	my $xpath = `which x264`;
	chomp $xpath;
	return ($xpath ne '');
}

sub prep_frames {
	my $imager = shift;
	my $t = shift;
	my $outputfile = shift;
	# look from one hour after midnight to one hour before midnight; this avoids getting the wrong NAS directory around times when we switch NASes
	my $tend = $t->new() + ONE_DAY * (23/24);
	$t += ONE_DAY / 24;
	# get input directories for starting and ending day
	print "looking for image dirs for " . $t->strftime('%Y-%m-%d %H:%M:%S') . ' and ' . $tend->strftime('%Y-%m-%d %H:%M:%S') . "\n";
	my $inputdir = $imager->imageDir($t) . '/' . $t->strftime('%Y%m%d');
	my $idir2 = $imager->imageDir($tend) . '/' . $tend->strftime('%Y%m%d');
	print "\t$inputdir\n\t$idir2\n";
	# construct a list of the files to work with
	# we want files in one of these two directories whose names end in '_prev.jpg' and with size >1KB (to avoid some corrupt images)
	# the stuff with sort and cut is to get them in the right order, since `find` isn't very good about that, compared to, say, ls (which in turn isn't good for getting only nonzero file sizes)
	my $filelist = `find "$inputdir" "$idir2" -iname "*_prev.jpg" -size +1k -printf '%P\t%p\n' | sort | cut -d '\t' -f 2`;
	my @files = split(/\n/,$filelist);
	printf "found %u images\n", $#files + 1;

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
	printf "linked %u frames\n", $index-1;
	return $outputdir;
}

sub readConf {
# read the named parameter out of the conf file
# if no parameter is named, returns a hashref with the entire conf file
my $cf_raw = shift();
my $cf = abs_path($cf_raw);
my $propname = shift();

open(CFH, $cf) or die "couldn't read conf file \"$cf_raw\"";
my $data = {};

my $line;
my $readMulti = 0;
my $multiKey = '';
while($line = <CFH>) {
	chomp $line;
	# remove comments
	$line =~ s/#.*//;
	# remove empty space
	$line =~ s/^\s+|\s+$//g;
	if($line eq '') { next; }
	# try reading multivalue key/value set
	if(length($line)>1 && $line =~ /^\$\$/) {
		if($line =~ /^\$\$\s*(\w+)\s*$/ && $1 ne '') {
			# put ourselves in multi-value mode, and save the key
			$multiKey = $1;
			$readMulti = 1;
			# check whether this key's already in use
			if(exists $data->{$multiKey}) { print STDERR "Multiple copies of key $multiKey; will be using the last value specified"; }
			# create an empty arrayref to save values in
			$data->{$multiKey} = [];
		} else {
			$readMulti = 0;
		}
	} elsif($readMulti) { # midway through a multiline read
		if($line !~ /^.?%/) { # skip %% lines, which were important to bryan's original java code but are legacy
			# matlab code strips leading/trailing space here, but I think that'salready been done
			# add non-empty values to the list
			if($line ne '') {
				push @{$data->{$multiKey}}, $line;
			}
		}
	} elsif(length($line)>1 && $line =~ /^%%/) {
		print STDERR "Found %% outside multiline value in config file.  File may be corrupt, or you need to remember not to use matlab style comments\n";
	} else { # regular line
		# split into a key and a value
		(my $key, my $value) = split(/\s+/, $line);
		if(exists $data->{$key}) { print STDERR "Multiple copies of key $key; will be using the last value specified"; }
		$data->{$key} = $value;
	}
}

close(CFH);
return $data;

}

sub siNormalizePath {
# resolve $ refs and convert to absolute paths
my $p = shift;

if( $p =~ /^\$/ ) {
	if( $p =~ /^\$\$\// ) {
		# $$ is a shortcut to the forecast code root
		die "looking up \$\$ paths is currently unsupported";
	} else {
		$p =~ s/^\$([^\/\\]+)//;
		my $drive = getDriveByID($1);
		$p = $drive.$p;
	}
} else {
	$p = abs_path($p);
}

return $p;
}

sub getDriveByID {
# lookup drive by drive ID
my $driveid = shift;

# get a list of mounted drives
open(MNTLST, "-|", "mount") or die "couldn't call mount";
my $drive;
my $cid;
while($drive = <MNTLST>) {
	(undef, $drive) = split(/ on /, $drive);
	($drive, undef) = split(/ (?:type|\()/, $drive);
	# look on each mountpoint for a file called drive_id
	if(-f "$drive/drive_id") {
		open(DID,"$drive/drive_id");
		$cid = <DID>; chomp $cid;
		close(DID);
		# check if the requested id matches the one from the file
		if(lc($cid) eq lc($driveid)) { last; }
	}
	# on no match, leave drive path empty
	$drive = '';
}
close(MNTLST);

return $drive;

}

sub siImager {
# load the imager data structure for the named imager
my $name = shift;
my $loc;
my $imager = undef;
for $loc (@{$setup->{siImagerSearchPath}}) {
	$loc = siNormalizePath($loc);
	if(-f "$loc/$name/imager.conf") {
		$imager = readConf("$loc/$name/imager.conf");
		last;
	}
}
if(defined $imager) {
	$imager->{name} = $imager->{NAME};
	$imager->{latitude} = $imager->{LATITUDE}+0;
	$imager->{longitude} = $imager->{LONGITUDE}+0;
	$imager->{altitude} = $imager->{ALTITUDE}+0;
	if(ref $imager->{IMAGE_DIR} ne 'ARRAY') { #assume it's a single string, need to add a date for it:
		my $old_value_tmp = $imager->{IMAGE_DIR};
		$imager->{IMAGE_DIR} = [$imager->{IMAGE_DIR}, '0'];
	}
	for (my $i = 1; $i <= $#{$imager->{IMAGE_DIR}}; $i += 2) {
		if($imager->{IMAGE_DIR}->[$i] == 0) {
			$imager->{IMAGE_DIR}->[$i] = gmtime(0); # long long ago
		} else {
			$imager->{IMAGE_DIR}->[$i] = Time::Piece->strptime($imager->{IMAGE_DIR}->[$i], "%Y-%m-%d %H:%M:%S");
		}
	}
	bless($imager,'siImager');
}
return $imager;
}

package siImager;

sub imageDir {
	my $imager = shift;
	my $day = shift;
	my @idirs = @{$imager->{IMAGE_DIR}};
	my $i;
	for ($i = 1; $i <= $#idirs; $i += 2) {
		# for now, we rely on the days being in order
		# that means the first time one of the lookup dates is greater than the day, we return the previous path
		if($idirs[$i] > $day) {
			$i-=3; # -2 for going to the time step that's actually right, then -1 for going to the path, which comes before the timestamp
			last;
		}
	}
	# if we didn't find a date that exceeds the desired one explicitly, that means the latest path is still valid
	if($i > $#idirs) { $i -= 3; }
	return main::siNormalizePath($idirs[$i]);
}
