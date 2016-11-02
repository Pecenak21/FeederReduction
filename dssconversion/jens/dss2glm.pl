#!/usr/bin/perl

# Script to read an openDSS file (give it a single master file, it will read any redirects in the master file) and generate a (hopefully) equivalent GridLAB-D model.  At this time, the GridLAB-D model is output to STDOUT, so you should probably redirect it to a file.

use strict;
use File::Spec;

my $file = shift;
our $classes = {};

# Start by reading the openDSS model from a file
# the reader returns a hash reference to a model
my $model = readDSS($file);
# opendss waits to construct a list of busses until after the model is specified and a certain command is run.  We need to have a list before we can export
# OpenDSS does this in response to specific commands, but we'll do it here to make sure that it happens before we try to export
constructBusses($model);

# next, we convert to a GridLabD style model.  This is where all the hard work happens
my $gldmodel = convertDSSToGLD($model);

# Finally, write the model data out in GLM format. Should be pretty easy
writeGLM($gldmodel);

exit 0;

sub readDSS {#(filename, [model_to_add_to]) 
	my $filename = shift;
	if(!(my $model = shift) ) {
		$model = {};
	}
	(my $v, my $d, $filename) = File::Spec->splitpath($filename);
	$d = File::Spec->catpath($v,$d,'');
	chdir($d) if ($d ne '');
	open(my $ifh, $filename) or die "failed to open file \"$filename\"";
	
	my $nl;
	my $lastobj = '';
	while($nl = <$ifh>) {
		# separate the command from the rest of the line; skip comments and empty lines
		(my $cmd, $nl) = split(/\s+/,$nl,2);
		if ($cmd eq "" || $cmd =~ /^\/\// || $cmd =~ /^\!/) { next; }
		# remove any trailing comments
		($nl, undef) = split(/\s*(\!|\/\/)/,$nl);
		$cmd = lc($cmd);
		if(defined $nl) { #this block silences a warning which occurs when a command has no arguments (e.g. solve, clear)
			$nl =~ s/[\r\n]+$//;
		} else {
			$nl = '';
		}
		if($cmd eq "redirect") {
			print STDERR "reading file ${nl}\n";
			my $cwd = `pwd`;
			readDSS($nl,$model);
			chdir $cwd;
			$lastobj = '';
		} elsif($cmd eq "clear") {
			print STDERR "Clearing Model!\n";
			$model = {};
		} elsif($cmd eq "new") {
			$lastobj = newObj($nl);
			my $class = $lastobj->{class};
			unless( defined $model->{$class} ){
				$model->{$class} = {};
			}
			if( defined $model->{$class}->{$lastobj->{name}} ) {
				print STDERR "two objects named $class:$lastobj->{name}\n";
			} else {
				$model->{$class}->{$lastobj->{name}} = $lastobj;
			}
		} elsif($cmd eq "set") {
			# only one of these, so it's low priority
			$lastobj = '';
			print STDERR "set $nl\n";
		} elsif($cmd eq "buscoords") {
			# only one of these as well
			$lastobj = '';
			# For now we just make sure we have the full path to the bus coords file, and we'll actually read it in later when we're constructing the bus list
			$cmd = `pwd`; chomp $cmd;
			$model->{buscoords} = "$cmd/$nl";
		} elsif($cmd eq "calcvoltagebases") {
			$lastobj = '';
			print STDERR "figure out if we need to calc base voltages\n";
		} elsif($cmd eq "~" || $cmd =~ /~/) {
			next if $lastobj eq '';
			$nl =~ s/\s*=\s*/=/g;
			$nl =~ s/\s+([\)\]\}])/$1/g;
			my @line = split /[\s,]+/, $nl;
			addToObj($lastobj,@line);
		} else {
			$lastobj = '';
			print STDERR "unknown cmd: $cmd\n";
		}
	}
	
	close($ifh);
	return $model;
}

sub newObj {#(nl)
	my $nl = shift;
	my $lastobj = {};
	$nl =~ s/\s*=\s*/=/g;
	$nl =~ s/\s+([\)\]\}])/$1/g;
	my @line = split /[\s,]+/, $nl;
	my $classname;
	if($nl =~ /object=([\w\.]+)/i) {
		$classname = $1;
	} else {
		$classname = shift(@line);
	}
	($lastobj->{class}, $lastobj->{name}) = split(/\./, $classname,2);
	$lastobj->{class} = lc($lastobj->{class});
	addToObj($lastobj,@line);
	return $lastobj;
}

sub addToObj{#($lastobj,@line)
	my $lastobj = shift;
	my $prop;
	my $pname;
	my $lastname = 'undef';
	my @fieldorder = ("r1", "x1", "r0", "x0", "c1", "c0", "normamps", "emergamps", undef);
	my $twdg = undef;
	while ($prop = shift) {
		($pname, $prop) = split /=/,$prop;
		if(!defined $prop) {
			$prop = $pname;
			$pname = undef;
			foreach my $i (0..$#fieldorder) {
				if($lastname eq $fieldorder[$i]) {
					$pname = $fieldorder[$i+1];
					last;
				}
			}
			if(!defined($pname)) {
				print STDERR "property that comes after \"$lastname\" skipped. val:\"$prop\"\n";
				next;
			}
		}
		if(($prop =~ /[\(\[\{]/) && ($prop !~ /[\)\]\}]/)) {
			my $tstr;
			my $done = 0;
			while($tstr = shift(@_)) {
				$prop .= " ".$tstr;
				if($tstr =~ /[\)\]\}]/) { $done = 1; last; }
			}
			print STDERR "missing closing ) of list.  so far: \"$prop\"\n" unless $done;
		}
		$pname = lc($pname);
		$lastname = $pname;
		if(lc($pname) eq 'wdg' && $lastobj->{class} eq 'transformer') {
			$twdg = $prop-1; # convert to an index
			next;
		}
		
		if (defined $twdg) {
			unless(defined($lastobj->{$pname})) {
				$lastobj->{$pname} = [];
			}
			$lastobj->{$pname}->[$twdg] = $prop;
		} elsif (!defined($lastobj->{$pname})) {
			$lastobj->{$pname} = $prop;
		} else {
			print STDERR "duplicate property \"$pname\"\n";
		}
	}
}

sub constructBusses {#(model)
	my $model = shift;
	my $busses = {};
	print STDERR "creating busses...\n";
	foreach my $class (keys %$model) {
		next if $class eq 'buscoords';
		foreach my $obj (values %{$model->{$class}}) {
			my $bus_n;
			if( defined $obj->{busses} ) {
				$bus_n = $obj->{busses};
			} elsif (defined $obj->{bus}) {
				$bus_n = $obj->{bus}
			} elsif (defined $obj->{bus1}) {# lines, capacitors, and maybe some other objects do it this way
				$bus_n = [$obj->{bus1}];
				push(@$bus_n,$obj->{bus2}) if(defined, $obj->{bus2});
			} elsif (defined $obj->{bus2}) {#hopefully this never happens, but if some objects only have a bus2...
				$bus_n = [$obj->{bus2}];
			} else { next; }
			if(ref($bus_n) ne 'ARRAY') {
				$bus_n = [split(/ +/,$bus_n)];
			}
			foreach my $bus_name (@{$bus_n}) {
				# remove any leading or trailing grouping operators
				$bus_name =~ s/[\(\{\[\]\}\)]//g;
				# get phases and name separately
				my @phases = split(/\./,$bus_name);
				$bus_name = shift @phases;
				# save the phases for later
				unless(defined $busses->{$bus_name}) { $busses->{$bus_name} = {phases => []}; }
				push(@{$busses->{$bus_name}->{phases}}, @phases);
				unless(defined $busses->{$bus_name}->{usedby}) { $busses->{$bus_name}->{usedby} = []; };
				push(@{$busses->{$bus_name}->{usedby}},$obj);
			}
		}
	}
	$model->{busses} = $busses;
	return unless (defined $model->{buscoords});
	
	# load bus coordinates too
	my $filename = $model->{buscoords};
	delete $model->{buscoords};
	print STDERR "Load bus coordinates from $filename...\n";
	open(my $ifh, $filename) or die "failed to open file \"$filename\"";
	while(my $nl = <$ifh>) {
		$nl =~ s/[\r\n]+$//;
		# skip comments and blank lines
		if ($nl =~ /^\s?$/ || $nl =~ /^\/\// || $nl =~ /^\!/) { next; }
		(my $bus_n, my $lat, my $lon) = split(/\s?,\s?/,$nl);
		if(!(defined $bus_n && defined $lat && defined $lon)) {
			print STDERR "bad line in buscoords file:\n\t$nl\n";
		} elsif( !defined $busses->{$bus_n} ) {
			print STDERR "\tunknown bus $bus_n\n";
		} else {
			$busses->{$bus_n}->{latitude} = $lat;
			$busses->{$bus_n}->{longitude} = $lon;
		}
	}
	close($ifh);
}

# this function tries to convert an opendss model to gridlabd format
# of note, it is intentionally destructive to the input data so that it's easy to see what things it didn't know how to deal with.
sub convertDSSToGLD {#(dss)
	my $dss = shift;
	my $gld = {};
	# I was going to use a hashmap of functions by class, but I'm now going to go through them in a specific order so that I can make sure I get to some of them before others
	# start with the loads (no particular reason)
	if(defined($dss->{load})) {
		# first call the conversion function
		$gld->{load} = convertLoads($dss->{load});
		# then delete the loads so we don't get a nasty warning about them later
		delete $dss->{load};
	}
	foreach my $class (keys %$dss) {
		print STDERR "I don't know how to deal with the \"$class\" class yet.\n";
	}
	return $gld;
}

sub writeGLM {#(model)
	my $model = shift;
	print "module powerflow {\n\tsolver_method FBS; //pick a solver\n}\n\n";
	foreach my $class (keys %$model) {
		foreach my $obj (values %{$model->{$class}}) {
			print "object $obj->{class} {\n\tname $obj->{name};\n";
			delete $obj->{class};
			delete $obj->{name};
			foreach my $key (keys %$obj) {
				my $value = $obj->{$key};
				# quote the value if necessary
				if($value =~ /[\s,;]/) { $value = '"'.$value.'"'; }
				print "\t$key $value;\n";
			}
			print "}\n";
		}
	}
}

sub convertLoads {#(loads)
	my $loads = shift;
	my $gldloads = {};
	foreach my $l (values %$loads) {
		my $nl = {name => $l->{load}, class => 'load'};
		delete $l->{class}; delete $l->{name};
		$nl->{parent} = $l->{bus1}; delete $l->{bus1};
		#$nl->{
		foreach my $key (keys %$l) {
			print STDERR "I don't know how to deal with load property \"$key\" yet\n";
		}
		$gldloads->{$nl->{name}} = $nl;
	}
	return $gldloads;
}