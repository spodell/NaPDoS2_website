#!/usr/bin/perl 
# select_tab_instances.pl 
# Sheila Podell 
# February 18, 2011

#	Takes in a file in tab format 
#	selects first X instances of each category in column 0
# (assumes file has been pre-sorted by some other criteria
# in cases of multiple occurrences for col 0)

use warnings; 
use strict;

my $num_tablines = 0;
my $num_hits = 0;
my $max_num_instances = 1;

# Get filenames & open filehandles 
	unless (scalar @ARGV == 2)
	{
		print "Usage: $0 tab_filename max_num_instances\n";
		exit (0);
	}
	
	my $in_filename = $ARGV[0];
	$max_num_instances = $ARGV[1];
	open (INPUT,"$in_filename") or die "can't open input file\n$!\n";	
			
# Check for keys
	my %keywords = ();
	while (<INPUT>) 
	{
		next if ($_ =~ /^\s/);
		chomp;
		#$num_tablines++;
		my @tmp = split "\t", $_;
		$keywords{$tmp[0]}++;
		if ($keywords{$tmp[0]} <= $max_num_instances)
		{
			print  "$_\n";
		}
	}	
	
#user output		
	close INPUT;
	print "\n";
			
__END__
