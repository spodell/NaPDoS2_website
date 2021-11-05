#!/usr/bin/perl 
# select_tab_quant.pl 
# Sheila Podell 
# December 23, 2004

#	Takes in a file in tab format. Selects all lines
#	where a particular field meets a quantitative
# 	threshold

#   this version writes to STDOUT

use warnings; 
use strict;

# Get filenames & open filehandles 
	unless (scalar @ARGV == 3)
	{
		print "Usage: $0 tab_filename column_number[0-256] threshold\n";
		exit (0);
	}
	
	my $in_filename = $ARGV[0];
	my $selected_column = $ARGV[1];
	
	unless ($selected_column =~ /^\d+$/ 
			&& $selected_column > -1 
			&& $selected_column < 257)
	{
		print "$ARGV[1] is invalid column range (must be 0-256)\n";
		exit (0);
	}
	
	my $threshold = $ARGV[2];
	my $out_filename = substr ($in_filename, 0, 10);
	$out_filename .= ".select";
	open (INPUT,"$in_filename") or die "can't open input file\n$!\n";	
	#open (OUTPUT, ">$out_filename") or die "can't open outfile $out_filename\n $!\n";	

#set-up data structures
	
	my $num_tablines = 0;
	my $num_hits = 0;

# Check for keys	
	while (<INPUT>) 
	{
		next if ($_ =~ /^\s/);
		chomp;
		$num_tablines++;
		my @tmp = split "\t", $_;
		if ($tmp[$selected_column] && $tmp[$selected_column] > $threshold)
		{
			$num_hits++;
			print "$_\n";			
		}				
	}	
	
#user output
	print STDERR "$num_hits hits from $num_tablines lines written to $out_filename\n";		
	close INPUT;
	#close OUTPUT;	
		
__END__
