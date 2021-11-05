#!/usr/bin/perl 
# select_non-overlapping_tophits.pl 
# Sheila Podell 
# January 17, 2021

# Purpose: from blastp output, select non-overlapping matches that may be from same query
# simple tophits would select only one match per query sequence
# note: tricky parameter: how many mismatches and gaps are reasonable? 
	#excludes sequences withan arbitrary number of mismatches (set as global parameter)
	# default = 200
	# can be over-ridden at command line ARGV[2]
	
# Takes as input 
	# a tab-delimited text file of blast output (m8 or -outfmt 6)
	# a number for max coordinate overlap (e.g 5 residues)
# Returns as output
	# selected lines of tab-delimited text
	# tophits for each query sequence with non-overlapping coordinates

use warnings; 
use strict;

# Get filenames & open filehandles 
	unless (scalar @ARGV == 2)
	
	{
		print "Usage: $0 tab_filename max_coord_overlap\n";
		exit (0);
	}
	
	my $in_filename = $ARGV[0];
	my $max_overlap = $ARGV[1] || 10;
	my $max_mismatches = $ARGV[2] || 200;
	my $max_gap_difference = $ARGV[3] || 3;
	
# validate input
	unless (-s "$in_filename")
	{
		print "No input data in $in_filename \n";
		print "Usage: $0 tab_filename max_coord_overlap\n";
		exit (0);
	}
	unless ($max_overlap =~ /^\d+$/)
	{	
		print "Maximum overlap needs to be a positive integer (default = 10) \n";
		print "Usage: $0 tab_filename max_coord_overlap\n";
		exit (0);
	}
	
# sort the lines in order of bitscore
	my $sorted_input = `sort  -nr -k12 $in_filename`;
	my @bitscore_sorted_lines = split "\n", $sorted_input;
	my $current_query = "";
	my @ordered_query_objects = ();
	my @occupied_positions = ();
	my %query_objects = ();	# key = query_name # value = object reference (contains all lines for that query)
	my $line_num = 0;
	my $count = 0;

# create an object for each query containing associated set of bitscore-ordered lines
	foreach (@bitscore_sorted_lines) 
	{
		next if ($_ =~ /^\s/);
	# get rid of DOS line endings	
		$_ =~ s/\r/\n/sg;		
		chomp;
		
		my @params = split "\t", $_;
		my $num_params = scalar @params;
		my $query_id = $params[0];
		unless (defined $query_id && $num_params == 12)
		{
			print STDERR "Incorrect input format, line $.:\n$_\n";
			exit(0);
		}
		$line_num ++;
			
	# add line to query object, including a list of lines sorted by descending bitscore	
		if (!(exists $query_objects{$query_id}))
		{
			$current_query = &new_query("record", \@params);
			$query_objects{$query_id} = $current_query;
			push @ordered_query_objects, $query_id;	
		}
		else 
		{
			$current_query = $query_objects{$query_id}; 
		}
		# add current line to $records{$key} object	
		push @{$current_query->{all_lines}}, "$_";				
	}

# select non-overlapping lines from query objects according to start sites
	my @saved_lines = ();
	foreach my $query (@ordered_query_objects)
	{
		my @query_lines = @{$query_objects{$query}->{all_lines}};
		my $hit_num = 1;	
	
	#debug - check if lines are still in correctly sorted order (desc bitscore)
		 #print STDERR "\nANALYZING: $query query lines order: \n";
	
	# clear out occupied positions from previous query	
		@occupied_positions = ();
	
	# first line should always be saved (lines already been sorted by decreasing bitscore) 
	# use shift to remove it from previous analysis
		my $tophit_line = shift @query_lines;
		push @saved_lines, $tophit_line;
		#print STDERR "Line 98 saving tophit line $tophit_line\n";

	# store position of saved line
		my @tmp1 = split "\t", $tophit_line;
		my $tophit_start = $tmp1[6];
		my $tophit_end = $tmp1[7];
		my $tophit_mismatches = $tmp1[4];
		my $tophit_gaps = $tmp1[5];
				
	# reverse start and end, if in opposite orientation
		if ($tmp1[6] > $tmp1[7])
		{
			$tophit_start = $tmp1[7];
			$tophit_end = $tmp1[6];
		}
		#print STDERR "tophit start=$tophit_start, end=$tmp1[8], mismatches=$tophit_mismatches gaps =$tophit_gaps\n";	
			
	# adjust positions outward by overlap factor			
		my $prev_min = $tophit_start - $max_overlap;
		if ($prev_min < 0)
		{
			$prev_min = 0;
		}
		my $prev_max = $tophit_end + $max_overlap;
	
	# fill occupied_positions array with current line
		for (my $i=$prev_min; $i<=$prev_max ; $i++)
		{
			$occupied_positions[$i] = 1;
		}		 
  
	# identify additional tophits from remaining lines from query with different start sites
	foreach my $line (@query_lines)
	{
		my @tmp = split "\t", $line;
		my $current_mismatches = $tmp[4];
		my $current_gaps = $tmp[5];
		my $current_bitscore = $tmp[11];
		my $current_start = $tmp[6];
		my $current_end = $tmp[7];	
	
	# reverse start and end, if in opposite orientation
		if ($tmp[6] > $tmp[7])
		{
			$current_start = $tmp1[7];
			$current_end = $tmp[6];
		}				
		my $current_min = $current_start - $max_overlap;
		if ($current_min < 0)
		{
			$current_min = 0;
		}
		my $current_max = $current_end + $max_overlap;
		my $mismatch_difference = $tophit_mismatches - $current_mismatches;
		my $gaps_difference = $current_gaps - $tophit_gaps; 	
		
	# check to see if current start and stop positions are already occupied
	# get rid of alignments that have an unreasonable number of gaps
		if (
				(defined $occupied_positions[$current_start])
				 || (defined $occupied_positions[$current_end])
				 || ($gaps_difference > $max_gap_difference)
			)
		{
			#print "found occupied positions $current_min,$current_max \n";
			next;
		}		
		
		elsif ($mismatch_difference < $max_gap_difference) #($tophit_mismatches  + $max_overlap)) 
		{		 
				push @saved_lines, $line;
			# fill occupied_positions array with current line
				for (my $i=$current_min; $i<=$current_max ; $i++)
				{
					$occupied_positions[$i] = 1;
				}
		}	
	}	
 }

# send saved lines to stdout
	foreach my $saved (@saved_lines)
	{
		print "$saved\n";	
	}	
	
#####################
# SUBROUTINES
#####################

sub new_query {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  my $line = join "\t", @properties;
  my @all_lines = ();
  my @selected_lines = ();
  $self->{query_id} = $properties[0];
  $self->{all_lines} = \@all_lines;
  $self->{selected_lines} = \@selected_lines;
  #$self->{count_query_lines} = scalar @$self->{all_lines};

  #$self->{max_line} = $line; 
  #$self->{select_value} = $properties[$selected_column]; ;

  return($self)
}

sub new_line {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  $self->{line} = join "\t", @properties;
  $self->{query_id} = $properties[0];
  $self->{start} = $properties[6];
  $self->{end} = $properties[7];
  $self->{bitscore} = $properties[11];
  return($self)
}
		
__END__		
	
	