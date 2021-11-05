#!/usr/bin/perl
# randomize_lines.pl 
# Sheila Podell 
# September 29, 2008

#	takes in a list of lines 
#	randomizes order for output (to STDOUT) 
# 	suitable for really long lists

use warnings; 
use strict;


# get number of lines that need to be randomized
	 unless (scalar @ARGV == 1)
 	{
 		print STDERR "Usage: $0 filename\n";
 		exit(0);
	}	
	my @lines =();
 	open (INPUT,"$ARGV[0]") or die "can't open input file $ARGV[0]\n$!\n"; 		 	
 	while (<INPUT>) 
 	{ 
		chomp;
		push @lines, $_;
 	}
 	close INPUT;
 	my $num_lines = scalar @lines;
		
# generate list of random index numbers
	my @list = &get_rand_nums($num_lines);	
	
#print to stdout
	foreach (@list)
	{
		print "$lines[$_]\n";
	}
	

####################
# SUBROUTINES
####################
sub get_rand_nums
{
	my ($num_to_get) = @_;
	my $i;
	
# create an array of empty slots to hold the numbers
# actual numbers go in the ordered list
# if rand number corresponds to a full slot, 
# get the next empty one below it. 
# when get to first position @array[0] in the array of slots,
# negative numbers will start filling from the last position @array[-1]
	my @slots = ();
	for ($i=0;$i< $num_to_get;$i++)
	{
		$slots[$i] = -1;
	}
	
	my @ordered_list = ();	
	my $next;
	my $safety = $num_to_get;
	my $counter = 0;
	my $num_success = 0;

	srand;
		
	while (scalar @ordered_list < $num_to_get)
	{ 
		$next = int(rand ($num_to_get -1));
		my $last_array_pos = $next - $num_to_get;
		$counter++;
		if ($slots[$next] < 0)
		{
			$slots[$next] = 1;
			push @ordered_list, $next;
			$num_success++;
		}
		else	
		{
			for ($i = $next; $i > $last_array_pos; $i--)
			{
				# if slot is unfilled, fill it, put position # into ordered list 
				if ($slots[$i] < 0)
				{					
					my $actual_position = $i;
					$slots[$actual_position] = 1;
					if ($actual_position < 0)
					{
						$actual_position = $num_to_get + $i;
					}
					
					push @ordered_list, $actual_position;
					$num_success++;
					$i = $last_array_pos;
				}		
			}			
		}			
	}
	return @ordered_list;
}
