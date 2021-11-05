#!/usr/bin/perl
#get_subsequence.pl
# Sheila Podell 
# May 30, 2005
# revised August 11, 2009

# This version correctly processes scaffolds that include "N" characters

# takes in a fasta format file
# uses tab delimited input file to get subsequences at specified coordinates:
# sourceID start stop select_id[optional]

# outputs fasta format file with subsequences

use warnings; 
use strict;

# get filenames
	unless (scalar @ARGV ==2)
	{
		print STDERR "Usage: $0 id_list fasta_filename\n";
		exit(0);
	}
	
# get selection parameters - note that each id can point to more than one record
# because hash value is an array reference

	my $select_file = $ARGV[0];
	my %select_ids = (); 
	my $num_requests = 0;
	
	open (SELECT, "<$select_file") or die "couldn't open $select_file, $!\n";
	while (<SELECT>)
	{
		chomp;
		$num_requests++;
		my @temp = split /\t/, $_;
		next unless (defined $temp[1] && $temp[1] =~ /\d+/ && $temp[2] =~ /\d+/);
		my $record = &new("record", \@temp);
		my $id = shift @temp;
		unless (exists $select_ids{$id})
		{
			$select_ids{$id} = [$record];
			next;
		}
		my $array_reference = $select_ids{$id};
		push @$array_reference, $record;
	}	
	close SELECT;
	

# go thru sequence db file, select all records that match select list 
	my $infilename = $ARGV[1];
	open (INFILE, "<$infilename") or die "couldn't open $infilename, $!\n";
	my $select = 0;		# 0=false, 1= true
	my $count_hits = 0;
	my %hits = ();		# db filter - only required seqs selected
	my $current = "";
		
	while (my $line = <INFILE>)
	{
		if ($line =~ /^\>(.+)/)
		{
			last if ($count_hits == $num_requests);
			$select = &check_name($1);
			if ($select > 0)
			{
				$current = $line;
				chomp $current;
				$current =~ s/>//;
				$hits{$current} = "";
				$count_hits++;
			}
			next;
		}				
		if ($select > 0)
		{
			chomp $line;
			$hits{$current} .= $line;
		}		
	}	
	close INFILE;

# output	
	foreach (sort keys %hits)
	{
		my @header = split / /, $_;
		my $arrayref = $select_ids{$header[0]};
		
	#debug
		unless (exists $select_ids{$header[0]})
		{
			print "nothing in %select_ids for ${_}\n";
			next;
		}
						
		foreach my $seq (@$arrayref)
		{			
			my $begin = $seq->{start};
			my $end = $seq->{stop};
  			my $length = $seq->{len};
  			my $strand = $seq->{direction};
  			print ">$seq->{name}";
  			print " [positions $begin-$end";
  			
  		# substring starts at zero, coordinates start at 1
			my $substring = uc(substr($hits{$_},$begin -1, $length));		

		# print reverse complement if string is on opposite strand
			if ($strand <0)
			{
				print ", reverse complement";
				$substring =~ tr/ACGTN/TGCAN/;
				$substring = reverse $substring;
			}
			$substring = &pretty_seq($substring);
			print "]\n$substring\n";
		}		
	}
	
#################################
# SUBROUTINES
#################################

sub new {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  # error check on input line 
	if (scalar @properties < 3)
	{
		return 0;
	} 
  my ($start, $stop, $direction);
  $self->{query} = $properties[0]; 
  #figure out direction
	if ($properties[2] > $properties[1])
		{
			$start = $properties[1];
			$stop =$properties[2];
			$direction = 1;
		}
		else
		{
			$start = $properties[2];
			$stop =$properties[1];
			$direction = -1;
		}

  $self->{start} = $start; 
  $self->{stop} = $stop;
  $self->{len} = $stop - $start +1;
  $self->{direction} = $direction;
  if (defined  $properties[3])
  {
  	 $self->{name} = "$properties[3] $self->{query}";
  }
  else
  {
  	$self->{name} = $self->{query};
  }
  return($self)
}

sub check_name
{
	my ($header) = @_;
	my $result = 0;
	chomp $header;
	my @seq_words = split (" ", $header);
	my $seq_id = $seq_words[0];
		
	if (exists $select_ids{$seq_id})
	{
		$result = 1;
	}
	return $result;		
}


sub pretty_seq
{
	my ($seq) = @_;
	my $pretty_seq = "";
	my $count = 0;
	my $line_length = 60;
	$seq =~ s/ //g;
	while ($count < length $seq)
	{
		$pretty_seq .= substr ($seq, $count, $line_length);
		$pretty_seq .= "\n";
		$count += $line_length;
	}
	chomp $pretty_seq;
	return $pretty_seq;
}

__END__
