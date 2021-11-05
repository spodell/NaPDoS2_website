#!/usr/bin/perl
# size_limit_seqs.pl
# Sheila Podell 
# January 15, 2020

# inputs:
#	a file containing multiple fasta format sequences
# 	minimum sequence size to select
# outputs selected records to a new file
#
# suitable for really large (gigabyte) fasta source files
# and large numbers of retrieval sequences

use warnings; 
use strict;

# validate input
	unless (scalar @ARGV ==2)
	{
		print STDERR "Usage: $0 fasta_filename minimum_seq_length\n";
		exit(0);
	}

	my $infilename = $ARGV[0];
	my $min_seq_size = $ARGV[1];
	my $outfilename = "$infilename"."_selected";	
	unless (defined $infilename && -s $infilename)
	{
		print STDERR "Usage: $0 fasta_filename num_seqs\n";
		print STDERR "Invalid fasta filename.\n";
		exit (0);
	}
	unless ($min_seq_size =~ /^\d+$/ && $min_seq_size > 0)
	{
		print STDERR "Usage: $0 fasta_filename num_seqs\n";
		print STDERR "Invalid sequence size.\n";
		exit (0);
	}	

	my $num_records = `grep -c ">" $infilename`;
	chomp $num_records;	
	unless ($num_records > 0)
	{
		print STDERR "No fasta records found in $infilename.\n";
		exit (0);	
	}
	print STDERR "$infilename contains $num_records fasta records\n";
	print STDERR "Minimum size requested = $min_seq_size\n";
	

# go thru infile (once only), select all records that match criterion 
	open (INFILE, "<$infilename") or die "couldn't open $infilename for reading, $!\n";
	open (OUTFILE, ">$outfilename") or die "couldn't open $outfilename for writing, $!\n";
	my $count_hits = 0;
	my $line_num = -1;
	my $record_num = -1;
	my %seq_objects =();	#key = ID number value = sequence object
	my $current = ""; #current sequence object
	
	while (my $line = <INFILE>)
	{	
	# get rid of DOS line endings	
		$line =~ s/\r/\n/sg;
		if (eof)
		{
			$record_num++;
			chomp $line;
			$line =~ s/ //g;
			$current->{sequence} .= uc "$line";	
			my $sequence = $current->{sequence};
 			my $seq_length = length $sequence;		
		    if (defined $current)
			{				 				
 				if (defined $sequence && $seq_length > $min_seq_size)
 				{
 					my $outseq = &pretty_seq("$sequence");			
 					print OUTFILE "$current->{header}";
 					print OUTFILE "\n";
 					print OUTFILE "$outseq\n";
 					$count_hits++;
				}				
 			}	
		}		
		elsif ($line =~ /^\>(.+)/)
		{
			$record_num++;
			if (defined $current && $. >1)
			{				
				my $sequence = $current->{sequence};
				my $seq_length = length $sequence;
				if (defined $sequence && $seq_length > $min_seq_size)
				{
					my $outseq = &pretty_seq("$sequence");			
					print OUTFILE "$current->{header}";
					print OUTFILE "\n";
					print OUTFILE "$outseq\n";
					$count_hits++;				
				}
				%seq_objects =();
			}						
			my @tmp = split " ", $line;
		# clear previous seq_object from memory
			%seq_objects =();
			$current = new_seq("record",\@tmp);
			next;
		}				
		else
		{
			chomp $line;
			$line =~ s/ //g;
			$current->{sequence} .= uc "$line";	
		}		
	}

	print STDERR "$count_hits/$record_num records selected.\n";	
	close INFILE;
	close OUTFILE;
		
#################################
# SUBROUTINES
#################################

sub new_seq {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  my $header = join " ", @properties;
  $self->{header} = $header;  
  $properties[0] =~ s/\>//;
  $self->{id} = shift @properties;   
  $self->{annotation} = join " ", @properties;
  $self->{sequence} ="";
  $self->{seq_length} ="";  
  return($self)
}

sub pretty_seq
{
	my ($seq) = @_;
	my $pretty_seq = "";
	my $count = 0;
	my $line_length = 70;
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


