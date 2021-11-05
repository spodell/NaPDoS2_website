#!/usr/bin/perl
# transeq_to_multifasta.pl
# Sheila Podell 
# August 9, 2008

#	takes as input an amino acid sequence file in fasta format 
#	Each sequence may contain multiple coding regions, separated
# 	by stop codons (e.g. 6 frame output from EMBOSS TRANSEQ)

#	splits each fasta sequence on stop codons (*)


#	removes any sequences with zero length (annotation only)
# 	and any sequences shorter than specified on command line 
#	(by default, removes any with 8 characters or fewer)

# 	labels each saved sequence with ID number according to sequential
# 	position on parent contig

#	outputs selected sequences to STDOUT

use warnings; 
use strict;

# get in & out file handles
	unless (@ARGV == 1 || @ARGV == 2)
	{
		print "usage: $0 <fasta_filename> <minimum_length>\n";
		exit (0);
	}	
	my $infilename = $ARGV[0];
	my $min_len =  $ARGV[1] || 8;
	open (INFILE, "<$infilename") or die "couldn't open $infilename, $!\n";	
	my $flag = 0;
	my @seq_list = ();  # list of sequence objects
	my $current = "";
	
# get sequence, id, annotation
	while (my $line = <INFILE>)
	{
		next if $line =~ /^\s+$/;
		chomp $line;		
		if ($line =~ />(.+)/) 
		{								
			my @tmp = split " ", $line;
			$current = new_seq("record",\@tmp);
			push (@seq_list, $current);
			$flag = 1;
		}													
		else
		{			
			$current->{sequence} .= uc"$line";
		}	
	}
	close INFILE;
	
# split each input sequence into subsequences, if contain internal stops
	my @subseq_list = ();	# list of sequence objects
	foreach my $next (@seq_list)
	{		
		my @subseqs = split /\*/, "$next->{sequence}";
		my $subseq_num = "0";
		foreach my $subseq (@subseqs)
		{
			next unless (defined $subseq && length $subseq > $min_len);
			my $new_annotation = "$next->{id}".".$subseq_num"." $next->{annotation}";
		# remove translation frame info, so all translated aa seqs will have same id num
			my $new_id = "$next->{id}"."_$subseq_num";
			# if ($new_id =~ /(.+)_\d/)
# 			{
# 				$new_id = $1;
# 			}
			
			my $new_header = "$new_id"." $new_annotation";
			my @tmp = split " ", $new_header;
			$current = new_seq("record",\@tmp);
			push (@subseq_list, $current);
			$current->{sequence}=$subseq;
			$subseq_num++;
		}
	}

	foreach my $next (@subseq_list)
	{
		my $seq = &pretty_seq($next->{sequence});
		print ">$next->{header}\n$seq\n";			
	}
	
###################################################
# SUBROUTINES
###################################################
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
