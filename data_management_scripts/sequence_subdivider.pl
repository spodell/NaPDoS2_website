#!/usr/bin/perl
# sequence_subdivider.pl
# Sheila Podell 
# March 28, 2020

# takes as input a sequence file in fasta format
# divides up into pieces of given length


# Usage: $0 <options> fasta_filename
# 	Options:
# 	  -p	prefix	(default = "CNTRL")
# 	  -n	read length	(default = 500 nt)
# 	  -g	approximate gap length (default = 0 nt)
# 	  -d	coverage depth (default = 1X)
# 	  -o	overlap length (default = 50 nt)


use warnings; 
use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;

# Default values
	my($USAGE, $input_fasta, $prefix, $readsize, 
		$gap_target, $depth, $overlap, $mated, $qual, $xml);
	
	$prefix = "CNTRL";
	$readsize = 500;
	$gap_target = 0;
	$depth = 1;
	$overlap = 50;
	$mated = 0;	
	
	my $message = qq(\tUsage: $0 <options> fasta_filenames
	Options:
	  -p	prefix	(default = "CNTRL")
	  -n	subsequence length	(default = 500 residues)
	  -d	coverage depth (default = 1X)
	  -o	overlap length (default = 50 nt)     
)."\n";

# get input parameters from command line
	GetOptions( 
				"p=s" => \$prefix,
				"n=i" => \$readsize,
				"g=i" => \$gap_target,
				"d=i" => \$depth,
				"o=i" => \$overlap,
				'q' => \$qual,
				'x' => \$xml,
        );
        		
	$input_fasta = $ARGV[0];
	&validate_input;
	my $skip_interval = int($gap_target/$readsize) +1;
	my $read_span = int($gap_target/$readsize);
	my $actual_gap = $read_span * $readsize;
	my $insert_size = 2 * $readsize + $actual_gap;
	print STDERR "  Creating overlapping subsequences from file $input_fasta: 
  Read length = $readsize, Overlap = $overlap residues.\n";
	
# convert source sequences to objects
	my @input_seqs = &get_seq_objects($input_fasta);	
	unless (scalar @input_seqs >0)
	{
		print STDERR "ERROR: Couldn't find any fasta sequences\n";
		exit(0);
	}

# prepare output files
	my $fasta_outfile = "$prefix".".seq";
	open (FASTA_OUT, ">$fasta_outfile") or die "can't open FASTA_OUT file $fasta_outfile\n$!\n";
	
# create read objects for each sequence using sliding window			
	my $input_seq_num = 0;
	my $mate_count = 0;
	my $read_count = 0;
	for (my $i=0; $i< $depth; $i++)
	{
		foreach my $input_seq (@input_seqs)
		{
			$input_seq_num++;
			my $source_seq = $input_seq->{sequence};
			my $seq_id = $input_seq->{id};
			my @tiled_reads = &create_reads($source_seq, $overlap, $seq_id, $input_seq_num, $i);
			$read_count += scalar @tiled_reads;
		
		# write read info to output files
			foreach my $read (@tiled_reads)
			{			 
			# collect read info	
				my $seqname = $input_seq->{id};
				my $readstart = $read->{start};
				my $readend = $read->{end};
				my $read_dir = $read->{dir};
				my $readnum = $read->{readnum};
				my $read_id = $read->{id};
			#don't want mates (reverse complement) for protein subsequences
				my $mate_id = undef;
				my $read_seq = &pretty_seq($read->{sequence});
				
			# get rid of suffix for template_id	
				my @tmp = split '\.', $read_id;
				my $template_id = $tmp[0];
				unless (defined $template_id)
				{
					print STDERR "no template ID for $read_id\n";
					exit(0);
				}
				
				unless (defined $mate_id)
				{
					$mate_id = "unmated";
				}
				
				my $header =  ">$read_id  Parent seq_id $seqname, pos. $readstart-$readend, $read_dir";			
			#debug	
				#print STDERR "$header\n";
			
			# print .seq file entry
				print FASTA_OUT "$header\n$read_seq\n";
		
			# print qual file entry
				my $qual_value = 40;
			# this means likelihood of random chance is 1e-4, so all qual values are extremely good
				if (defined $qual)
				{			
					print QUAL_OUT "$header\n";
					my @seqlines = split '\n', $read_seq;
					foreach my $line (@seqlines)
					{
						$line =~ s/A/$qual_value /g;
						$line =~ s/T/$qual_value /g;
						$line =~ s/C/$qual_value /g;
						$line =~ s/G/$qual_value /g;
						
					# qual value for ambiguous nucleotides = 10	
						$line =~ s/N/10 /g;
					# get rid of terminal space	
						chop $line;
						print QUAL_OUT "$line\n";			
					}
				}
			# print XML file entry	
			
		# note: can't specify a right clip, because read length will vary 
		# unmated reads at 5 prime end will be shorter than readlength, especially if
		# depth > 1 is used
				if (defined $xml)
				{
					print XML_OUT qq(
		<trace>
			<CENTER_NAME>Artificially generated read</CENTER_NAME>
			<CENTER_PROJECT>Source filename $input_fasta | parent sequence $seqname</CENTER_PROJECT>
			<CLIP_LEFT>0</CLIP_LEFT>
			<CLIP_RIGHT></CLIP_RIGHT>
			<INSERT_SIZE>$insert_size</INSERT_SIZE>
			<INSERT_STDEV>0</INSERT_STDEV>
			<PROGRAM_ID>readmaker.pl version 1.0</PROGRAM_ID>
			<SEQ_LIB_ID>$prefix</SEQ_LIB_ID>
			<SOURCE_TYPE>Artificially generated</SOURCE_TYPE>
			<TEMPLATE_ID>$template_id</TEMPLATE_ID>
			<TI>$read_id</TI>
			<TRACE_END>$read_dir</TRACE_END>
			<TRACE_FORMAT>SCF</TRACE_FORMAT>
			<TRACE_NAME>$read_id</TRACE_NAME>
			<TRACE_TYPE_CODE>WGS</TRACE_TYPE_CODE>						
		</trace>);
				}
			}			
		}
	}
# cleanup
	if (defined $xml)
		{
			print XML_OUT qq(</trace_volume>);
		}
		close FASTA_OUT;
		close QUAL_OUT;
		close XML_OUT;

# output success message
	print STDERR "  Created $read_count reads for $prefix from $input_fasta\n";
	
##################################
# SUBROUTINES
##################################
sub validate_input
{
	if($USAGE || !($input_fasta)) 
	{
		print STDERR $message;
		exit(0);
	} 
	
	unless (defined $input_fasta && -s "$input_fasta")
	{		
		print STDERR "\tERROR: Couldn't find fasta sequence file(s) $input_fasta\n";
		print STDERR $message;
		exit(0);
	}
	
	unless ($readsize > 0)
	{
		print STDERR "\tERROR: Readsize must be > 0\n";
		exit(0);
	}
	
	unless ($gap_target >= 0)
	{
		print STDERR "\tERROR: Gap size must be greater than or equal to zero.\n";
		exit(0);
	}
	
	foreach my $file (@ARGV)
	{
		unless (defined $file && -s "$file")
		{
			print STDERR "unable to open fasta sequence file $file\n";
			print STDERR $message;
			exit(0);
		}
	}
}

sub get_seq_objects
{
	my ($infilename) = @_;
	my $current;
	my @seq_object_list = ();
	open (INPUT, "<$infilename") or die "couldn't open $infilename, $!\n";
	
	while (my $line = <INPUT>)
	{
		next if $line =~ /^\s+$/;
		chomp $line;		
		if ($line =~ />(.+)/) 
		{								
			my @tmp = split " ", $1;
			$current = new_seq("record",\@tmp);
			push @seq_object_list, $current;
		}													
		else
		{			
			$current->{sequence} .= uc"$line";
		}	
	}	
	close INPUT;
# check to be sure sequence is nucleic acid	
	my $seqstring = $current->{sequence};
	my $seq_length = length ($seqstring);
	# my $num_ATCGN = (uc$seqstring =~ tr/ATCGN//);
# 	unless ($num_ATCGN eq length $seqstring)
# 	{
# 		warn "$current->{id} doesn't look like nucleic acid\n";
# 		return undef;
# 	}
	return @seq_object_list;	
}
sub new_seq {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  $self->{id} = shift @properties;
  
  $self->{annotation} = join " ", @properties;
  unless (defined $self->{id})
  {
  	warn "no id found for sequence $self->{seq_number}\n";
  	$self->{id} = "sequence $self->{seq_number}";
  }
  
  $self->{sequence} ="";
  return($self)
}
  
sub new_read {
  my ($className, $start, $end, $sequence, $readnum, $parent_seq_id, $parent_seq_num) = @_;
  my $self = {};
  bless $self, $className;
  $self->{start} = $start;
  $self->{end} = $end;
  $self->{dir} = "forward";
  $self->{sequence} = $sequence;
  $self->{readnum} = $readnum;
  $self->{parent_seq_id} = $parent_seq_id;
  $self->{parent_seq_num} = $parent_seq_num;
  $self->{name} = "$prefix"."_S$parent_seq_num"."R$readnum";
  $self->{id} = "$self->{name}".".b1";
  
  return($self)
}
  
sub create_reads
{
	my ($seq, $olap, $parent_seq_id, $parent_seq_num, $depth_iteration) = @_;
	my $len = length $seq;
	my @forward_reads = (); # list of read objects
			
# sliding window	
	my ($winstart, $winend, $winsize, $window, $current);
	my $readnum = 0;
	my $offset = 0;
	for (my $i = $depth_iteration +1; $i <= $len-($readsize-1); $i = $i + $readsize - $olap) 
	{  	
    	$readnum++;
    	$winstart = $i;

    	if ($depth_iteration > 0)
    	{
    		$offset = int($readsize/($depth_iteration +1));
    		$winstart = $i+ $offset;    	   	
    	}
    	$winend = $winstart + ($readsize-1);
    	$winsize = $winend - $winstart;
    	$window = substr($seq, $winstart, $readsize);   	
    	
    	next if length $window < $readsize;
    	$current = &new_read("read", $winstart, $winend, $window, $readnum,
    							$parent_seq_id, $parent_seq_num);    	
    	push @forward_reads, $current;
    	
    	# debug
    		#print STDERR "winstart=$winstart winend=$winend\n";
    }    
    # Make one final read to include partial window, if necessary   
    if ($winend < $len)    
    {
    	$readnum++;
    	$winstart = $winend -$olap;
    	$winend = $len;
    	my $winlen = abs($winstart - $winend);
    	$window = substr($seq, $winstart, $winlen);
    	$current = &new_read("read", $winstart, $winend, $window, $readnum, $parent_seq_id, $parent_seq_num);
    	push @forward_reads, $current;    	
	# debug
    	#print STDERR "winstart=$winstart winend=$winend\n";
    }        
	return @forward_reads;
}

sub get_mates
{
	my ($arrayref, $skip) = @_;
	my @read_objects = @$arrayref;
	my %matepairs = ();	# key = forward_read, value = reverse_read (or undef)
	
	for (my $i=0; $i < scalar @read_objects; $i++)
	{		
		my $forward = ($read_objects[$i]);						
		my $reverse = $read_objects[$i+$skip];
	# don't overwrite anything already defined as mate by look ahead	
		next if (defined $forward->{mate});
		
	# keep forward as singleton if no reverse exists	
		unless (defined $reverse)
		{
			$matepairs{$forward} = undef;
			next;
		}
			
	# get reverse_comp of sequence	
			my $seq = $reverse->{sequence};
			my $reverse_seq = uc reverse ($seq);
			$reverse_seq =~ tr/ACGT/TGCA/;
			$reverse->{sequence} = $reverse_seq;
			$reverse->{dir} = "reverse";
			$reverse->{mate}=$forward;
			$forward->{mate}=$reverse;
			$matepairs{$reverse} =$forward;
			$matepairs{$forward} = $reverse;
	# assign id based on matepair
		$reverse->{id} = "$forward->{name}".".g1";
	
	}		
	return %matepairs;	
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

