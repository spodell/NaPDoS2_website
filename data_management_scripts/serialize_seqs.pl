#!/usr/bin/perl
# serialize_largeseqs.pl
# Sheila Podell 
# March 21, 2012

#	Breaks file containing multiple fasta format sequences
#	into sequential subfile chunks, each containing a specified 
#   number of the original sequences. 
#	Writes each chunk to output file in fasta format.

# suitable for really large (multigigabyte) fasta source files

use warnings; 
use strict;
use File::Basename;

# get cmd line params
	unless (scalar @ARGV ==2)
	{
		print STDERR "Usage: $0 fasta_filename num_seqs_per_subfile\n";
		exit(0);
	}

# get number of input fasta sequences (records), number of subfile requests
	my $infilename = $ARGV[0];
	my $num_seqs_per_subfile = $ARGV[1];
	unless ($num_seqs_per_subfile =~ /^\d+$/ && $num_seqs_per_subfile > 0)
	{
		print STDERR "Usage: $0 fasta_filename num_seqs_per_subfile\n";
		print STDERR "Invalid number of sequences requested.\n";
		exit (0);
	}
	
	my $num_records = `grep -c ">" $infilename`;
	chomp $num_records;
	print STDERR "$infilename contains $num_records fasta sequences\n";
	print STDERR "Number of sequences/subfile = $num_seqs_per_subfile\n";	
	unless ($num_records > 0 && $num_records >= $num_seqs_per_subfile)
	{
		print STDERR "Can't do this.\n";
		exit (0);	
	}	
	
# get basename for input file (strip out path if it was included in cmd line)
	my $base_name = basename($infilename);

# calculate number of records per subfile, create array of subfile names
	my @subfile_names = ();
	my ($start_num, $end_num, $filename);
	for (my $i = 0; $i < $num_records; $i+= $num_seqs_per_subfile)  
	{
		$start_num = $i+1;
		$end_num = $i + $num_seqs_per_subfile;
		
		if ($end_num > $num_records)
		{		
			$end_num = $num_records;
		}
		my $nextname = "$base_name".".$start_num"."-"."$end_num";
		push @subfile_names, $nextname; 
	}
	my $num_subfiles = scalar @subfile_names;;
	print STDERR "Total number of subfiles = $num_subfiles\n";
	
# add a reality check in case number of files is really large
	if ($num_subfiles > 100)
	{
		print STDERR "Are you sure you want to create this many new files? (y/n)\n";
		my $response = "";
		chomp ($response = <STDIN>);
		unless ($response eq "y" || $response eq "yes")
		{
			die "aborting program\n";		
		}
	}
# warn user about disk space for large input files
	&check_filesize($infilename);	

# go thru infile (once only) to create output files
	open (INFILE, "<$infilename") or die "couldn't open $infilename, $!\n";
	my $record_num = 1; # start with 1 instead of zero to match filename output
	
	my $current_seq_obj = "";
	my $current_outfile_name = shift @subfile_names;	
	my %sequence_objects = ();	#key = record number  # value = object reference
	
	open (OUTPUT, ">$current_outfile_name") or die "can't open output file $current_outfile_name";
	
	while (my $line = <INFILE>)
	{
		next if $line =~ /^\s+$/;
		chomp $line;
				
		if ($line =~ /^\>(.+)/)
		{			
		# if there's already a finished sequnce, write it to file	
			if (exists $sequence_objects{$record_num})
			{
				print OUTPUT "$current_seq_obj->{header}\n";
				print OUTPUT "$current_seq_obj->{sequence}\n";
				$record_num++;				
			}
			
		# start a new  object for the current line
			$current_seq_obj =new_seq("record", $line);
			my $obj_num = $record_num;
			$sequence_objects{$record_num} =$current_seq_obj;
		
		# clear previous sequence out of memory
			my $prev_record = $record_num - 1;
			delete $sequence_objects{$prev_record};			
					
		# determine whether a new filehandle needs to be started
			my $current_seq_range = $current_outfile_name;
			$current_seq_range =~ s/$base_name.//;
			my @tmp2 = split "-", $current_seq_range;		
			my $current_file_firstseq = $tmp2[0];
			my $current_file_lastseq = $tmp2[1];
			
			if ($record_num > $current_file_lastseq)
			{
				close OUTPUT;
				$current_outfile_name = shift @subfile_names;
				open (OUTPUT, ">$current_outfile_name") or die "can't open output file $current_outfile_name";				
			}		 
		}				
		else
		{
		# add next line to sequence attribute of current seq object
			$current_seq_obj->{sequence} .= "$line\n";
		}		
	}
	
	# process final sequence in file
		print OUTPUT "$current_seq_obj->{header}\n";
		print OUTPUT "$current_seq_obj->{sequence}\n";

# Feedback & cleanup
	print STDERR "Finished processing $record_num sequences\n"; 
	close OUTPUT;	
		
#################################
# SUBROUTINES
#################################

sub new_seq {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  $self->{header} = $param;
  $self->{sequence} ="";  
  return($self)
}

sub check_filesize
{
	my ($name) = @_;
	my $filesize = -s $name;	
	my $filesize_KB = $filesize/1000;
	my $filesize_MB = $filesize_KB/1000;
	my $filesize_GB = $filesize_MB/1000;
	my $filesize_txt = "$filesize_KB";
		$filesize_txt .= " KB";
	if ($filesize_MB >1)
	{
		$filesize_txt = sprintf("%.1f", $filesize_MB);
		$filesize_txt .= " MB";
	}
	if ($filesize_GB >1)
	{
		$filesize_txt = sprintf("%.1f", $filesize_GB);
		$filesize_txt .= " GB";
	}
	print STDERR "input file size = $filesize_txt\n";
	if ($filesize_GB > 10)
	{
		print STDERR "Are you sure you want to use this much disk space for new files? (y/n)\n";
		my $response = "";
		chomp ($response = <STDIN>);
		unless ($response eq "y" || $response eq "yes")
		{
			die "aborting program\n";		
		}
	}
}

__END__


