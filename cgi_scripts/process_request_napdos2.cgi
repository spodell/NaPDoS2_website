#!/usr/bin/perl
# process_request_napdos2.pl
# Sheila Podell
# December 11, 2019

# manages initial user input from NaPDos (natural product domain seeker) web server
# validates & classifies input sequences for further processing
# assigns job_id, creates job directory
# saves state (job parameters) for delayed retrieval by other scripts
# sends output for short jobs immediately to screen
# creates delayed accession page for long jobs
# modified to disallow queries submitted too close together in time (< minimimum seconds)
# modified 2018 to disallow queries exceeding maximum estimated process time (15 mins)
# modified deny file size max to 30 MB 

use warnings; 
use strict;
use DBI;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use Benchmark;
use File::Basename;

# global variable
	my $gene_max_filesize = 500 * 1024 *1024;	# max size in MB for immediate processing
	my $genome_max_filesize = 50 * 1024 *1024;
	my $deny_filesize = 500* 1024 *1024;;	# 30 MB too big to run online with Vers 1

	$CGI::POST_MAX=1024 * 1024 * 501;  # max 600 MB posts - avoid denial of service attacks
	my $queue_type = "small";
	my $max_seqs = 500000;	# maximum num input sequences
        my $max_time = 54000;   # maximum time in seconds (= 15 mins)
    my $dirsize_limit = "4000000000"; #4GB max temp directory size
    my $daily_diskspace_limit = "10000000000";  # 10 GB for individual user
       my $daily_query_limit = "200";  # number of  queries for individual user 
# server-specific directories
	my $webpath = "https://npdomainseeker.sdsc.edu";	
	my $unixpath = "/misc/www/projects/npdomainseeker";	
	#my $config_filename = "$unixpath/cgi-bin/napdos2019_sdsc.cfg";
	my $config_filename = "$unixpath/cgi-bin/napdos_dev_sdsc.cfg";
		
# relative paths
	my $temp_dir = "$unixpath/tmp";
	my $temp_dir = "$unixpath/tmp2";
	my $logfile = "$unixpath/tmp2/pksdb.log";
	my $custom_page_header = "$unixpath/napdos_templates/napdos2_header";	
	my $custom_page_footer = "$unixpath/napdos_templates/napdos2_page_footer";

# initialize CGI
	my $cgi = CGI->new;
	my $ref_seq_file = $cgi->param('ref_seq_file');
	if (defined $cgi->param('config_filename'))
	{
		my $new_config = $cgi->param('config_filename');
		$config_filename = "$unixpath/cgi-bin/$new_config";
	}	
	
	my $query_type = $cgi->param('query_type') || "pcr"; # alternatives: genome cds aa
		my $compare = "DNA";
		if ($query_type eq "aa")
		{
			$compare = "protein";
		}
	my $domain_type =  $cgi->param('domain_type') || "KS";	# alternative is "C"
	if ($ref_seq_file =~ /C/)
		{			
			$domain_type = "C";
		}	

	my $input = $cgi->param('Sequence');
	my @input_lines = ();
# CGI preserves evil DOS & MAC line endings - must convert these to UNIX for textbox	
	if (defined $input && length $input > 0)
	{
		$input =~ s/\r/\n/sg;
		@input_lines = split/\n/, $input;
	}
	my $INFILE = $cgi->param('seqfile');
	
# Advanced parameters	
	#my $ks_hmm_evalue = $cgi->param('ks_hmm_evalue') || "1e-5";
	#my $ks_min_matchlength = $cgi->param('min_matchlength') || "200";
	#my $c_hmm_evalue = $cgi->param('c_hmm_evalue') || "1e-5";
	#my $c_min_matchlength = $cgi->param('min_matchlength') || "200";
	my $min_matchlength = $cgi->param('min_matchlength') || "200";
	my $path_blast_evalue = $cgi->param('path_blast_evalue') || "1e-8";
	my $max_hits =  $cgi->param('max_hits') || "1";
	#my $tree_hits =  $cgi->param('tree_hits') || "5";
	my $job_id = $cgi->param('job_id');	# won't normally be getting this from CGI
	
# Prepare to write HTML response (site-wide templates)
	my $header = &get_file_text($custom_page_header);	
	$header =~ s/\$webpath/$webpath/g;	
	my $footer = &get_file_text($custom_page_footer);

# Start writing response to screen
	print $cgi->header;
	print $cgi->start_html(-title=>'Search Results',
							-style=>{'src'=>"../css/main.css"});
	print $header;	
				
# Create Job_id to keep track of subsequent processing
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;
	my $dttime_year = sprintf "%04d",($year);  ## four digits to specify the year
	my $dttime_mon  = sprintf "%02d",($mon);      ## zeropad months
	my $dttime_mday = sprintf "%02d",$mday;           ## zeropad day of the month   
	my $current_date = "$dttime_year"."$dttime_mon"."$dttime_mday";
	 unless (defined $job_id && length $job_id > 1)
	 {
	 	$job_id = "$$"."$yday"."$sec";
	 }
	
# track user
	my $start_time = `date`;
	chomp $start_time;
	
	my $user = "unknown_addr";	
	my $remote_host = $ENV{'REMOTE_HOST'};
	my $remote_addr = $ENV{'REMOTE_ADDR'};
	my $remote_user = $ENV{'REMOTE_USER'};
	if (defined $remote_addr)
	{
		$user = "$remote_addr";
	}
	if (defined $remote_host)
	{
		$user .= " $remote_host";
	}	
	if (defined $remote_user)
	{
		$user .= " $remote_user";
	}	
		
# open or create site logfile (with write permissions to webserver)	
	unless (-e $logfile && -w $logfile)
	{
		open (SITE_LOG, ">$logfile") or warn "couldn't open new logfile $logfile\n $!\n";
		print SITE_LOG "New Logfile $start_time\n";
		close SITE_LOG;
		#&send_err("can't write logfile=$logfile<br />\n");
	}
			
	my @log_elements = (
		$user,
		$start_time,
		"job_id: $job_id",
		);
	my $log_printstring = join "\t", @log_elements;

# check most recent user & timestamp - prevent overload
	my $minimum_wait_time = 20; # seconds
	my $num_per_minute = 60/$minimum_wait_time;
	my @prev_log_entry = &get_prev_log_entry();
	my $prev_user = $prev_log_entry[0];
	my $prev_time = $prev_log_entry[1];	
	my $time_diff = &get_time_diff($start_time, $prev_time);
	
	open (SITE_LOG, ">>$logfile") or warn "couldn't open logfile $logfile for appending\n $!\n"; 
	print SITE_LOG "$log_printstring\n";
	
	if ($time_diff > -1  && $time_diff < $minimum_wait_time)
	{
		my $message =  "Time since previous user query = $time_diff seconds. <br />";
		if ($user eq $prev_user)
		{
			$message .= "  Please allow at least one minute between requests. <br />";
		}
		else
		{
			$message .= "  Site is currently too busy to process your request. Please try again later. <br />";
		}
		&send_err ("$message");
	}
	
	
# check available disk space
	my $dirsize = `du -s $temp_dir |cut -f1`;	
	chomp $dirsize;
	my $dirsize_string = &human_readable_filesize($dirsize); # add KB, MB, or GB		
	if ($dirsize > $dirsize_limit)
	{
		#print "tmp file disk usage = $dirsize_string<br />";
		my $message .= "  Site is currently too busy to process your request. Please try again later. <br />";
		&send_err ("$message");
	}
	
# Validate/parse sequence input format (@sequences is list of sequence object references)
	# check input sequence for illegal characters, size exceeding max
	my $input_filesize = 0;	#in kilobytes, assume 1kb per character
	my $size_MB = 0;
		
	if (defined $INFILE && length $INFILE >0)
	{
		&check_filename($INFILE);
		$input_filesize = &check_filesize($INFILE);	
	# getting the lines into an array with this subroutine removes dos or mac line endings by chomp
		@input_lines = &get_seqs_from_file($INFILE);					
	}
	elsif (defined $input)
	{
		my $input_char = length $input;		
		$input_filesize = $input_char;
		if ($input_filesize > $deny_filesize)
		{
			$size_MB = $input_filesize / 1000;
			&send_err("File size of $size_MB MB is too large to process online. 
			Please <a href=/contact.html>contact</a> us if you would like to request a custom analysis.");
		}
	}
	unless ($input_filesize > 0 && scalar @input_lines > 0)
	{
		&send_err("No sequences entered.");
	}
	
# check $user file sizes for current date (note that logfile records file sizesin bytes)
  my $user_diskspace = &tally_user_diskspace($user,$start_time);  
  my $diskspace_string = &human_readable_filesize($user_diskspace);
   		#print "Daily disk space usage for $user = $user_diskspace = $diskspace_string. <br />";
   #my $daily_diskspace_limit = "2000000";  # 20 MB
  if ($user_diskspace > $daily_diskspace_limit)
  {
  	&send_err("Files uploaded in past 24 hours ($diskspace_string) have exceed daily user limit.");
  }
 
# check $user number of daily file entries
	my $user_queries = &tally_user_queries($user,$start_time); 
	#print "num daily queries for user $user = $user_queries <br />";
	if ($user_queries > $daily_query_limit)
	{
  	   &send_err("Number of queries submitted in past 24 hours ($user_queries) have exceed daily user limit.");
  	} 
 
	# get number of fasta sequences
		my $num_input_seqs = 0;
		my $first_seq = "";
	# check first to be sure it's fasta format
		unless (($input_lines[0] =~ /^>.+/) && !($input_lines[1] =~ /^>.+/))
		{
			&send_err("Input data does not appear to be in FASTA format. Please check to be sure first line
			is not blank.");
		}
		foreach (@input_lines)
		{
			if ($_ =~ />/)
			{
				$num_input_seqs++;
				next;
			}	
			if ($num_input_seqs ==1)
			{
				$first_seq .= "$_";
			}	
		}
		if ($num_input_seqs > $max_seqs)
		{
			&send_err("Number of input sequences ($num_input_seqs) exceeds maximum ($max_seqs).
			<br />Please consider our <a href=/napdos2/prefiltering.html> suggestions for pre-filtering</a> 
			your input file before submission.");
			#<br />Please <a href src=../contact.html>contact</a> if you would like to request a custom analysis.");
		}
	my $first_seq_len = length $first_seq || 0;
	
	# could check DNA/protein type and sequence length here as well, to save time
		my $min_nucl_seq_len = $min_matchlength * 3 ; 
		if (($query_type eq "genome") && ($first_seq_len < $min_nucl_seq_len))
		{			
			&send_err ("Input DNA seq length \($first_seq_len nt\) is too 
			short for minimum domain detection align length \($min_matchlength aa\).<br /><br />
			<u>Potential solutions:</u>
			<ol>
			<li> Enable shorter input sequence lengths using  
			<a href = \"https://npdomainseeker.sdsc.edu/napdos2/run_analysis_v2.html#advanced_settings\">Advanced Settings.</a>
			</li>
			<br />
			<li> Select longer sequences by size filtering, using 
			<a href = \"http://npdomainseeker.sdsc.edu//napdos2/prefiltering.html\">File Size Management tools.</a>
			</li>
			</ol>
			<br />");
		}
# Assign a queue type (small=1 cpu, large = 4cpu)
	$queue_type = &assign_queue();
	
# save parameters for later recall
	my @suffixlist = (".cfg");
        my $db_name = fileparse($config_filename, @suffixlist);
        my $napdos_version ="napdos2";
		
	my %saved_params = (
		"ref_seq_file" => "$ref_seq_file",
		"query_type" => "$query_type",
		"first_seq_len" => "$first_seq_len",
		#"ks_hmm_evalue" => "$ks_hmm_evalue",
		#"ks_min_matchlength" => "$ks_min_matchlength",
		"min_matchlength" => "$min_matchlength",
		#"c_hmm_evalue" => "$c_hmm_evalue",
		#"c_min_matchlength" => "$c_min_matchlength",
		"path_blast_evalue" => "$path_blast_evalue",
		"max_hits" => "$max_hits",
		#"tree_hits" => "$tree_hits",
		"job_id" => "$job_id",
		"input_filesize" => "$input_filesize",
		"queue_type" => "$queue_type",
		"user_id" => "$user",
		"start_time" => "$start_time",
		"domain_type" => "$domain_type",
		"num_input_seqs" => "$num_input_seqs",
		"compare" => "$compare",
		"config_filename"  => "$config_filename",
		"db_name"  => "$db_name",
		);
		
# more info to log
	if ($ref_seq_file =~ /C/)
	{
		$domain_type = "C"
	}
	
	$log_printstring = "\tDOMAIN_TYPE=$domain_type\n";
	$log_printstring .= "\tQUERY_TYPE=$query_type\n";
	$log_printstring .= "\tINPUT_FILESIZE=$input_filesize\n";
	$log_printstring .= "\tNUM_INPUT_SEQS=$num_input_seqs\n";
	$log_printstring .= "\tFIRST_SEQ_LEN=$first_seq_len\n";
	$log_printstring .= "\tREF_SEQFILE=$ref_seq_file\n";
	$log_printstring .= "\tDB_NAME=$db_name\n";
	print SITE_LOG "$log_printstring";			
		
# create sub directory in tmp_dir with job_id, open job-specific log file
# BEWARE: umask codes are the opposite of unix command line codes
# (lower numbers have MORE permissions with umask command)		
		
		umask 0001;
		mkdir "$temp_dir/$job_id";
		#print ("The current umask is: ", umask(), "<br />");
		my $job_logfile_name = "$temp_dir/$job_id/$job_id".".log";		
		open (JOBLOG, ">$job_logfile_name") or warn "can't write job logfile\n $!\n";
		print JOBLOG "JOB ID $job_id \t $start_time\n\n";
		print JOBLOG "JOB PARAMETERS\n";
		foreach my $key (sort keys %saved_params)
		{
			unless (defined $saved_params{$key})
			{
				print JOBLOG "undefined parameter $key\n";
				print STDERR "undefined parameter $key\n";
				&send_err("undefined parameter $key\n");
			}			
			print JOBLOG "\t$key=$saved_params{$key}\n";
		}
		close JOBLOG;
				
# open job-specific params file
	my $params_file_name = "$temp_dir/$job_id/$job_id".".params";		
		open (PARAMS, ">$params_file_name") or warn "can't write job params file\n $!\n";
		foreach my $key (sort keys %saved_params)
		{		
			print PARAMS "$key=$saved_params{$key}\n";
		}
		close PARAMS;
		
		my $input_fastaname = "$temp_dir/$job_id/$job_id".".fasta";

	# write user sequences to temp file for later retrieval
		open (FASTA, ">$input_fastaname")or warn "can't write temp fasta\n $!\n";
		foreach (@input_lines)
		{	
			chomp;
			print FASTA "$_\n";
		}		
		close FASTA;
			
	# create access page for delayed retrieval		
		my $results_text = &start_html_page();
		
		$results_text .= $custom_page_header;
		$results_text .= &html_file_output($query_type);
		
		my $results_page = "../tmp2/$job_id/index.html";
		open (RESULTS, ">$results_page") or warn "can't write results page\n $!\n";
		print RESULTS 	"$results_text";
		close RESULTS;
			
# estimate processing time		
	#my $est_time = &get_est_time();	
	my $est_time = &get_est_time_v2();
	#print "$est_time";
		

# User Start processing	
	my $output_html = qq(
		<h2>Job ID $job_id</h2>);
	
	$ output_html .= qq(<p>$est_time </p>);
				
	if ($query_type eq "pcr" || $query_type eq "cds" || $query_type eq "aa")
	 {	 	
	
		if ($queue_type eq "small")
		{
			
			$output_html .= &call_blast_html();			
		}
		
	 }
	 elsif ($query_type eq "genome")
	 {	 	
	 	$output_html .= &call_cand_search_html();
	 }
	 else
	 {
	 	&send_err ("Unrecognized query type: $query_type\n.");
	 }
	 	 		
# Finish HTML response
	print $output_html;
	my $parameter_smry = qq(
	<br/>
	<br/>
	<hr />
	<h3>Search Parameters</h3>
	<ul>
	<li>query type = $query_type</li>
	<li>number of input sequences = $num_input_seqs</li>
	<li>domain type = $domain_type</li>
	);
	
	
	if ($domain_type eq "KS")
	{
		$parameter_smry .=qq(
		<li>minimum match length = $min_matchlength aa</li>
		);
	}
	elsif ($domain_type eq "C")
	{
		$parameter_smry .=qq(
		<li>minimum match length = $min_matchlength aa</li>
		);
	}	
	$parameter_smry .= qq(    
		<li>min e-value for db match = $path_blast_evalue</li>		
		<li>output table matches per query = $max_hits</li>);
		#<li>tree neighbors per query = $tree_hits</li>
	
	my $napdos_version = "NaPDoS v2";
	$parameter_smry .= qq(
		<li>reference sequences = $ref_seq_file</li>
		<li>database version = NaPDoS2_v13b</li>
	
	</ul>);
	
	print "$parameter_smry";

	print "$footer";
	
############################################
# SUBROUTINES
############################################
sub send_err
{
	
	my ($err_msg) = @_;
	print qq(<p>Sorry, your request could not be processed.</p>);
	print "$err_msg";
	if (defined $footer && length $footer > 0)
	{
		print $footer;
	}
	else
	{
		print $cgi->end_html;
	}
	print SITE_LOG "\tERROR: $err_msg\n";
	exit (0);
}

sub start_html_page
{
my $text =qq(<?xml version="1.0" encoding="iso-8859-1"?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US"><head><title>Search Results</title>
<link rel="stylesheet" type="text/css" href="$webpath/css/main.css" />
</head><body>);
return $text;

}

sub html_file_output
{
	#my ($qtype) =@_;
	
	my $pathway_result_page = "$job_id"."_pathways.html";
	my $genomic_locations_page = "$job_id"."_genomic_locations.html";
	my $data_dir = "$webpath/tmp2/$job_id/";

# link to a script that checks to see if job is done or not
	my $text  = qq( <h2>Job ID $job_id</h2>
		Click on the following link to begin analysis.
		<br /><br />				 	
		<form enctype="multipart/form-data" action="$webpath/cgi-bin/access_results_v2.cgi" method="post">
			<input type="hidden" name="job_id"  value=$job_id />
			<input type="hidden" name="query_type"  value="$query_type" />
			<button class = "link_button" type="submit" >JOB ID: $job_id</button>		
		</form>
		<br />
		<p><b>Note</b>: there may be some delay before results for large jobs are available.</p>
		<p>You may bookmark <a href="$webpath/tmp2/$job_id/">this page</a> for data retrieval after the job is complete.</p>
		);
	
# general footer	
	$text .= qq(
	<p></p>
	</div>
	</div> <!--End Main-->
	
	<div id="footer"></div> <!--End Footer-->
	</body>
	</html>			
	);
	
	return $text;
}

sub get_seqs_from_file
{
	# returns an array
	my ($INFILE) = @_;
	my @input_lines = ();	
	my %uniq_lines = ();	# key = line text # value - # times seen 
	my %duplicated_lines = ();	
	while (my $line = <$INFILE>)
	{
	# get rid of DOS line endings	
		$line =~ s/\r/\n/sg;
		chomp $line;
		push (@input_lines, $line);
	# check for illegal characters in header
		if ($line =~ /\t/)
		{
			&send_err("Input file contains one or more tab characters, which must be replaced with spaces.<br />");
		}
		if 	($line =~ /^>.+/)
		{
			if (exists $uniq_lines{$line})
			{
				$duplicated_lines{$line}++;
			}			
			$uniq_lines{$line}++;
		}
	}
	if (scalar keys %duplicated_lines > 0)
		{
			my $msg = "ERROR: Non-unique id numbers detected.<br />";
			foreach my $key (keys %duplicated_lines)
			{
				my $value = $duplicated_lines{$key};
				$value++;	
				$key =~ s/>//;
				$msg .= "&nbsp;&nbsp;&nbsp;&nbsp;$key &nbsp; ($value)<br />";
			}
			&send_err("$msg<br />")
		}			
	return @input_lines;
}

sub check_filesize
{
	my ($INFILE) = @_;
	my @stats = stat($INFILE);
	my $filesize = $stats[7];
	my $filesize_string = $stats[7]/1000000;
	$filesize_string .=" MB";
	my $limit_string = $deny_filesize/1000000;
	$limit_string .=" MB";	
	if ($filesize > $deny_filesize)
		{
			&send_err("Input file size ($filesize_string) exceeds maximum ($limit_string).
			<br />Please consider our <a href=/napdos2/prefiltering.html> suggestions for pre-filtering</a> 
			your input file before submission.");
			#<br />Please <a href src=../contact.html>contact</a> us to request a custom analysis.");
		}
	return $filesize;
}

sub human_readable_filesize
{
	my ($filesize) = @_;
	my $human_readable_size = $filesize; # in bytes
	if ($filesize < 1000)
	{
		$human_readable_size .=" bytes";
	}
	elsif ($filesize > 1000 && $filesize < 1000000)
	{
		
		$human_readable_size = $filesize/1000;
		$human_readable_size .=" KB";
	}	
	elsif ($filesize > 1000000 && $filesize < 1000000000)
	{
		$human_readable_size = $filesize/1000000;
		$human_readable_size .=" MB";
	}
	elsif ($filesize > 1000000000)
	{
		$human_readable_size = $filesize/1000000000;
		$human_readable_size .=" GB";
	}
	return $human_readable_size;
}

sub check_filename
{
	my($name) = @_;
	my $namelength = length $name;
	if ($namelength > 64)
	{		
		&send_err("File name length ($namelength) exceeds max number of characters (64)");
	}
	
	unless ($name =~ /^[_A-Za-z\d#\.\-\/\s]+$/gi)
	{
		$name = "illegal_characters";
		&send_err("Illegal characters in file name.<br> ");
	}
	if ($name =~ /[\.]{2}\// or $name =~ /[\|:?`\>\<\;]/)
	{		
		$name = "illegal_characters";
		send_err("Illegal characters in file name.<br> ");		
	}
}

sub assign_queue
{
	my $queue_type = "small";
	if ($query_type eq "pcr" || $query_type eq "cds" || $query_type eq "aa")
	{
		if ($input_filesize > $gene_max_filesize)
		{
			$queue_type = "large";
		}
		else
		{
			$queue_type = "small";
		}
		
	}
	elsif ($query_type eq "genome")
	{
		if ($input_filesize > $genome_max_filesize)
		{
			$queue_type = "large";
		}
		else
		{
			$queue_type = "small";
		}
	}
	else
	{
		&send_err("undefined query type");
	}
	return $queue_type;
}
	
sub read_config_file
{
    my ($filename) = @_;
    open (INPUT, $filename)  or die "can't open MySQL config file $filename\n$!\n";
	my %dh_config = ();
	while (<INPUT>)
	{
		next if ($_ =~ /^\s+$/);
		next if ($_ =~ /^#/);
		if ($_ =~ /\[(.+)\]=\s*(.+)/)
		{
			$dh_config{$1}=$2;		
		}			
	}
	close INPUT;			
    return %dh_config;
 }  
 
sub get_file_text
{
	# returns a string
	my ($text_sourcefile, $destination) = @_;
	my @text_lines = ();
	 open (INPUT, "$text_sourcefile")  or warn "can't open $text_sourcefile\n$!\n";
	 while (<INPUT>)
	 {
	 	chomp;
	 # convert DOS line endings to UNIX
	 	$_ =~ s/\r/\n/sg;
	 	push @text_lines, $_;
	 }	 
	 close INPUT;
	 my $custom_text = join "\r", @text_lines;
	 return $custom_text;		 	
}

sub call_blast_html
{
	my $text = qq(
		<form enctype="multipart/form-data" action="../cgi-bin/domain_blast2_diamond.cgi" method="post">
			<input type="hidden" name="ref_seq_file"  value=$ref_seq_file />
			<input type="hidden" name="config_filename"  value=$config_filename />
			<input type="hidden" name="max_hits"  value=$max_hits />);
			#<input type="hidden" name="tree_hits"  value=$tree_hits />
		$text .= qq(<input type="hidden" name="compare"  value= $compare />
			<input type="hidden" name="job_id"  value=$job_id />
			<input type="hidden" name="num_input_seqs"  value= $num_input_seqs />
			<input type="hidden" name="min_matchlength"  value=$min_matchlength />
			<input type="hidden" name="path_blast_evalue"  value=$path_blast_evalue />
			<button type="submit" >SUBMIT&nbsp;&nbsp;JOB</button>);	
		return $text;
}

sub call_cand_search_html
{
	my $text = "";
	
	 
	my $immediate_analysis = qq(	 	
		<form enctype="multipart/form-data" action="../cgi-bin/genome_search2_diamond.cgi" method="post">);
		
	$immediate_analysis .= qq(
			<input type="hidden" name="ref_seq_file"  value=$ref_seq_file />
			<input type="hidden" name="config_filename"  value=$config_filename />
			<input type="hidden" name="max_hits"  value=$max_hits />);
			#<input type="hidden" name="tree_hits"  value=$tree_hits />
		$immediate_analysis .= qq(<input type="hidden" name="compare"  value= $compare />
			<input type="hidden" name="job_id"  value=$job_id />
			<input type="hidden" name="num_input_seqs"  value= $num_input_seqs />
			<input type="hidden" name="min_matchlength"  value=$min_matchlength />
			<input type="hidden" name="path_blast_evalue"  value=$path_blast_evalue />
			<button type="submit" >SUBMIT&nbsp;&nbsp;JOB</button>		
		</form>
		);
	 	
	 my $part3 = qq(<p>Note: estimated time to completion: approximately ___ (mins/hours).</p>	
	 	
	 	<hr />

	 	<h3>E-mail Notification</h3>	 	
	 	<p>If you would like to be notified by e-mail when this job is complete, please enter your address below.</p>	 	
	 	<form action="email.cgi"  method="post">
		<p><input type="text" name="uemail" style="width: 50%"></input></p>		
		</form>

	 	<p>&nbsp;</p>
	 	<p>&nbsp;</p>
	 	<p>&nbsp;</p>	 	
	 	<p>$start_time<p>);
	 	
	 $text .= $immediate_analysis;	
	 	
	# $text .= $part3;
	 return $text;
}

sub get_est_time
{
	my $result_text = "";
	
	#my $est_sec = int($input_filesize/40000) + 2; #seconds, based on rough empirical tests (overestimates a little)
# revised estimate based on diamond instead of NCBI blast - 50X faster	
	my $est_sec = int($input_filesize/200000) + 2; #seconds, based on rough empirical tests (overestimates a little)

	# factor in minimum sequence length
	my $seqlength_factor = 1;
	if ($domain_type eq "KS")
	{		
		$seqlength_factor = 200/$min_matchlength;
		
	}
	if ($domain_type eq "C")
	{		
		$est_sec *= 200/$min_matchlength;
	}
	$est_sec = int($est_sec * $seqlength_factor);
	
# convert to minutes for bigger jobs
	my $est_min = $est_sec / 60;
	my $est_time = sprintf ("%.0f", $est_min);
		
	if ($est_sec < 60)
	{
		$result_text = qq(Estimated processing time: ~<b>$est_sec seconds</b>.);
	}
	
	elsif ($est_sec > 60)
	{
				
		$result_text = qq(Estimated processing time: ~<b>$est_time minute);
		if ($est_time > 1)
		{
			$result_text .= qq(s);
		}
		$result_text .= qq(</b>.);
	}
# abort if exceeds max processing time	
	elsif ($est_sec > $max_time)
	{
		my $max_mins = $max_time/60;
		&send_err("Estimated time exceeds maximum web-server allowance of $max_mins minutes.");
	}
	return $result_text;
}

sub get_est_time_v2
{
	my $result_text = "";
	
# revised estimate based on diamond 	
	my $est_sec = 200/$min_matchlength;; #seconds, based on rough empirical tests 
	if ($query_type eq "aa")
	{		
		$est_sec = int(2 + $input_filesize/3000000); # input_filesize is in bytes
		#print "line 778 query_type=$query_type sec = $est_sec input_filesize = $input_filesize num_input_seqs=$num_input_seqs <br />";
	
	# try to predict cases that users have selected for high numbers of positive domains
	# these will take longer to filter non-overlapping hits because have much higher numbers of blast matches 
		if ($num_input_seqs > 100 && $min_matchlength > 100 && $input_filesize < 1000000)
		{
			$est_sec += 20;
		}
		
		elsif ($num_input_seqs > 1000 && $input_filesize < 2000000)
		{
			$est_sec += 20;
		}
	
	}
	elsif ($query_type eq "genome")
	{		
		$est_sec = int(3 + $input_filesize/2000000 + $num_input_seqs/700);
		if ($first_seq_len >1000000) # genome or well-assembled metagenome
		{
			$est_sec += 8;
		}
		if ($num_input_seqs > 10000) # metagenome or PCR
		{
				$est_sec = $est_sec /2.6;
		}
	# try to identify PCR candidates, that could have high percentage of positives 
	# versus genomes and metagenomes, which won't (affects time to get trimmed nucleic coords)
		if ($min_matchlength < 100 
				&& $first_seq_len < 300  
				&& $num_input_seqs > 1000)
 		{
 				#print "looks like PCR<br />";
 				#$est_sec = $est_sec /2.6;
 				$est_sec = $num_input_seqs*0.025;
 		}
	# reduction in time for short PCR sequences
	 	my $size_seq_ratio = ($input_filesize / $num_input_seqs) /1000;
	 	if ($num_input_seqs > 1000 
	 	     && $min_matchlength < 100 
	 	     && $size_seq_ratio < 0.4) 
	 	{
	 		$est_sec =  $est_sec * 0.25;
	 	}	
 		 	
	}
	else 
	{
		print "couldn't determine compare type";
	}	

# convert to minutes for bigger jobs
	my $est_min = $est_sec / 60;
	my $est_time = sprintf ("%.0f", $est_min);
		
	if ($est_sec < 60)
	{
		my $rounded_sec = sprintf ("%.0f", $est_sec);
		$result_text = qq(Estimated processing time: ~<b>$rounded_sec seconds</b>.);
	}
	
	elsif ($est_sec > 60)
	{
				
		$result_text = qq(Estimated processing time: ~<b>$est_time minute);
		if ($est_time > 1)
		{
			$result_text .= qq(s);
		}
		$result_text .= qq(</b>.);
	}
# abort if exceeds max processing time	
	elsif ($est_sec > $max_time)
	{
		my $max_mins = $max_time/60;
		&send_err("Estimated time exceeds maximum web-server allowance of $max_mins minutes.");
	}
	return $result_text;
}
sub get_prev_log_entry
{
	my ($start_time) = @_;
	open (SITE_LOG,"<$logfile") or warn "couldn't open logfile $logfile for reading\n $!\n";
	my $previous_user = "";
	my $previous_start_time = "";
	while (<SITE_LOG>)
	{
		if ($_ =~ /job_id:/)
		{
			my $line = $_;
			my @previous_data = split "\t", $line;
			$previous_user = $previous_data[0];
			$previous_start_time = $previous_data[1];
		}
	}	
	close (SITE_LOG,"<$logfile" );
	my @prev_data = ($previous_user, $previous_start_time);
	return @prev_data;
}

sub tally_user_diskspace # in bytes
{
	my ($current_user_id, $current_timestamp) = @_;
# parse current date
	my @date_terms = split " ", $current_timestamp;
	my $current_date = "$date_terms[0]"." $date_terms[1]"." $date_terms[2]";

	my $tally_user_diskspace = 0;
	my $log_user_id = "$current_user_id";
	my $log_date = "";
	open (LOGFILE,"<$logfile") or warn "couldn't open logfile $logfile for reading\n $!\n";
	while (<LOGFILE>)
	{
		if ($_ =~ /job_id:/)
		{
			my $line = $_;
			my @previous_data = split "\t", $line;
			$log_user_id = $previous_data[0];
			my $log_timestamp = $previous_data[1];
			my @tmp = split " ", $log_timestamp;
			$log_date = "$tmp[0]"." $tmp[1]"." $tmp[2]";
		}
		if ($_ =~ /INPUT_FILESIZE=(\d+)/ 
		&& ("$log_user_id" eq "$current_user_id")
		&& ("$current_date" eq "$log_date" )
		)
		{		
			$tally_user_diskspace += $1;
		}
	}	
	close (LOGFILE);
	return $tally_user_diskspace;
}
sub tally_user_queries 
{
	my ($current_user_id, $current_timestamp) = @_;
# parse current date
	my @date_terms = split " ", $current_timestamp;
	my $current_date = "$date_terms[0]"." $date_terms[1]"." $date_terms[2]";

	my $tally_user_diskspace = 0;
	my $tally_user_queries = 0;
	my $log_user_id = "$current_user_id";
	my $log_date = "";
	open (LOGFILE,"<$logfile") or warn "couldn't open logfile $logfile for reading\n $!\n";
	while (<LOGFILE>)
	{
		if ($_ =~ /job_id:/)
		{
			my $line = $_;
			my @previous_data = split "\t", $line;
			$log_user_id = $previous_data[0];
			my $log_timestamp = $previous_data[1];
			my @tmp = split " ", $log_timestamp;
			$log_date = "$tmp[0]"." $tmp[1]"." $tmp[2]";
			if ( ("$log_user_id" eq "$current_user_id")
			&& ("$current_date" eq "$log_date" ))
			{		
				$tally_user_queries++;
			}
		}
		
	}	
	close (LOGFILE);
	return $tally_user_queries;
}


sub get_time_diff
{
	my ($t1, $t2) = @_;
	unless (defined $t1 && defined $t2)
	{
		return (-1);
	}
		
	my ($t1_weekday, $t1_month, $t1_monthday, $t1_hh_mm_ss, $t1_timezone, $t1_year) = split " ", $t1;
	my ($t2_weekday, $t2_month, $t2_monthday, $t2_hh_mm_ss, $t2_timezone, $t2_year) = split " ", $t2;
	my ($t1_hh, $t1_mm, $t1_ss) = split ":", $t1_hh_mm_ss;
	my ($t2_hh, $t2_mm, $t2_ss) = split ":", $t2_hh_mm_ss;
	
		
	 if (($t1_year == $t2_year) && 
 		($t1_month == $t2_month) && 
 		($t1_monthday == $t2_monthday) && 		
 		($t1_hh == $t2_hh)  && 
 		($t1_mm == $t2_mm) )		
	{
			
		my $sec_diff = $t1_ss - $t2_ss;
		return ($sec_diff);
		
	}
	else 
	{
		return (-1);
	}	
}

__END__
	
# command line debug
	my $cmd_line = 0;
	{
		if ($cmd_line > 0)
		{
			$INFILE = "/Users/spodell/Desktop/PKS_db_2010/test_input_sequences/test_input.faa";
			$ref_seq_file = "$unixpath/fasta/All_KS_seqs_17Mar09.faa";
		}	
	}

# qq(
# 	<h3>Analysis Results</h3>	
# 	 	<p>Your results have been entered into the analysis queue.<br />
# 		There may be some delay before results for large jobs are accessible.</p>);
# 		
		
	#  my $delay_message =qq (	<p><a href ='$results_page'><b><u>Job ID: $job_id</u></b></a></p>
# 	 	<p>&nbsp;</p>);
	 	
	# revise tmp_dir
 		umask 0001;
 		#-rw-rw-r-- 1   28641 w28641 56 Aug 17 12:52 test_text_file2	 	
		my $test_cmd = "whoami";
		$read_test = `$test_cmd`;
		print  "script user = $read_test <br><br>";
		
		$test_cmd = "ls /projects/napdos/tmp2/ ";
                $read_test = `$test_cmd`;
                print "$test_cmd <br>";
                print  "$read_test <br><br>";

                mkdir("/projects/napdos/tmp2/newdir_perlcmd");
                 mkdir("/projects/napdos/tmp2/newdir_perlcmd/subdir");
                $test_cmd = "ls /projects/napdos/tmp2/ ";
                $read_test = `$test_cmd`;
                print "$test_cmd <br>";
                print  "$read_test <br><br>";				
				
