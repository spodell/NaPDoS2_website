#!/usr/bin/perl
# process_align_request_v2019.cgi
# Sheila Podell
# September 10, 2015

# manages user input after protein-level blast search is completed
# sends requests to screen (trimmed sequences or raw blast)
# or sets up another request form (for alignment or tree)
# giving estimated processing time, additional options for output format

use warnings; 
use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);

# global variables
	my $max_hit_num = 50000;
	my $max_aln_num = 100;

# server-specific directories
	my $webpath = "https://npdomainseeker.sdsc.edu";	
	my $unixpath = "/misc/www/projects/npdomainseeker";	
					
# relative paths
	my $temp_dir = "$unixpath/tmp2";
	my $logfile = "$unixpath/tmp/pksdb2.log";
	my $custom_page_header = "$unixpath/napdos_templates/napdos2_header";	
	my $custom_page_footer = "$unixpath/napdos_templates/napdos2_page_footer";
	my $group_membership_table = "$unixpath/templates/group_numbers.tab";	
	
# initialize CGI
	my $cgi = CGI->new;
	my $job_id = $cgi->param('job_id');
		my $data_dir = "$temp_dir/$job_id/";
	my @hit_select_list =  $cgi->param('list');
	my $result_type = $cgi->param('result_type');
	my $align_format = $cgi->param('align_format');
	my $blastout = $cgi->param('blast_output');
	my $align_group = $cgi->param('align_group');
	my $num_matched = $cgi->param('num_matched');
	my $script_source = $cgi->param('script_source');
	# see also section below %saved_params
	
# Prepare to write HTML response (site-wide templates)
	my $header = &get_file_text($custom_page_header);	
# replace placeholder in templates (called $webpath) with webpath from current script
	$header =~ s/\$webpath/$webpath/g;	
	my $footer = &get_file_text($custom_page_footer);

# Start writing response to screen
	print $cgi->header;
	print $cgi->start_html(-title=>'Search Results',
							-style=>{'src'=>"../css/main.css"});
	#wait to print header, so some options can use plain text format

# check to see if need to run rest of script or not	
	my $blast_file = "$data_dir/$job_id".".blastp"; # full path!
	
	if ($result_type eq "raw_blast")
	{
		# check to see if file exists
		my $output = &get_file_text("$blast_file");
		print qq(<pre>$output</pre>);
		print $footer;
		exit(0);
	}
	else
	{						
		#print $napdos_template;
		print $header;
	}

# Get parameter info from job_id file, if available, instead of CGI input
	my $param_info_file = "$data_dir/$job_id.params";
	my %saved_params = &read_params_file("$param_info_file");
		my $config_filename = $saved_params{config_filename};	
		my $query_type = $saved_params{query_type}; #possible values: genome aa 
		my $domain_type = $saved_params{domain_type};
		my $input_filesize= $saved_params{input_filesize};
		my $num_input_seqs= $saved_params{num_input_seqs};
		my $first_seq_len= $saved_params{first_seq_len};
		my $tree_hits= $saved_params{tree_hits} || "5"; # number of tree hits per query
		my $min_matchlength = $saved_params{min_matchlength} || "200";
		my $matchlength = 200;
		my $num_ref_seqs = 0;
		if ($domain_type eq "KS")
		{
			$matchlength = $min_matchlength;
			$num_ref_seqs = 414;
		}
		elsif ($domain_type eq "C")
		{
			$matchlength = $min_matchlength;
			$num_ref_seqs = 172;
		}				

# figure out number of reference sequences
	my $num_queries_selected = scalar @hit_select_list;
	my $num_closest_relatives = $num_queries_selected * $tree_hits;
	
	my $group_membership_text = "";
	my %group_membership_index = ();
	
	if (-s $group_membership_table)
	{
		$group_membership_text = &get_file_text($group_membership_table);
		my @lines = split "\r", $group_membership_text;
		foreach my $line (@lines)
		{
			chomp $line;
			my @tmp = split "\t", $line;			
			$group_membership_index{$tmp[0]} =$tmp[1];
		}
	}	
	
	if (exists $group_membership_index{$align_group})
	{
		$num_ref_seqs = $group_membership_index{$align_group};
	}
	
# Start response depending on request type	
	my $output_html = qq(
		 <h2>Job ID $job_id</h2>
 		);
	 
	my $num_closest = $num_queries_selected * $tree_hits; 		
	my $est_time = &get_est_time($num_queries_selected, $num_ref_seqs);
	my $closest_relative_est_time = &get_est_time($num_queries_selected, $num_closest);
	my $full_db_est_time = &get_est_time($num_queries_selected, $num_ref_seqs);	
 	
 # check number of queries selected
 	unless (scalar $num_queries_selected > 0)
 	{
 		&send_err("No sequences selected");
 	}
 	
	if ($num_queries_selected >$max_hit_num)
	{	
		&send_err(qq(Number of hits found ($num_queries_selected) exceeds maximum 
		allowable using this online server ($max_hit_num). <br /><br />
		It may be necessary to select a domain class subset, <a href ="$webpath/napdos2/prefiltering.html" target="_blank">
		reduce file input size</a> or use <a href ="$webpath/napdos2/run_analysis_v2.html" target="_blank">
		customized parameters</a> to reduce the number of saved matches.<br /><br />	
		Please <a href ="$webpath/contact.html">contact us</a> 
		if you wish to inquire about custom analysis arrangements for large data sets.<br />));
	}
					
# get result type, send appropriate parameter choices for next step	
	my $genomic_search_filename = "$data_dir/$job_id" ."_genomic_locations.html"; 	
	
	if ($result_type eq "genome_location")
	{
		if (-s $genomic_search_filename)
		 {
		 	$output_html = &get_file_text("$genomic_search_filename");
		}
		else
		{
			&send_err ("unable to retrieve stored html file $genomic_search_filename");
		}
	}	
	elsif ($result_type eq "trimmed_nucs")
	{						 	
		$output_html .= &call_align_html(\@hit_select_list);			
		$output_html .= qq(
		<h4>Trimmed sequence format:</h4>);
												
		$output_html .= qq(<table cellpadding = "10" cellspacing = "10">
		<tr>
			<td>
				<input type = "radio" name = "trim_type" value = "amino" checked = "checked" ></input>&nbsp; Amino acid 
			</td>);
		 		
		if ($query_type eq "genome")	
		{	
			$output_html .= qq(<td>
				<input type = "radio" name = "trim_type" value = "nucleic"></input>&nbsp; Nucleic acid 
				</td>);
				$output_html .= qq(</tr>)		
		}
		else
		{
			$output_html .= qq(<td></td></tr>);
		}		
		$output_html .= qq(</table>);	
		$output_html .= qq(<button class = "button" type="submit" >GET TRIMMED SEQS</button>);
		$output_html .=  qq(</form>\n);
		if ($num_matched > 2000 && $query_type eq "genome")
		{
			$output_html .= qq(<br /><p><b>Note:</b> Due to computational limitations, 
			if the initial number of domain matches detected exceeds 2000, <br />
			 nucleotide sequence output will be untrimmed, regardless of the 
			 number selected for further analysis.</p>);
		} 			 		
	}
	elsif ($result_type eq "msf") 
	{	 
		 if ($num_queries_selected > $max_aln_num)
		{
			&send_err ("Number of queries selected ($num_queries_selected) exceeds  
			number of matches for alignment and tree building ($max_aln_num)<br />");
		}
	 	
		 	$output_html .= qq(<p>$est_time </p>);
		 	$output_html .= &call_align_html(\@hit_select_list);
		 	$output_html .= qq(	<button class = "button" type="submit" >GET ALIGNMENT</button>);
		 	$output_html .= qq(		 	
		 	);
		 	$output_html .=  qq(</form>\n);
		 			
		# modify time estimate secs based on # of seqs
	}
	elsif ($result_type eq "newick" || $result_type eq "tree")
	{						
			#$output_html .= qq(<p>$est_time </p>);
		####################
		# Need to insert new section into what was formerly &call_align_html subroutine	
		# otherwise, the align_group doesn't get changed from hidden value passed from original stored params
		
		if ($num_queries_selected > $max_aln_num)
		{
			&send_err ("Number of queries selected ($num_queries_selected) exceeds  
			number of matches for alignment and tree building ($max_aln_num)<br />");
		}
		
		$output_html .=qq(
		<form enctype="multipart/form-data" action="$webpath/cgi-bin/align_seqs2_v2.cgi" method="post">\n
			<input type="hidden" name="job_id"  value="$job_id" />\n
			<input type="hidden" name="result_type"  value="$result_type" />
			<input type="hidden" name="script_source"  value="$script_source" />		
			<input type="hidden" name="num_matched"  value="$num_matched" />
			); 		 
		 if ($align_group =~ /KS/ || $align_group =~ /C/)
		{
		$output_html .=  qq(
		<br />
		<h3>Select Preferred Reference Alignment:</h3>	
		<table cellpadding = "5" cellspacing = "5" border = "0">	
	
		<tr>
			<td>&nbsp;</td>
			<td>
				<input type = "radio" name = "align_group" value = "closest" checked ="checked" ></input>
			</td>				
			<td>
				<b>Closest database matches </b><br /> 
				$num_queries_selected queries + closest database references<br />
				$closest_relative_est_time<br />
			</td>
		</tr>
		<tr>
			<td>&nbsp;</td>
			<td>
				<input type = "radio" name = "align_group" value = "$align_group"  ></input>
			</td>				
			<td>
				<b>$align_group domains </b><br /> 
				 $num_queries_selected queries + all $num_ref_seqs database references<br />
				 $full_db_est_time<br />
			</td>			
	</table>
	<br />
	<br />
		);
		
		}
		else
		{
			$output_html .=  qq(
			<input type="hidden" name="align_group"  value="$align_group" />
			);

		
		}
		
		foreach my $line (@hit_select_list)
		{
			$output_html .= qq(
			<input type="hidden" name = "list" value = "$line" />\n);
		
		}	
	####################

			$output_html .= qq(	<button class = "button" type="submit" >SUBMIT JOB</button>);
		 	$output_html .= qq(
		 	
		 	
		 	);
		 	$output_html .=  qq(</form>\n); 
	}
	else 
	{
		&send_err("Unrecognized request type\n");
	}


	
# Finish HTML response
	
	print $output_html;
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
	#print SITE_LOG "\n\tERROR: $err_msg\n";
	exit (0);
}



sub read_params_file
{
    my ($filename) = @_;
    open (INPUT, "$filename")  or warn "can't open params file $filename\n$!\n";
    
	my %param_input = ();
	while (<INPUT>)
	{
		chomp;		
		next if ($_ =~ /^\s+$/);
		next if ($_ =~ /^#/);
		my @tmp = split "=", $_;
		$param_input{$tmp[0]}=$tmp[1];					
	}
	close INPUT;
	
	unless (scalar keys %param_input > 1)
	{
		&send_err("No saved parameters recognized for job $job_id");
	}
    return %param_input;
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
	 	push @text_lines, $_;
	 }	 
	 close INPUT;
	 my $custom_text = join "\r", @text_lines;
	 return $custom_text;		 	
}

sub call_align_html
{
	my ($arrayref) = @_;
	my @list = @$arrayref;
	my $text = qq(
		<form enctype="multipart/form-data" action="$webpath/cgi-bin/align_seqs2_v2.cgi" method="post">\n
			<input type="hidden" name="job_id"  value="$job_id" />\n
			<input type="hidden" name="align_group"  value="$align_group" />
			<input type="hidden" name="result_type"  value="$result_type" />
			<input type="hidden" name="align_format"  value="$align_format" />
			<input type="hidden" name="script_source"  value="$script_source" />
			<input type="hidden" name="num_matched"  value="$num_matched" />
			);
			#<input type="hidden" name="queryfile"  value=$hmm_cand_file />
			
		foreach my $line (@list)
		{
			$text .= qq(
			<input type="hidden" name = "list" value = "$line" />\n);
		
		}					
		
		return $text;
}

sub get_est_time
{
	my ($num_queries, $num_ref) = @_;
	
	my $result_text = "";
	my $est_sec = 5;
	
	# factor in result type and num sequences to be analyzed
	if ($result_type eq "newick")
	{
		$est_sec = $est_sec + 10 + ($num_queries+ $num_ref)*0.25;
	}
	elsif ($result_type eq "msa")
	{
		$est_sec += $num_queries *2;	
	}
	
	#factor in minimum sequence length
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
	return $result_text;
}
__END__
