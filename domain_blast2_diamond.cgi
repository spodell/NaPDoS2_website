#!/usr/bin/perl
# domain_blast2_diamond.cgi
# Sheila Podell
# February 1, 2021

#	Takes in CGI input FASTA sequence(s) 
# 	Accepts sequence from files as well as textbox (files over-ride text) 
#	Parses, verifies fasta format of input
#	Cuts off if filesize larger than 500,000 sequences 
	# $limit, set in get_seqs_from_file subroutine )

# 	Runs BLAST search as unix system call to command lline
#	Returns results as html page in choice of formats

use warnings; 
use strict;
use DBI;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use Benchmark;
use Scalar::Util qw(looks_like_number);

# global variables
	my $num_processors = 2;

# initialize CGI 
	use CGI qw(:standard);
	use CGI::Carp qw(fatalsToBrowser);
	my $cgi = CGI->new;
	my $ref_seq_file = $cgi->param('ref_seq_file');
	my $max_hits =  $cgi->param('max_hits');
	my $compare = $cgi->param('compare'); #"protein" or "DNA"
	# this isn't giving correct information for parent file type - always says "protein
	my $evalue = $cgi->param('evalue') || "1e-8";
	my $num_nuc_matches = $cgi->param('num_nuc_matches') || "0";
	my $out_format = $cgi->param('output') || "table";	# "table" or ??
	my @input_lines = ();
	my $num_input_seqs = $cgi->param('num_input_seqs') || 0;
	my $job_id = $cgi->param('job_id');
	my $config_filename = $cgi->param('config_filename');
	
# check $max_hits value to be sure it's a number 
	unless (looks_like_number($max_hits)) 
	{
  		&send_err("Illegal characters in max_blast_hits variable:<br> $max_hits ");
	} 
	
	if ($max_hits =~ /\&/)
	{
		&send_err("Illegal characters in max_blast_hits variable:<br> $max_hits ");
	}

# figure out domain type, if not explicitly sent
	my $domain_type ="KS";
	if ($ref_seq_file =~ /C/)
		{			
			$domain_type = "C";
		}
	
# server-specific directories
	my $webpath = "https://npdomainseeker.sdsc.edu";
	my $unixpath = "/misc/www/projects/npdomainseeker/napdos2/";	
	my $blast_dir = "$unixpath/bin/blast-2.2.18_linux-x64/bin";
	my $diamond_dir = "$unixpath/bin/diamond-0.9.29/";		
	my $config_filename = "$unixpath/cgi-bin/napdos_dev_sdsc.cfg";	# should get this from calling script or web page
	my $transeq_dir = "$unixpath/bin/EMBOSS-6.6.0/emboss/";
	#my $hmm_prog_dir = "$unixpath/bin/hmmer-3.2/bin/";	
	
# relative paths
	my $pks_script_dir = "$unixpath/scripts/";
	my $database_fasta = "$unixpath/fasta/$ref_seq_file";
	my $temp_dir = "/projects/napdos/tmp2/";
	my $data_dir = "$temp_dir/$job_id/";
	my $template_dir = "$webpath/templates/";
	my $webpath_pathway_template_dir = "$webpath/pathway_templates2/";
	my $unixpath_pathway_template_dir = "$unixpath/pathway_templates2/";
	my $site_log = "$temp_dir/pksdb.log";
	my $custom_page_header = "$unixpath/napdos_templates/napdos2_header";	
	my $custom_page_footer = "$unixpath/napdos_templates/napdos2_page_footer";

# Prepare to write HTML response
	my $header = &get_file_text($custom_page_header);	
	$header =~ s/\$webpath/$webpath/g;	
	my $footer = &get_file_text($custom_page_footer);
	my $pathway_result_page = "$job_id"."_pathways.html";
	my $search_detail_page = "$job_id"."_search_results.html";
	
# Start writing HTML response
	print $cgi->header;
	print $cgi->start_html(-title=>'Search Results',
							-style=>{'src'=>"$webpath/css/main.css"});
	
# print to file as well as screen
	open (PATHWAYS_SMRY_HTML, ">$data_dir/$pathway_result_page") or warn "couldn't open $pathway_result_page for writing\n $!\n";
	#my $pathways_html_file_text = &start_html_file();
	#print PATHWAYS_SMRY_HTML "$pathways_html_file_text";
	
	open (SEARCH_DETAIL_OUTPUT, ">$data_dir/$search_detail_page") or warn "couldn't open $search_detail_page for writing\n $!\n";
	my $search_html_file_text = &start_html_file();
			
# open logfiles to append messages
	open (SITE_LOG, ">>$site_log") or warn "couldn't open logfile $site_log for appending\n $!\n";

# check_filename exits if finds illegal characters (|?`> : ; < ..)
	my $INFILE = $cgi->param('infile') || "$data_dir/$job_id".".fasta";	
	
	if (length $INFILE > 2)
	{
		&check_filename($INFILE);
		@input_lines = &get_seqs_from_file($INFILE);
	}
	
# check to make sure all required accessory programs are available 
	&check_programs;
				
# Validate/parse input format (@sequences is list of sequence object references)
	my $seqref = \@input_lines;
	my @sequences = &get_seqs($seqref);		
	my $sequence_type = &check_sequence_type(\@sequences);
	my $program = "blastp"; #default
	
# make sure input sequence is appropriate (nt vs aa)
	if ($compare eq "protein")
	{
		unless ($sequence_type eq "protein")
		{
			&send_err("Data input problem.<br> Fasta format amino acid sequence required for protein comparison.");
		}
		$program = "blastp";
	}
	if ($compare eq "DNA")
	{
		unless ($sequence_type eq "nucleic")
		{
			&send_err("Data input problem.<br> Fasta format nucleic acid sequence required for DNA comparison.");
		}
		$program = "blastx";
	}

	unless (-s $database_fasta)
	{		
		&send_err("Couldn't find blast database reference file $ref_seq_file.<br> Please try clearing your browser cache.");
	}
	
# Get parameter info from job_id file, if available, instead of CGI input
	my $param_info_file = "$data_dir/$job_id.params";
	my %saved_params = ();
	my $path_blast_evalue = "$evalue";
	my $min_matchlength = 200;
	unless (-s "$param_info_file")
	{
		&send_err("Can't find input data for JOB $job_id $param_info_file");
	}
	
		%saved_params = &read_params_file("$param_info_file");
		my $query_type = $saved_params{query_type}; #genome or aa
		my $domain_type = $saved_params{domain_type};
		$ref_seq_file = $saved_params{ref_seq_file};
		$config_filename = $saved_params{config_filename};
		$max_hits = $saved_params{max_hits} || 1;
		$evalue = $saved_params{evalue} || "1e-8";
		$job_id = $saved_params{job_id};
		$min_matchlength = $saved_params{min_matchlength} || 200;
		$path_blast_evalue = $saved_params{path_blast_evalue} || "1e-8";
		$evalue = $path_blast_evalue;
	
# write user input into PROTEIN fasta (either orig aa or nucleotide ORFs) 
# exclude sequences that don't meet minimum size requirements
	my $pos_filename = "$data_dir/$job_id.fa";
	my %query_hit_tally = ();	# key = query sequence_id  # value = num maches
	open(FASTA_POS, ">$pos_filename") or die "can't open temp file $pos_filename for writing.\n $!";
	
	foreach (@sequences)
	{
		my $seqlen = length ($_->{sequence});
		#print "annot=$_->{annotation} seqlen=$seqlen<br>";
		next unless ($seqlen > $min_matchlength);		
		$query_hit_tally{$_->{ID}} = 0;
		print FASTA_POS "$_->{annotation}\n";
		print FASTA_POS "$_->{sequence}\n";
	}	
	close FASTA_POS;
	
# Get genome domain info file, if exists
	my $filename = "$data_dir/$job_id"."_genomic_locations.tab";
	my $genome_filetext = &get_file_text($filename);

	my %genomic_cand_info = ();	# key = cand_id #value = tab-delimited text from Domain locations table

	if (defined $filename && -s $filename)
	{
		%genomic_cand_info = &read_genome_tabfile($filename); # key = domain number ("domain_id), value = rest of line
	}
	# note: this file will not exist for large PCR data sets
	
# system call to blast to get output, keep track of time
	my $t_blast_start = new Benchmark;
	my $blast_format = "";
	
	# this isn't used with diamond
	if ($out_format eq "table")
	{
		$blast_format = " -m8 ";
	}
		
# diamond db conversion 
	my $diamond_db = $database_fasta;
	$diamond_db =~ s/.faa/.dmnd/;
	my $blastp_result_filename = "$job_id".".blastp";

# run a blastp search of amino acid sequences agains database		
## save enough matches in raw blastp to support later tree building (default = 5)
	# diamond default max is 25 (not enough for some sequences that are well represented in db)
	# --max-target-seqs 0 gets ALL matches meeting align criteria
	# for well conserved sequences, might need 400+ matches per domain to get them all
	my $blastp_raw_filename = "$blastp_result_filename".".raw";
	my $blast_output = "";	
   if (-s "$data_dir/$blastp_raw_filename") # if this was already done for a genome
	#if ((-s $genome_filetext) || -s "$data_dir/$blastp_raw_filename") # if this was already done for a genome
	{
		$blast_output = &get_file_text("$data_dir/$blastp_raw_filename");
	#	print "found pre-existing raw blastp output file $data_dir/$blastp_raw_filename<br />";
	}	
	else 
	{
		my $blastp_command2 = "$diamond_dir/diamond $program -d $diamond_db -q $pos_filename -e $evalue --max-target-seqs 0 -p $num_processors --quiet -o $data_dir/$blastp_raw_filename";	
		#print "$blastp_command2 <br />";
		$blast_output = `$blastp_command2`;
	}
	
	my $num_lines_raw_blast_output = `wc -l $data_dir/$blastp_raw_filename` || 0;
	chomp $num_lines_raw_blast_output;	
	
# filter out lines with alignments that don't meet minimum length requirement
	my $blastp_filtered_filename = "$job_id"."_filtered.blastp";
	my $filter_blastp_cmd = "$pks_script_dir/select_tab_quant_greater_than_v2.pl $data_dir/$blastp_raw_filename 3 $min_matchlength > $data_dir/$blastp_filtered_filename";  
  #print "$filter_blastp_cmd<br />";
   `$filter_blastp_cmd`;
	`mv -f $data_dir/$blastp_filtered_filename $data_dir/$blastp_raw_filename`;

# get tophits with non-overlapping coordinates
# because simple tophits won't get all valid hits with long sequences with > 1 domain if limit # of matches 		
	my $num_lines_filtered_blast_output = 0;
	my $overlap_distance = 10; # num of residues for alignment "slop" to detect a single domain
	my $max_mismatch = "200";
	if ($num_lines_raw_blast_output > 0)
	{			
		# don't re-do filtering if was already done by genome_search2_diamond.pl
		if ((-s $genome_filetext) && -s ("$data_dir/$blastp_result_filename")) # if this was already done for a genome
		{
			$blast_output = &get_file_text("$data_dir/$blastp_result_filename");
			#print "found pre-existing filtered blastp output file $data_dir/$blastp_result_filename<br />";
		}
		else
		{
			#print "Line 259 num_lines_raw_blast_output=$num_lines_raw_blast_output<br>
			#unable to find raw blast output $data_dir/$blastp_result_filename <br>";
			#my $tophits_cmd = "$pks_script_dir/select_tab_instances.pl $data_dir/$blastp_raw_filename $max_hits >$data_dir/$blastp_result_filename";
			#$blast_output = `$tophits_cmd`;
		}
				
		unless ($query_type eq "genome") #don't overwrite results from genome_search2_diamond.pl 
		#unless ((-s $genome_filetext) && -s "$data_dir/$blastp_result_filename") #don't overwrite results from genome_search2_diamond.pl 
		{			
			my $tophits_cmd = `$pks_script_dir/select_non-overlapping_tophits.pl $data_dir/$blastp_raw_filename $overlap_distance >$data_dir/$blastp_result_filename`;						
		# debug - works at command line
			my $debug = "$pks_script_dir/select_non-overlapping_tophits.pl $data_dir/$blastp_raw_filename $overlap_distance >$data_dir/$blastp_result_filename";
			#print "$debug<br>";
		}
		
	# In-place removal of any extra blank lines
	# This screws up all C-domain candidates - don't do it.
	# note: `sed -i '/^$/d' $filename works from command line but not from	CGI script
		#my $blank_line_removal_cmd = qq(perl -i -n -e "print if /S/" $data_dir/$blastp_result_filename);
		#my $debug2 = qq(perl -i -n -e "print if /S/" $data_dir/$blastp_result_filename);
		#print "$debug2<br>";		 
		#my $result = `$blank_line_removal_cmd`;
		
		$num_lines_filtered_blast_output = `wc -l $data_dir/$blastp_result_filename` || 0;
		chomp $num_lines_filtered_blast_output;
	}	
	
# don't use command line output from diamond directly, because it has extra lines 
	# diamond v0.9.29.130 | by Benjamin Buchfink <buchfink@gmail.com>
	# Licensed under the GNU GPL <https://www.gnu.org/licenses/gpl.txt>
	# Check http://github.com/bbuchfink/diamond for updates.
	# here is the solution:
	#my $blank_line_removal = `sed -i '/^$/d' $data_dir/$blastp_result_filename`;

	my $t_blast_finish = new Benchmark;
	my $blast_time = timediff ($t_blast_finish, $t_blast_start);
	print SITE_LOG "\tDOMAIN BLAST SEARCH diamond:", timestr($blast_time, 'nop'), "\n";
		
# error message if no hits
	if ($num_lines_filtered_blast_output < 1) 
	{		
			print $cgi->start_html;
			print qq(<br><font size = "+1"><b>No matches found</font></b><br /><br />);
			print "Search parameters:<br>";
			print "&nbsp;&nbsp;&nbsp;domain type: $domain_type<br />";
			print "&nbsp;&nbsp;&nbsp;comparison type: $compare<br />";
			print "&nbsp;&nbsp;&nbsp;input sequence type: $sequence_type<br />";
			print "&nbsp;&nbsp;&nbsp;blastall program: $program<br />";
			print "&nbsp;&nbsp;&nbsp;minimum match length: $min_matchlength<br />";
			print "&nbsp;&nbsp;&nbsp;e-value cutoff: $evalue<br /><br />";
			print "<p>You may want to try modifying search parameters 
			to allow shorter minimum alignmnent lengths and/or a less stringent 
			e-value cutoff.</p>";									
		#debug
			#print "$blast_command<br>";	
	}
	else
	{			
		$blast_output = &blastout_from_file("$data_dir/$blastp_result_filename");
			
	}
		
	my @blast_headers = ("Query_id",
	"Subject_id", "% identity", "align_length", "mismatches", "gap_openings",
	"q_start", "q_end", "s_start", "s_end", "e-value", "bit_score");
my $blast_header_line = join "\t", @blast_headers;

# get query and subject id numbers. Count number of blast matches per query
	my @blast_lines = split /\n/, $blast_output;
	my @query_ids = ();
	my @match_ids = ();	
	
	foreach my $line (@blast_lines)
	{
		chomp $line;
		my @tmp = split '\t', $line;
		push @query_ids, $tmp[0];
		push @match_ids, $tmp[1];
		$query_hit_tally{$tmp[0]}++;
	}

# get info on queries that had no blast matches
	my $total_num_input = scalar (keys %query_hit_tally);
	my @unmatched_ids =();
	foreach (keys %query_hit_tally)
	{
		if ($query_hit_tally{$_}< 1)
		{
			push @unmatched_ids, $_;
		}
	}
	my $num_unmatched = scalar @unmatched_ids;
	
	my $num_matched = $total_num_input - $num_unmatched;
	#my $filename = "$data_dir/$job_id"."_genomic_locations.tab";
	
# could get this from KS domains list for genome
	# 316073017_blastp_cand_coords this is the correct number of matches
		
	#my $genome_filetext = &get_file_text($filename);
	#my $num_matched = `wc -l `;
	# num_matched = grep -c ">" 316073017_trimmed_blast_KSdomain_cands.faa		
	#my $stat_string = qq(BLAST matches found for $num_matched/$total_num_input input sequences.);	
		
# connect to mySQL database (connection parameters from config file)
    my %config_params= &read_config_file($config_filename);
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_port = $config_params{db_port};
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";	
	
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	if (exists $config_params{db_port} && defined $db_port)
	{
		$db_path .= "$db_port";
	}
	my @row = ();
	my $element;
	my $firstline;
	my $sql;
	my $sth;
	
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password);
	unless (defined $dbh)
	{
		&send_err("Couldn't connect to database");
	}		

# construct SQL query to get pathnames 
	$firstline = "select domains.domain_name, path_name, domain_type, domain_class, domain_subclass
	from pathways, domains 
	where pathways.path_id = domains.path_id
	and domains.domain_name in (";	
	my @result_table=&sql_select_list(\@match_ids, $firstline);
					
# convert db info to hash table, key=domain name, value = rest of tab-delimited string 
	my %domain_pathways;
	foreach (@result_table)
	{
		chomp;
		my @tmp = split '\t', $_;
		
		my $key =shift @tmp;
		$domain_pathways{$key} = \@tmp;
	}
	
# construct SQL query to get all group names and num members for current domain being searched
	# note that group_name consists of class, subclass, separated by "!" character	
	$firstline = qq(select group_name, count(memberships.domain_name) as num_members
from domains, memberships
where domains.domain_name = memberships.domain_name
and domain_type = "$domain_type"
group by group_name
order by num_members desc);	
	my @group_info = &sql_select_simple($firstline);	
	
# construct SQL query to get group membership tally for matches

# this needs to be separately from pathways tally, otherwise will fail to join
	# if there is no group associated with match
	# but currently (10/14/20) there are no groupnames being found
	
	$firstline = "select group_name, domains.domain_name 
		from domains, memberships 
		where domains.domain_name = memberships.domain_name
		and domains.domain_name in (";	
	my @matched_groups_table=&sql_select_list(\@match_ids, $firstline);
# convert to hash
	my %matched_groups_table = ();	# key = groupname # value = num hits in current data set
	foreach (@matched_groups_table)
	{
		chomp;
		my @tmp = split '\t', $_;
		my $key =$tmp[0];
		$matched_groups_table{$key}++;
	}
	
# pass search tab results to display function to get html page		
	my $tab_results_table = &display_tab_results($blast_output, \%domain_pathways, \%matched_groups_table, \@group_info);	
	
	my $tabout_filename = "$job_id"."db_search_results.tab";
	&text_to_file($tab_results_table, "$temp_dir/$job_id/$tabout_filename");
	
#	print tally info to screen
	my $subclass_tally_html_filename = "$job_id"."subclass_tally.html";  # duplicates PATHWAYS_SMRY_HTML
	my $subclass_tally_html = &start_html_file();
	$subclass_tally_html .= &check_form_javascript();
	$subclass_tally_html .= &subclass_tally("$temp_dir/$job_id/$tabout_filename");
	print PATHWAYS_SMRY_HTML "$subclass_tally_html\n";
	&text_to_file($subclass_tally_html, "$temp_dir/$job_id/$subclass_tally_html_filename");
	print $subclass_tally_html;
	
# problem: subclass_tally_html includes an extra copy of header, already put into CGI output
# and this is the wrong header!
		 
#	store search results for future display
	my $display = &display_search_results_html($blast_output, \%domain_pathways, \%matched_groups_table, \@group_info);	                                

	my $search_html_header = &start_html_file(); 
	$search_html_file_text .= "$display";
	$search_html_file_text .= "$footer";	
	print SEARCH_DETAIL_OUTPUT $search_html_file_text;	

# Finish HTML response (custom ending)
	{
		print "$footer";
	}
	
# clean up
	print $cgi->end_html;
	$sth->finish();
	$dbh->disconnect();
	
	close SITE_LOG;
	#close JOB_LOG;
	close PATHWAYS_SMRY_HTML;
	close SEARCH_DETAIL_OUTPUT; 
	
############################################
# SUBROUTINES
############################################

sub send_err
{	
	my ($err_msg) = @_;
	print $cgi->start_html;
	print qq(<p>Sorry, your request could not be processed.</p>);
	print $err_msg;	
	
	#print JOB_LOG "\nERROR: $err_msg\n";
	print SITE_LOG "\n\tERROR: $err_msg\n";
	
	print $cgi->end_html;
	close SITE_LOG;
	close JOB_LOG;
	close PATHWAYS_SMRY_HTML;
	unlink "$data_dir/$pathway_result_page";
	exit (0);
}

sub get_seqs_from_file
{
	my ($INFILE) = @_;
	&check_filename($INFILE);
			
		
	my @input_lines = ();
	my @stats = stat($INFILE);
	my $filesize = $stats[7];
	my $filesize_string = $stats[7]/1000000;
	$filesize_string .=" MB";
	my $limit = 500000000; #150 MB
	my $limit_string = $limit/1000000;
	$limit_string .=" MB";
	
	if ($filesize > $limit)
		{
			send_err("Input file size ($filesize_string) domain_blast2 exceeds maximum ($limit_string).");
		}
	open (FILEHANDLE, "$INFILE") or warn "can't open debug filehandle\n";	
	while (my $line = <FILEHANDLE>)
	{
		chomp $line;		
		push (@input_lines, $line);
	}			
	close FILEHANDLE;
	return @input_lines;
}

sub check_filename
{
	my($name) = @_;
	unless ($name =~ /^[_A-Za-z\d#\.\-\/\s]+$/gi)
	{
		send_err("Illegal characters in file name:<br> $name ");
	}
	if ($name =~ /[\.]{2}\// or $name =~ /[\|:?`\>\<\;]/)
	{		
		send_err("Illegal characters ( $& ) in file name:<br> $name ");		
	}
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
 
sub check_programs
{	
	my @program_list = (
		"$diamond_dir/diamond",	
		"$transeq_dir/transeq",
		"$pks_script_dir/transeq_to_multifasta.pl",
		"$pks_script_dir/getseq_multiple.pl",
		"$pks_script_dir/get_subsequence.pl",
		"$pks_script_dir/select_tab_quant_greater_than_v2.pl",
		"$pks_script_dir/select_tab_instances.pl",
		"$pks_script_dir/select_non-overlapping_tophits.pl",
	);	
	
	foreach my $prog(@program_list)
	{
		unless (-s "$prog" && -x "$prog")
		{
			&send_err("ERROR: can't find accessory program $prog.");
			exit(0);
		}
	}							
}
 
sub read_genome_tabfile
{
    my ($filename) = @_;
    open (INPUT, $filename)  or warn "can't genome domain table $filename\n$!\n";
	my %tabfile_lines = ();	# key = candidate id # value = rest of line
	while (<INPUT>)
	{
		next if ($_ =~ /^\s+$/);
		next if ($_ =~ /^#/);
		next if ($_ =~ /frameshift/);
		chomp;
		
		my @tmp = split "\t", $_;
		my $cand_id = $tmp[6];
		unless (defined $cand_id && length $cand_id >0)
		{
			warn "can't find cand_id name in read_genome_tabfile, line $.";
			next;
		}
		if ($cand_id =~ /(.+)\*$/)
		{
			$cand_id = $1;
		}
		
		$tabfile_lines{$cand_id} = $_;
					
	}
	close INPUT;			
    return %tabfile_lines;
 }

sub get_seqs 
{
	my ($seqfile_name) = @_;
	my @seqlines = @$seqfile_name;

	my $counter = -1;
	my @seqs = ();	#a list of sequence object references
	my $limit = 500000;
	
	if (scalar @seqlines <2)
	{
		send_err("no sequences entered");
	}
			
	foreach my $line (@seqlines)
	{
		chomp $line;
		next if ($line =~ /^\s+$/);		
		if ($line =~ /^>.+/) 
		{					
			chomp $line;
			$counter++;
			if ($counter > $limit)
			{
				send_err("too many sequences ($counter) - exceeds limit");			
			}
			my %seq_obj = new_seq();
			my $seq_ref = \%seq_obj;			
			push (@seqs, $seq_ref);
			$seqs[$counter]->{annotation} = $line;
			my @tmp = split " ", $line;
			$tmp[0] =~ s/>//;			
			$seqs[$counter]->{ID} = $tmp[0];
			next;								
		}
		if ($counter > -1)
		{			
			$seqs[$counter]->{sequence} .= uc"$line";
			$seqs[$counter]->{sequence} =~ s/ //g; #remove any internal spaces
			$seqs[$counter]->{sequence} =~ s/\r//g; #remove any internal linefeeds
			$seqs[$counter]->{sequence} =~ s/\f//g; #remove any internal carriage returns
		}					
	}
		
	unless (scalar @seqs > 0)
	{
		my $err = "Problem with input file format<br>";
		$err .= "First line: $seqlines[0]<br>";
		$err .= "<br>If uploading files using Internet Explorer 5.x on a Macintosh,"; 
		$err .= " try a different browser.";
		send_err($err);
	}	
	
	if ($counter < 0)
		{
			send_err("no sequences found");
		}
	return @seqs;
}

sub new_seq 
{
	my %sequence = (
		annotation => "",
		sequence => "",
	);	
	return %sequence;
}

sub check_sequence_type
{
	my ($seqlist_ref) = @_;
	my @seqlist = @$seqlist_ref;
	my ($seqstring, $seqtype);
	foreach (@seqlist)
	{
		$seqstring .= $_->{sequence};
	}
	
	my $num_ATCGN = (uc$seqstring =~ tr/ATUCGN//);
		my $num_aminos = (uc$seqstring =~ tr/DEFHIKLMPQRSVXY//);
	
	if  ($num_ATCGN eq length $seqstring)
	{
		$seqtype = "nucleic";
		return $seqtype;
	}			
	if ($num_aminos > 1)
	{
		$seqtype = "protein";
		return $seqtype;
	}
	
	return $seqtype;
}

sub blastout_from_file
{
	my ($blastfile) = @_;
	my $text = "";
		
# debug
	#print "blast_result_file=$blastfile<br />";			
	# print lines			
	open(BLASTOUT, "<$blastfile") or die "can't open diamond blast infile reading.\n $!";
	
	while(<BLASTOUT>)
	{
		chomp;
		$text .= "$_\n";
	}
	
	close BLASTOUT;
	return $text;
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

sub text_to_file
{
	my ($text, $outfilename) = @_;
	open(TEXTOUT, ">$outfilename") or die "can't open temp file $outfilename for writing.\n $!";

# sort in numerical order (by 1st column, if there are domain id numbers)
#	my @textlines = split "\n", $text;	
#	my $header = shift @textlines;
#	my @sorted = sort {$a <=> $b} @textlines;
#	my $sorted_text = join "\n", ($header, @sorted);
#	print TEXTOUT "$sorted_text\n";
	
	print TEXTOUT "$text\n";
	close TEXTOUT;
}

sub construct_sql_list
{
	my ($array_ref) = @_;
	my @array = @$array_ref;
	my $outstring = "";
	my %uniq = ();
	
	foreach (@array)
	{		
		unless (exists $uniq{$_})
		{
			$uniq{$_}++;
			$outstring .= qq('$_' ,);
		}
	}	
	# remove terminal comma
	$outstring = substr($outstring, 0, -1);
	return $outstring;	
}

sub sql_select_simple 
{
	# takes simple select statement, returns array of tab-delimited lines
	
	my ($firstline) = @_;
	my $sql = $firstline;	
	my @results_tbl = ();
	
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	my @row = ();
	while (@row = $sth->fetchrow_array)
	{				
		my $row_txt;
		foreach my $element (@row)
		{
			$row_txt .= "$element\t";
		}
		push @results_tbl, "$row_txt";
	}
	
	#print "sql_select_simple:<br />   $sql<br />";	
	return @results_tbl;

}
sub sql_select_list
{		
	my ($list, $firstline) = @_;
	my @select_list = @$list;
	my $total_num = scalar @select_list;
	my $sql = $firstline;
	
	my @results_tbl = ();
		
	unless ($total_num >0)
	{
		print STDERR "couldn't find any selections in list.\n";
		exit(0);
	}
	my $count = 0;			
	foreach (@select_list)
	{
		$count++;
		chomp;
		$_ =~ s/'//g;		# get rid of internal quotes
		$sql .= "'$_',";
		
	# send large select statements in pieces to avoid exceeding MySql max_allowed_packet size	
		if ($count % $max_lines_per_packet == 0)
		{
			$sql = substr($sql, 0, -1); # get rid of extra newline and comma, if present
			$sql .= ");\n";
						
			$sth = $dbh->prepare($sql)
    			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
			$sth->execute()
    			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			@row = ();
			while (@row = $sth->fetchrow_array)
			{				
				my $row_txt;
				foreach my $element (@row)
				{
					$row_txt .= "$element\t";
					#print "<br />line 866 next sql out=$row_txt<br />";
				}
				push @results_tbl, "$row_txt";
			}
			$sql = $firstline;
		}
		elsif ($count == $total_num)
		{
			$sql = substr($sql, 0, -1); # get rid of extra newline and comma, if present
			$sql .= ");\n";
			
			 $sth = $dbh->prepare($sql)
     			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
 			$sth->execute()
     			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			@row = ();
			while (@row = $sth->fetchrow_array)
			{
		 		 my $row_txt = join "\t", @row;
				push @results_tbl, "$row_txt";
			}
		}		
		
	}		
	return @results_tbl;
}


sub start_html_file
{
my $webpath = "https://npdomainseeker.sdsc.edu";

# start html, general page  header
my $text =qq(<?xml version="1.0" encoding="iso-8859-1"?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US"><head><title>Search Results</title>
<link rel="stylesheet"  href="/css/main.css" type="text/css" media="screen"  />
</head><body>
<script src="../scripts/sorttable.js" type="text/javascript"></script>);

# why is this script not showing up?
# general napdos page  header, links relative to nested tmp/job_id directory
	$text .= qq(
	<div id="header" class="centered">
<a href="$webpath/napdos_home_v2.html"><img src="$webpath/image/NaPDoS_blue_2.png" alt = "NaPDoS logo" /></a>
<div id="subtitle">Natural Product Domain Seeker</div>
<div id="topmenu" class="centered">
<ul>
<li><a href="$webpath/napdos2/napdos_home_v2.html">Home</a></li>
<li><a href="$webpath/napdos2/quickstart.html">QuickStart</a></li>
<li><a href="$webpath/napdos2/run_analysis_v2.html">Run Analysis</a></li>
<li><a href="$webpath/napdos2/classification.html">Classification</a></li>
<li><a href="$webpath/napdos2/pathways_v2.html">BGCs</a></li>
<li><a href="$webpath/napdos2/contact_v2.html">Contact Us</a></li>
</ul>
<div class="clear"></div>
</div> <!--End Topmenu-->
</div><!--End Header-->
<div id="main" class="centered shadowmain">
<div id="maininner">);
	
	return $text;
}
sub display_tab_results
{
	#($blast_output, \%domain_pathways, \%matched_groups_table, $stat_string, \@group_info);
	my ($output, $hashref1, $hashref2, $arrayref) = @_;
	my @outlines = split /\n/, $output; 	
	
	my %db_fields = %$hashref1;	#key =pathname  #value= arrayref of variables (domain type, domain subtype 
	my %found_groups = %$hashref2;
 	my $tab_text = "";
	my @col_names = ("","cand_id","database match","percent identity", "align length","e-value");
	my @metadata_fields = ("BGC match", "domain class", "domain subclass");

	#my $num_metadata_fields = scalar @metadata_fields;		
	@col_names = (@col_names, @metadata_fields);
	
	foreach my $colname (@col_names)
	{
		$tab_text .= "$colname\t";
	}
	$tab_text .= "\n";
	
	#debug
		#$tab_text .= "query type=$query_type \n";
		#$tab_text .= "compare type=$compare \n";

		foreach my $line (@outlines)
		{			
			my @cols = split /\t/, $line;
			my $queryname = "$cols[0]"; 
			my $matchname = "$cols[1]"; #orf name - might not be unique
			my $pct_id = sprintf ("%.0f", "$cols[2]");			
			my $align_len = "$cols[3]";
			my $qstart = "$cols[6]";
			my $qend = "$cols[7]";
			my $evalue = "$cols[10]";
			my $domain_name = join "_", ($queryname,$qstart,$qend);
			$queryname = $domain_name;
			
			#my $query_match_pair = "$queryname".","."$matchname";
			my $blast_lines_arrayref = \@cols;
			# my $multidomain_queryname = $queryname; # contains orf number AND protein align coordinates
			# debug
				#$tab_text .= "domain_blast line 931 queryname=$queryname \n";
			 # if ($query_type eq "genome") 
#  			{
# 				$multidomain_queryname = join "_", ($queryname,$qstart,$qend); 
#  				$queryname = $multidomain_queryname;
#  				#$tab_text .= "domain_blast line 935 queryname=$queryname \n";
#  			}
			
		# check to see if there is any domain id tag from genome information table	
			my $domain_id = " ";
			if (defined $genomic_cand_info{$queryname})
			{
				my @tmp = split "\t", $genomic_cand_info{$queryname};
				$domain_id = $tmp[1];	
			}						
			my $textline = join "\t", ($domain_id,$queryname,$matchname,$pct_id,$align_len,$evalue);
			$tab_text .= "$textline\t";
		
		# collect group and pathway info, add to blast results	
			my $pathname = "";
			my $path_info = "";
			my ($arrayref, $domain_name,$class, $subclass);
			my @metadata_properties = ();
			if (exists $db_fields{$matchname})
			{
				$arrayref = $db_fields{$matchname};
				@metadata_properties = @$arrayref;
				$pathname = $metadata_properties[0];	
				$domain_name = $metadata_properties[1];
				$class = $metadata_properties[2];
				$subclass = $metadata_properties[3];
				if ($pathname =~ /unknown_.+/)
				{
					$pathname = "unknown";
				}
				$path_info = lc ("$pathname"."_info.html");
				unless (-s "$unixpath_pathway_template_dir/$path_info")
				{
					$path_info = "no_path_info".".html";
				}
			}
			else
			{
				$pathname = "no data";
				$path_info = "no_path_info".".html";
			}			
			
			$tab_text .= "$pathname\t$class\t$subclass\n";
		}
		
		my @textlines = split "\n", $tab_text;	
		my $header = shift @textlines;
		my @sorted = sort {$a <=> $b} @textlines;
		my $sorted_text = join "\n", ($header, @sorted);
		return $sorted_text;		
	}
sub display_search_results_html
{
	my ($output, $hashref1, $hashref2, $arrayref) = @_;
	my @outlines = split /\n/, $output;
	my $webpath = "https://npdomainseeker.sdsc.edu";
# if query type = genomic, sort outlines 	
	my %db_fields = %$hashref1;	#key =pathname  #value= arrayref of variables (domain type, domain sutype 
	my %found_groups = %$hashref2;
	my $v1_webpath ="https://npdomainseeker.sdsc.edu";
		
# fix in case there is more than one match per input sequence id	
	my $num_matchlines = scalar @outlines;  # problem if there is a blank line at end
	#print "display_search_results num_matchlines=$num_matchlines<br />";
	# my $non_blank_matchlines = 0;
# 	foreach (@outlines)
# 	{
# 		chomp;
# 		next if ($_ =~ /^\s+$/);
# 		$non_blank_matchlines++;
# 	}
# 	$num_matchlines = $non_blank_matchlines;
	
	$num_matched = $num_matchlines;
	#print "line 998 num_matched=$num_matched<br />";
	#print "num_matched =$num_matched<br />";
 
 my $html_text2 .= &check_form_javascript();
 
 	$html_text2 .= qq(<h2>Database Search Results</h2>
	$num_matchlines $domain_type domains were identified from $num_input_seqs input sequences.<br />);
		
	$html_text2 .= qq(<p>Use check boxes to select candidates for <b><a href="#further_analysis">further analysis</a></b> below. 
	<br />(options may take a few seconds to load for large match tables)
	Click on column headers to sort (multiple clicks toggle between ascending and descending order).</p>);

if ($num_matched > 100)
	{
			$html_text2 .= qq(<br />Note: further analysis options will not appear until table is fully loaded 
			(may take a while for large tables). 
			<br />Maximum number of matches for alignment and tree building is 100.);
	}
	else 
	{
			$html_text2 .= qq(Click on column headers to sort (multiple clicks toggle between ascending and descending order).);
	}

		
$html_text2 .= qq(<form  name="align_fasta" action="$webpath/cgi-bin/process_align_request_v2019.cgi" method="post">		
	<input type = "hidden" name = "job_id" value = "$job_id" />
	<input type = "hidden" name = "num_matched" value = "$num_matched" />
	<input type = "hidden" name = "script_source" value = "allseq" />);
	
	$html_text2 .= qq(<p><input class="button" type="button" value="Select All" onclick="this.value=check(this.form.list)" /></p>);	
	
$html_text2 .= qq(<!-- Start table of Individual sequence matches -->);	
if ($num_matched <500)
{
	$html_text2 .= qq(<script src="../../scripts/sorttable.js" type="text/javascript"></script>
	<table class = "sortable2" width = "100%" cellpadding = "5" cellspacing = "0" border = "1">
	<tr>);
}
else 
{
	$html_text2 .= qq(	
	<table class = "sortable2" width = "88%" cellpadding = "5" cellspacing = "0" border = "1">
	<tr>);
}
	 	
		my @col_names = ("","cand_id","database match","percent identity", "align length","e-value");			
		my @metadata_fields = ("BGC match", "domain class", "domain subclass",);

		@col_names = (@col_names, @metadata_fields);
		
		foreach my $colname (@col_names)
		{
			if (length $colname <1)
			{
				$html_text2 .= qq(<th class="unsortable">$colname</th>\n);
			}
			else
			{
				$html_text2 .= qq(<th>$colname</th>\n);
			}
		}
	$html_text2 .= qq(
		</tr>);	
		

# fill in table with data
		foreach my $line (@outlines)
		{			
		next if ($line =~ /database/i);
		
		$html_text2 .= qq(<tr>);		
			
		# $html_text .=subset of blast columns
			my @cols = split /\t/, $line;
			my $queryname = "$cols[0]";
			my $matchname = "$cols[1]";
			my $pct_id = sprintf ("%.0f", "$cols[2]");			
			my $align_len = "$cols[3]";
			my $qstart = "$cols[6]";
			my $qend = "$cols[7]";
			my $evalue = "$cols[10]";
			my $blast_lines_arrayref = \@cols;

		
		# wrap long names, so table doesn't get too wide	
			if (length "$queryname" > 20)
			{
			$html_text2 .= qq(<style>
			td {
					word-break: break-all;
			}
			</style>
			);
			}
			
		# check to see if there is any domain id tag from genome information table	
			my $domain_id = " ";
			if (defined $genomic_cand_info{$queryname})
			{
				my @tmp = split "\t", $genomic_cand_info{$queryname};
				$domain_id = $tmp[1];			
			}
			
		# fix for displaying multidomain orfs (works for genomes, but not aa queries	
			my $multidomain_queryname = join "_", ($queryname,$qstart,$qend);
			$queryname = $multidomain_queryname;
			
			$html_text2 .=qq(  <td align = "center"><input type = "checkbox" name = "list" value = "$line"></input></td>\n);			
			$html_text2 .=qq(  <td align = "left">$queryname);
				if (defined $domain_id && length $domain_id >1)
				{
					$html_text2 .= " [$domain_id]";
				}
			
			$html_text2 .=qq(</td>\n);		
			$html_text2 .=qq(  <td align = "left">$matchname</td>\n);
			$html_text2 .=qq(  <td align = "center">$pct_id</td>\n);
			$html_text2 .=qq(  <td align = "center">$align_len</td>\n);
			$html_text2 .=qq(  <td align = "center">$evalue</td>\n);				
		
		# collect group and pathway info, add to blast results	
			my $pathname = "";
			my $path_info = "";
			my $domain_class_ref = ""; #for URL link to background page
			my ($arrayref, $domain_name,$class, $subclass);
			my @metadata_properties = ();
			if (exists $db_fields{$matchname})
			{
				$arrayref = $db_fields{$matchname};
				@metadata_properties = @$arrayref;
				$pathname = $metadata_properties[0];
				$domain_name = $metadata_properties[1];
				$class = $metadata_properties[2];
				$domain_class_ref = $class;
					$domain_class_ref =~ s/ /_/g;
				$subclass = $metadata_properties[3];
				if ($pathname =~ /unknown_.+/)
				{
					$pathname = "unknown";
				}
				$path_info = ("$pathname"."_info.html");
				unless (-s "$unixpath_pathway_template_dir/$path_info")
				{
					$path_info = "no_path_info".".html";
				}
			}
			else
			{
				$pathname = "no data";
				$path_info = "no_path_info".".html";
			}			
			$html_text2 .=qq(  <td align = "left"><a href ="$webpath_pathway_template_dir/$path_info" target="_blank">$pathname</a></td>\n);
			$html_text2 .=qq(  <td align ="center">
				<a href="https://npdomainseeker.sdsc.edu/napdos2/background_v2.html#$domain_class_ref" target="_blank">$class</a></td>\n);
			$html_text2 .=qq(  <td align = "center">$subclass</td>\n);
			$html_text2 .=qq(</tr>\n);			
		}

# print rest of html form		
		$html_text2 .=qq(	
		</table>
		<!-- End Table Individual sequence matches -->	
	<p>	</p>
	<br />\n);
	

# offer opportunity to download tab-delimited file
	my $download_path = "$webpath/tmp2/$job_id/$tabout_filename";
	$html_text2 .= qq(
	<b>Right-click to <a href="$download_path" target="_blank">DOWNLOAD</a></b> this table in tab-delimited format. 
	<br />);
	
	
$html_text2 .=qq(

	<a name = "further_analysis"></a>	
			
	<h2>Options</h2>
	
<!-- Start Table 2 analysis options -->		
	<table cellpadding = "5" cellspacing = "5" border = "0">);
	
	my $genomic_search_filename = "$data_dir/$job_id" ."_genomic_locations.html";
	my $genome_loc_tabfile = "$data_dir/$job_id"."_genomic_locations.tab"; 	
	if (-s $genomic_search_filename && -s $genome_loc_tabfile)
	{
		$html_text2 .=qq(<tr>
		<td>&nbsp;</td>
		<td>
			<input type = "radio" name = "result_type" value = "genome_location" checked ="checked" ></input>
		</td>				
		<td>
			<b>View nucleotide coordinates </b>for all trimmed domain candidates.<br />
		</td>
	</tr>);	
	}
	
	$html_text2 .=qq(
	<tr>
		<td>&nbsp;</td>
		<td>
			<input type = "radio" name = "result_type" value = "trimmed_nucs" checked ="checked" ></input>
		</td>				
		<td>
			<b>Output selected sequences </b>in fasta format.<br />
		</td>
	</tr>);
	
# don't try to align if number is too high 
if ($num_matched > 100)
{
	$html_text2 .=qq(</table>
	<!-- End Table 2 analysis options -->
	<table width = "90%"> 		
 	<tr>
 		<td>
 			<p>
				<button type="submit" >GET RESULTS</button>
			</p>			
		</td>
 	</tr>
	</table>);
	return $html_text2;
}	
	
$html_text2 .=qq(	<tr>
		<td>&nbsp;</td>
		<td>
			<p><input type = "radio" name = "result_type" value = "msf"></input></p>
			<br /> 
		</td> 
		<td><b>Output Alignment </b> with closest database matches <br /> 
			Select alignment format: 
			<select name="align_format">
				<option value="fasta"> fasta </option>
				<option value="msf"> msf</option>		
				<option value="clw"> clustalw </option>				
			</select>
		</td>
		<td>&nbsp;</td>
	</tr>
	
	<tr>
		<td>&nbsp;</td>
		<td>
			<input type = "radio" name = "result_type" value = "newick" ></input>
		</td>			
		<td>
			<b>Construct tree </b>(candidate domains + blast matches + reference domains)<br />
		</td>
	</tr>
	
	</table>
	<!-- End Table 2 analysis options -->);



 $html_text2 .= qq(<table width = "90%"> 		
 	<tr>
 		<td>
 			<p>
				<button type="submit" >GET RESULTS</button>
			</p>			
		</td>
 	</tr>

	</table>);
		
# Assign pre-align group as all available KS or all available C domains
if ($domain_type eq "KS")
{
	$html_text2 .= qq(
 <input type = "hidden" name = "align_group" value = "all_classified_KS" ></input>
  );
 }
elsif ($domain_type eq "C")
{
	$html_text2 .= qq(
	<input type = "hidden" name = "align_group" value = "all_classified_C" ></input>
	);
}
else
{
	&send_err("Couldn't find domain_type in domain_blast2.cgi");
}
  
  
$html_text2 .= qq( </form>			
	);
	return $html_text2;
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
    return %param_input;
 } 

sub read_genome_tabfile2
# not used
{
	my ($filename) = @_;
    open (INPUT, $filename)  or warn "can't genome domain table $filename\n$!\n";
	my %tabfile_lines = ();	# key = candidate id # value = rest of line
	my %candidate_objects = (); # key = candidate id # value = object reference
	my $current;
	while (<INPUT>)
	{
		# skip header lines 
		next if ($_ =~ /seq/);
		next if ($_ =~ /^\s+$/);
		next if ($_ =~ /^#/);
		next if ($_ =~ /frameshift/);
		
		chomp;
		my @tmp = split "\t", $_;
		my $cand_id = $tmp[6];
		my $domain_id = $tmp[1];
		my $line_num = $domain_id;
	# enable later sorting by domain number	
		$line_num =~ s/KS//;
		$line_num =~ s/C//;
		unless (defined $cand_id && length $cand_id >0)
		{
			warn "can't find cand_id name in read_genome_tabfile, line $.";
			next;
		}
	# get rid of asterisk for frame shift
		if ($cand_id =~ /(.+)\*$/)
		{
			$cand_id = $1;
		}
		$tabfile_lines{$cand_id} = $_;
		my $line_text = $_;
		$current = &new_genomic_cand("domain", $cand_id,$line_num, $line_text, $domain_id);
		$candidate_objects{$cand_id} = $current;					
	}
	close INPUT;			
    return %candidate_objects;
}


sub new_genomic_cand
{
	my ($className, $cand_id, $line_num, $line_text) = @_;
 	my $self = {};
  	bless $self, $className;
  	 	
   $self->{cand_id} = $cand_id;
   $self->{line_num} = $line_num;
   $self->{line_text} = $line_text;
  	return($self);
}	
 	
sub subclass_tally
{
	my ($filename) = @_; #all match info
	my %group_tally = ();

# go thru infile, tally
    open (INPUT, $filename)  or warn "can't find genome db_hitfile $filename\n$!\n";
	my %group_tally = (); # key= "$class!$subclass", # value = tally
	my %group_input_ids = (); # key= "$class!$subclass", # value = newline separated string of id numbers
	my @outstring_headers = ("Class", "Subclass", "Num_matches");
	my $outstring = join "\t",@outstring_headers; # tab-delimited version of subclass tally
		$outstring .= "\n";
	my $html_text = ""; 
	my $cand_orfs_filename = "$job_id"."_cand_orfs.faa";
	my $search_details_filename ="$job_id"."_search_results.html";
	my $genomic_locations_html = "$job_id"."_genomic_locations.html";
	my $combined_results_pkg = "$job_id"."_combined_results.gz";
	my $matched_nucleotide_seqfile = "$job_id"."_parent_match_seqs.fna";
	my $subclass_tally_table = "$job_id"."subclass_tally.tab";
	my $webpath = "https://npdomainseeker.sdsc.edu";

# make sure $num_matched has been corrected if >1 match per input sequence	
	my $matchline_count = `wc -l $data_dir/$blastp_result_filename` || 0;
	
	if ($matchline_count =~ /(\d+)\s.+/)
	{
		$matchline_count = $1;
	}
	else 
	{
		$matchline_count = 0;
	}
	$num_matched = $matchline_count;
	 
# start html table
	$html_text .= qq(<h2>Domain Classification Summary</h2>);
	if ($min_matchlength < 200)
	{
		$html_text .= qq(
		<b>WARNING:</b><i>
		Minimum alignment length has been set below recommended threshold for reliable
		classification.<br /> Please interpret results with caution.</i>
		<br/>);
	}

# deal with unfixable differences in diamond output for how gaps are handled in
# blastx (nucleic acid query) versus blastp (trimmed orfs or aa query)
	$html_text .=	qq(<ul>
		<li>$num_matched $domain_type domains were identified from $num_input_seqs input sequences.</li>);
		# if (($num_unmatched >0) && (-s "$data_dir/$genomic_locations_html")) 
# 		{
# 			$html_text .=	qq(<li><p>NOTE: preliminary numbers of borderline-quality domain candidates identified in nucleic acid 
# 			sequences may not always agree exactly with those found in translated ORFs, due to differences in the way e-values are 
# 			calculated for the two data types and/or potential frame shifts or stop codons revealed by translation.</p></li>);	
# 		}
	# explain options	
		$html_text .=	qq(
		<li><p>Click on buttons below to view a detailed table of results, download domain sequences 
		in fasta format,<br /> 
		or build comparative trees.</p></li>
		</ul>
			
		<a href="$webpath/tmp2/$job_id/$search_details_filename" target="_blank">
		<button><b>VIEW ALL MATCHES</b></button></a><br  />);
        
	$html_text .= qq(<br /><hr />);

	# Table 2 - detailed info on individual classes
	$html_text .= qq(
	<form  name="subseq_search_smry" action="$webpath/cgi-bin/subseq_search_smry.cgi" method="post">
	<input type = "hidden" name = "job_id" value = "$job_id" ></input>
	<input type = "hidden" name = "compare" value = "$compare" ></input>
	<h3>Individual Domain Classes</h3>
	<ul>
	<li>Select one or more categories below to view a subset of matches.</li>
	</ul>);
	
	# add a note here about genomic sequences being untrimmed 
	if ((-s "$data_dir/$genomic_locations_html") && ($num_matched > 500))
		{
			# $html_text .= qq(Note that downloadable nucleotide sequences will be untrimmed if the total
# 			number of matches exceeds 500.<br />);
	   };
	
# don't want to show "Select All" button if there is only one group - because it doesn't work in this case	
	# problem: don't know how many groups there are yet.
		#my $num_groups = scalar (keys %group_input_ids);
		#$html_text .= qq( num_groups =	$num_groups<br />);
	#if ($num_groups >1)
	#{
		$html_text .= qq(<p><input class="button" type="button" value="Select All" onclick="this.value=check(this.form.list)"></input></p>);	
	#}
	
	$html_text .= qq(<table class = "sortable2" width = "50%" cellpadding = "5" cellspacing = "0" border = "1">
	 <tr>
 		<th style = "text-align: center"> Select </th>
		<th> Class </th>
 		<th> Subclass </th>
		<th style = "text-align: center"> Num <br />matches </th>
 	</tr>);
	
	# get groupname info
	while (my $line = <INPUT>)
	{
		chomp $line;
		next if ($line =~ /^\s+$/);
		next if ($line =~ /query/i); # skip header
		next if ($line =~ /database/i);
		my @tmp = split "\t", $line;
		next unless (defined $tmp[6] && defined $tmp[7]);
		my $class = $tmp[7] || "class_not_found";
		my $subclass = $tmp[8] || "subclass_not_found";
		my $groupname = "$class"."!"."$subclass";
		$group_tally{$groupname}++;
		$group_input_ids{$groupname} .= "$tmp[1]"."\n";
	}
		
	# create output in order of descending abundance 
	foreach my $key (sort {$group_tally{$b} <=> $group_tally{$a}} keys %group_tally)
	{	
		my @tmp = split "!", $key;
		#next if ($key =~ /domain/i); # skip header
		my $class = $tmp[0];
		next if ($class =~ /domain/i); # skip header
		my $domain_class_ref = $class;
			$domain_class_ref =~ s/ /_/g;
		my $subclass = $tmp[1];
		my $count = $group_tally{$key};
		my $id_list = $group_input_ids{$key};	
		my @seqlist = split "\n", $id_list;
		
		my $nospaces_txt = $key;
		$nospaces_txt =~ s/ /_/g;
		$nospaces_txt =~ s/!/-/g;
		my $subclass_id_num_file = "$job_id"."$nospaces_txt"."_ids";
		&text_to_file($id_list,"$data_dir/$subclass_id_num_file");
		
	# add class stats to tabfile	
		my @outlines = ($class, $subclass, $count);
			my $nextline = join "\t", @outlines;
		$outstring .= "$nextline\n";
	
	# add class stats to html output	
		$html_text .= qq(
		<tr>
		<td style = "text-align: center"><input type = "checkbox" name="list" value= "$nospaces_txt" ></input></td>
		<td><a href="$webpath/napdos2/classification.html#$domain_class_ref" target="_blank">$class</a></td>
		<td>$subclass</td>
		<td style = "text-align: center">$count</td>			
		</tr>	
		);	
	}
# end html
		$html_text .=qq(
		</table>
		<p>Right-click to <a href = "$webpath/tmp2/$job_id/$subclass_tally_table"><b>DOWNLOAD</b></a> this table in tab-delimited format.</p> 				
		<button type="submit" >VIEW A SUBSET</button>
		</form>);
		
		
# send html tab outstring to tabfile
		&text_to_file($outstring,"$data_dir/$subclass_tally_table");
		
		return $html_text;
}

sub check_form_javascript
{
	my $text = qq(<script type="text/javascript">
		<!-- Hide script from old browsers --
			function check_all(theForm) 
			{
		  		for (i=0; i<theForm.elements.length; i++) 
		  		{
					theForm.elements[i].checked = true;
		  		}
			}
			function check_none(theForm) 
			{
		  		for (i=0; i<theForm.elements.length; i++) {
				theForm.elements[i].checked = false; return false"> &nbsp;
		  }
		}
	// End of hiding -->
		</script>

			<script type="text/javascript" language="JavaScript">
									
			<!-- Begin
			var checkflag = "false";
			function check(field) {
			if (checkflag == "false") {
			for (i = 0; i < field.length; i++) {
			field[i].checked = true;}
			checkflag = "true";
			return "Unselect All"; }
			else {
			for (i = 0; i < field.length; i++) {
			field[i].checked = false; }
			checkflag = "false";
			return "Select All"; }
			}
			//  End -->
			</script>);
	return $text;
}		
__END__

my @fields = (
	"Query id",
	"Subject id",
	"pct identity",
	"align length",
	"mis-matches",
	"gap opens",
	"q_start",
	"q_end",
	"s_start",
	"s_end",
	"e-value",
	"bit-score",);
	my $text = join "\t", @fields;
	
Blast fields
0	Query_id
1	Subject_id
2	% identity
3	align_length
4	mismatches
5	gap_openings
6	q_start
7	q_end
8	s_start
9	s_end
10	e-value
11	bit_score

