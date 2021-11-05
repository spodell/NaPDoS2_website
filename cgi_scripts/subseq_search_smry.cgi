#!/usr/bin/perl
# subseq_search_smry.cgi
# Sheila Podell
# September 21, 2020

# takes as input:
	# a directory where previously NaPDoS calculated data is stored
	# within that directory, a list of previously identified query id numbers
	
	# gets header + selected lines from tabfiles based on id numbers
		# search results.tab
		# genome locations.tab
	# writes header + selected lines into new tabfiles
	# outputs a new detailed search results html file, with options for:
		# downloading tab-delimited files
		# retrieving trimmed fasta sequences (protein and/or nucleotide)
			# if > 500 matches, gets untrimmed only
		 
use warnings; 
use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use Benchmark;

# global variables	
	my $max_domain_hits = 2000;  # prevent server overload for blastx searches (domain coordinates)

# initialize CGI 
	use CGI qw(:standard);
	use CGI::Carp qw(fatalsToBrowser);
	my $cgi = CGI->new;
	my $ref_seq_file = $cgi->param('ref_seq_file');
	my $max_hits =  $cgi->param('max_hits');
	my $compare = $cgi->param('compare'); #"protein" or "DNA"
	my $out_format = $cgi->param('output') || "table";	# "table" or ??
	my $num_input_seqs = $cgi->param('num_input_seqs') || 0;
	my $job_id = $cgi->param('job_id');
	my @group_names = $cgi->param('list');
			
# figure out domain type, if not explicitly sent
	my $domain_type ="KS";
	if ($ref_seq_file =~ /C/)
		{			
			$domain_type = "C";
		}
		
# server-specific directories
	my $webpath = "https://npdomainseeker.sdsc.edu";
	my $unixpath = "/misc/www/projects/npdomainseeker/napdos2/";		
		
# relative paths
	my $pks_script_dir = "$unixpath/scripts/";
	my $database_fasta = "$unixpath/fasta/$ref_seq_file";
	my $temp_dir = "/projects/napdos/tmp2/";
	my $data_dir = "$temp_dir/$job_id/";
	my $template_dir = "$webpath/templates/";
	my $webpath_pathway_template_dir = "$webpath/pathway_templates2/";
	my $unixpath_pathway_template_dir = "$unixpath/pathway_templates2/";
	my $site_log = "$temp_dir/pksdb.log";
	my $custom_page_footer = "$unixpath/napdos_templates/napdos2_page_footer";	
	
# Start writing HTML response
	print $cgi->header;
	print $cgi->start_html(-title=>'Search Results',
							-style=>{'src'=>"/css/main.css"});
	my $display = &start_html_file();	
	
# Get parameter info from job_id file, if available, instead of CGI input
	my $param_info_file = "$data_dir/$job_id.params";
	my $query_type = "";
	my %saved_params = ();
	if (-s "$param_info_file")
	{
		%saved_params = &read_params_file("$param_info_file");
		$query_type = $saved_params{query_type};
		$ref_seq_file = $saved_params{ref_seq_file};
		$max_hits = $saved_params{max_hits} || 1;
		$job_id = $saved_params{job_id};	
	}
	else
	{
		#&send_err("Can't find input data for JOB $job_id $param_info_file");
	}
	
# open logfiles to append messages
	open (SITE_LOG, ">>$site_log") or warn "couldn't open logfile $site_log for appending\n $!\n";
	
# check to make sure all required accessory programs are available 
	&check_programs;
	
# identify files (that should already be) written by previous CGI scripts
	# amino acid data
		my $orig_aa_fasta_filename = "$data_dir/$job_id".".fa"; #4216427049.fa
		my $aa_blastp_tophits_tab_filename = "$data_dir/$job_id".".blastp"; #4216427049.blastp
		my $aa_blastp_raw_tab_filename = "$data_dir/$job_id".".blastp.raw"; #4216427049.blastp.raw
		my $aa_search_results_tab_filename ="$data_dir/$job_id"."db_search_results.tab"; #4216427049db_search_results.tab
		my $trim_aa_coords_filename = "$data_dir/$job_id"."_trim_aa_coords"; # does not pre-exist
		my $trimmed_aa_domains_fasta_filename = "$data_dir/$job_id"."_trim_aa.faa"; # does not pre-exist
		
	# nucleic acid data
		my $orig_nuc_filename = "$data_dir/$job_id"."_parent_match_seqs.fna"; #3056227934_parent_match_seqs.fna
			
		my $translated_orf_filename = "$data_dir/$job_id"."_cand_orfs.faa"; #2879427045_cand_orfs.faa		
		#my $trim_nuc_coords_filename = "$data_dir/$job_id"."_trim_cands_diamond.coords"; #2879427045_trim_cands_diamond.coords # has both orig & orfname
			# $trim_nuc_coords_filename will not exist for large data sets, match_num >2000
		my $selected_nuc_coords_filename = "$data_dir/$job_id"."_selected_trim_nuc_coords"; #will not exist for large data sets, match_num >2000
		my $nuc_blastx_filename	= "$data_dir/$job_id"."_trim_cands_diamond.blastx"; #will not exist for large data sets, match_num >2000
		my $trim_nuc_domains_filename = "$data_dir/$job_id"."_trimmed_nucleotide_"."$domain_type"."domain_cands.fna"; #2879427045_trimmed_nucleotide_KSdomain_cands.fna
		my $selected_trim_nuc_domains_filename = "$data_dir/$job_id"."_selected_nucleotide_domain_cands.fna"; # does not pre-exist	
		my $selected_id_filename = "$data_dir/$job_id"."_selected_ids"; # for large pcr files - no coordinates
		my $selected_domain_id_filename = "$data_dir/$job_id"."_selected_domain_ids";
	
# get rid of previous versions of selected coords and trimmed files within same job, if exist 
	if (-e "$selected_id_filename") #"$data_dir/$job_id"."_align_select_ids";
		{
			unlink $selected_id_filename ;
		}
	if (-e "$selected_domain_id_filename") #"$data_dir/$job_id"."_align_select_ids";
		{
			unlink $selected_domain_id_filename ;
		}	 

# get all id number files into a single array
	my $num_groupnames = scalar @group_names;
	unless ($num_groupnames >0)
	{
		&send_err ("no subgroups selected<br>");
	}
	my @id_array = ();
	my %id_hash = (); # key = domain number, value = count
	
	# get domain id numbers for all sequences that belong in each group
	foreach my $next_name (@group_names)
	{
		#print "subseq_search line 123 <br /> next groupname = $next_name <br />";
		my $id_filename = "$data_dir/$job_id"."$next_name"."_ids";
		my $id_txt = &get_file_text($id_filename);
		my @tmp = split "\n", $id_txt;
		foreach my $next_id (@tmp)
		{
			push  @id_array, $next_id;
			$id_hash{$next_id} = 1;		
		}
		
		unless (-s $id_filename)
		{
			#$display .= "couldn't find id file $id_filename <br />";
			#&send_err ("couldn't find id file $id_filename");
		}
	}

	my $num_matched = scalar @id_array;

# write id numbers to tmp file
	my $select_id_nums_filename = "$data_dir/$job_id"."_select_id_nums";
	my $select_id_nums_txt = join "\n", @id_array;
	$select_id_nums_txt .= "\n";  # make sure there is a newline at end 
	&text_to_file($select_id_nums_txt, $select_id_nums_filename);

# get trimmed aa fasta sequence files for the subset
	my $allseq_faa_filename = "$data_dir/$job_id".".fa";
	my $trim_aa_coords = "$data_dir/$job_id"."_select_aa_coords";
	my $trim_aa_seq_filename = "$data_dir/$job_id"."_select_seqs.faa";
	
##########	
# WARNING: egrep doesn't work with pipe characters (e.g. genbank ids) (selects them all)
# so sequences are selected with a foreach loop instead
	
	my $blastp_filename = "$data_dir/$job_id".".blastp";
	my $blastp_txt = &get_file_text($blastp_filename);
	my @blastp_lines = split "\n",$blastp_txt;
	my $selected_blastp_lines = "";
		
	foreach my $line (@blastp_lines)
	{
		my @tmp = split "\t", $line;
		my $orig_seq = $tmp[0];
		my $qstart = $tmp[6];
		my $qend = $tmp[7];
		my $domain_name = join "_", ($orig_seq,$qstart,$qend);
		if ($num_input_seqs > $max_domain_hits)
		{
			$domain_name =$orig_seq;
			#print "line 184 domain_name = $domain_name";
		}
			
		if (exists $id_hash{$domain_name})
		{		
			my $coords_line = join "\t", ($orig_seq ,$qstart,$qend, $domain_name);
			$selected_blastp_lines .= $coords_line; 
			#print "LINE 191 found selected domain_name=$domain_name <br />";			
		}
		else
		{
			#print "line 195 did not find selected domain_name=$domain_name <br />";
		}
	} 
	$selected_blastp_lines .= "\n";
	&text_to_file($selected_blastp_lines, $trim_aa_coords);	
	
	my $trimmed_fasta_txt =`$pks_script_dir/get_subsequence.pl $trim_aa_coords $allseq_faa_filename`;			
	&text_to_file($trimmed_fasta_txt, $trim_aa_seq_filename);
	
# Trimmed nucleotide sequences already have domain names
#75541525_trimmed_nucleotide_KSdomain_cands.fna
#75541525_trimmed_blast_KSdomain_cands.faa

###############
# But big PCR files <max_sequences use ORFNAMES - might need to fix this here
	#my $allseq_fna_filename = "$data_dir/$job_id"."_orig.fna"; 
###############
	
	my $trim_nuc_seq_filename = "$data_dir/$job_id"."_subset_select_seqs.fasta";
	my $genome_loc_tabfile = "$data_dir/$job_id"."_genomic_locations.tab";
	my $genomic_search_filename = "$data_dir/$job_id" ."_genomic_locations.html"; 
	my $full_search_tabfile = "$data_dir/$job_id"."db_search_results.tab";
	
# for genome search
	#	Query id	Database match	percent identity	align length	e-value	BGC match	domain class	domain subclass	
# KS18	gi|14794889|gb|AF357202.1|_3_24_7401_7826	amphotericin_KS08_cisAT	100	426	9.6e-248	amphotericin	type I modular cis-AT	no subclass


# for aa search
#Query id	Database match	percent identity	align length	e-value	BGC match	domain class	domain subclass	
 #	gi|307660187|gb|CM001016.1|_1_222	furanone_KS04_cisAT	100	222	4.1e-121	furanone	type I modular cis-AT	no subclass
 	
	my $subset_db_search_tabfile = "$data_dir/$job_id"."_subset_search_results.tab";
	my $selected_tab_lines = `head -1 $full_search_tabfile`;

# can't use egrep on genbank names - use foreach loop instead		
	my $allcand_tab_lines = &get_file_text($full_search_tabfile);
	my @alltab = split "\n", $allcand_tab_lines;
	my ($orig_seq,$domain_name);
	foreach my $next_line (@alltab)
	{
		my @tmp = split "\t", $next_line;
	
	# parsing is different for aa versus genomic tabfiles
	#because genomic tabfiles have extra first field with domain number 
	# also, for high volume genomic, $domain_id = $orig_seq_id 
	if ($query_type = "genomic")
	{
		$domain_name = $tmp[1];
	}
	elsif ($query_type = "aa")
	{
		$domain_name = $tmp[0];
	}
		
		if (exists $id_hash{$domain_name})
		{		
			$selected_tab_lines .= "$next_line\n";
			#print "line 279 found selected domain_name=$domain_name <br />";			
		}
		else
		{
			#print "line 283 did not find selected domain_name=$domain_name <br />";
		}	
	}	
		&text_to_file($selected_tab_lines, $subset_db_search_tabfile);	
	
# pass search tab results to html display function to get revised html version of search results page
	$display .= &display_html($selected_tab_lines);	                                
	
# Finish HTML response (custom ending)
	$display .= &get_file_text($custom_page_footer);
	print $display;
	
# clean up
	close SITE_LOG;
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
	
	print SITE_LOG "\n\tERROR: $err_msg\n";
	
	print $cgi->end_html;
	close SITE_LOG;
	#unlink "$data_dir/$pathway_result_page";
	exit (0);
}

sub check_programs
{	
	my @program_list = (
		#"$diamond_dir/diamond",	
		#"$transeq_dir/transeq",
		#"$pks_script_dir/transeq_to_multifasta.pl",
		"$pks_script_dir/getseq_multiple.pl",
		"$pks_script_dir/get_subsequence.pl",
		"$pks_script_dir/select_tab_quant_greater_than_v2.pl",
		"$pks_script_dir/select_tab_instances.pl",
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
	 my $custom_text = join "\n", @text_lines;
	 return $custom_text;		 	
}

sub text_to_file
{
	my ($text, $outfilename) = @_;
	open(TEXTOUT, ">$outfilename") or die "can't open temp file $outfilename for writing.\n $!";
	
	print TEXTOUT "$text";
	close TEXTOUT;
}

#sub get_selected_lines
#{
# 	my ($input_filename, $select_id_nums_filename) = @_;
# 	my $selected_lines = `egrep -f $select_id_nums_filename $input_filename`;
# 	$display .= "$selected_lines\n";

	#my ($input_filename, $id_hash_ref) = @_;
		#my %select_ids = %{$id_hash_ref};	
	# open(INPUT, "<$input_filename") or die "can't open file $input_filename for reading.\n $!";
# 	
# 	while(<INPUT>)
# 	{
# 		chomp;
# 		#next if ($_ =~ /domain\sid/);
# 		if ($_ =~ /domain\sid/)
# 		{
# 			$selected_lines .= "$_\n";
# 		}
# 		#$display .= "$_<br />";
# 		my @tmp = split " ", $_;
# 		my $current_id = $tmp[1];
# 		#$display .= "next id = $current_id<br />";
# 		if (exists $select_ids{$current_id})
# 		{
# 			$selected_lines .= "$_\n";
# 			$display .= "next id = $current_id<br />";
# 		}
# 	}
# 	
# 	close INPUT;
#	return $selected_lines;	
#}


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

sub start_html_file
{

# start html, general page  header
my $text = "";
# qq(<?xml version="1.0" encoding="iso-8859-1"?>
# <!DOCTYPE html
# 	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
# 	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
# <html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US"><head><title>Search Results</title>
# <link rel="stylesheet"  href="../css/main.css" type="text/css" media="screen"  />
# </head><body>);
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

sub display_html
{
	my ($output) = @_;
	my @outlines = split /\n/, $output;
	# first, blank name is for checkbox)
	my @col_names = ("","cand_id","database match","percent identity", 
	"align length","e-value", "BGC product match", "domain class", "domain subclass",);
	
	my ($queryname, $matchname, $pct_id, $align_len, $evalue, $pathway_product, 
	$pathway_product_link, $class, $domain_class_ref, $subclass);
	

	my $v1_webpath ="https://npdomainseeker.sdsc.edu";
 
# javascript addition
 my $html_text .= qq(

<script type="text/javascript">
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
			</script>
			
			);
			
 $html_text .= qq(<h2>Database Search Results</h2>
	Results below are for $num_matched $domain_type domain sequences, from $num_groupnames different categories.<br />);
	my $genomic_locations_html = "$job_id"."_genomic_locations.html";
	
	# if (-s "$genome_loc_tabfile") 
# 	{
# 		$html_text .=	qq(<p>WARNING: Numbers of borderline-quality domain candidates 
# 		identified by nucleic acid searches may <br />not always agree with numbers determined 
# 		by translated amino acid searches, due to differences <br />in blast algorithm alignment scores, 
# 		gaps due to frame shifts and/or stop codons.</p>);	
# 	}
	
	$html_text .= qq(<p>Use check boxes to select candidates for <b><a href="#further_analysis">further analysis</a></b> below. 
	(options may take a few seconds <br />to load for large match tables).</p>);
	

if ($num_matched > 100)
	{
			$html_text .= qq(<br />Note: further analysis options will not appear until table is fully loaded 
			(may take a while for large tables). 
			<br />Maximum number of matches for alignment and tree building is 100.);
	}
	else 
	{
			$html_text .= qq(Click on column headers to sort (multiple clicks toggle between ascending and descending order).);
	}

		
$html_text .= qq(<form  name="align_fasta" action="$v1_webpath/cgi-bin/process_align_request_v2019.cgi" method="post">		
	<input type = "hidden" name = "job_id" value = "$job_id" />
	<input type = "hidden" name = "num_matched" value = "$num_matched" />
	<input type = "hidden" name = "script_source" value = "subseq" />);
	
	$html_text .= qq(<p><input class="button" type="button" value="Select All" onclick="this.value=check(this.form.list)" /></p>);	
	
$html_text .= qq(<!-- Start table of Individual sequence matches -->);	
if ($num_matched <500)
{
	$html_text .= qq(<script src="../../scripts/sorttable.js" type="text/javascript"></script>
	<table class = "sortable2" width = "88%" cellpadding = "5" cellspacing = "0" border = "1">
	<tr>);
}
else 
{
	$html_text .= qq(	
	<table class = "sortable2" width = "100%" cellpadding = "5" cellspacing = "0" border = "1">
	<tr>);
}
	 	
# print headers, but not all will be sortable
	foreach my $colname (@col_names)
		{
			if (length $colname <1)
			{
				$html_text .= qq(<th class="unsortable">$colname</th>\n);
			}
			else
			{
				$html_text .= qq(<th>$colname</th>\n);
			}
		}
	$html_text .= qq(
		</tr>);
		
	$html_text .= qq(</tr><tr><style>
			td {
					word-break: break-all;
			}
			</style>);			
		
# go through each line
		foreach my $line (@outlines)
		{			
		next if ($line =~ /database/i);
		next if ($line =~ /Query/);
		# get rid of extra space at beginning of line
		if ($line =~ /^\s+(.+)/)
		{
			$line = $1;
		}
		$html_text .= qq(<tr>);		
			
			my @cols = split /\t/, $line;
			my $nucleotide_ks_domain = "";
			#if ($compare eq "DNA") # can't use this $compare type parameter doesn't seem to have made it here
		
		# only do this shift on genome sequences where the first column has KS id number instead of queryname
		# this will be either amino acids or very larger queries, where no $genome_loc_tabfile gets written
			if (-s $genomic_search_filename && -s $genome_loc_tabfile)
			{	
				#print "col1=$cols[0]<br />";
				$nucleotide_ks_domain = shift @cols;
			}
			
			
			my $queryname = "$cols[0]";
			my $matchname = "$cols[1]";			
			my $pct_id = sprintf ("%.0f", "$cols[2]");			
			my $align_len = "$cols[3]";
			my $evalue = "$cols[4]";
			my $pathway_product = "$cols[5]" || "unknown"; 
				if ($pathway_product =~ /unknown.+/ || $pathway_product eq "NULL" )
				{
					$pathway_product = "unknown";
				}
				my $pathway_product_link = ("$pathway_product"."_info.html");
				unless (-s "$unixpath_pathway_template_dir/$pathway_product_link")
				{
					$pathway_product_link = "no_path_info".".html";
				}		
			my $class = "$cols[6]";
			my $domain_class_ref = "$class"; # get rid of spaces so can specify a unix file name
				$domain_class_ref =~ s/ /_/g;
				my $subclass = "$cols[7]";
								
	# populate table lne		
			$html_text .=qq(<td style = "text-align: center"><input type = "checkbox" name = "list" value = "$line"></input></td>\n);			
			$html_text .=qq(<td style = "text-align: left">$queryname</td>\n);			
			$html_text .=qq(<td style = "text-align: left">$matchname</td>\n);
			$html_text .=qq(<td style = "text-align: center">$pct_id</td>\n);
			$html_text .=qq(<td style = "text-align: center">$align_len</td>\n);
			$html_text .=qq(<td style = "text-align: center">$evalue</td>\n);
	 		$html_text .=qq(<td style = "text-align: left"><a href =
	 			"$webpath_pathway_template_dir/$pathway_product_link" target="_blank">$pathway_product</a></td>\n);
 			$html_text .=qq( <td style = "text-align: center"><a href="https://npdomainseeker.sdsc.edu/napdos2/background_v2.html#$domain_class_ref" 
 				target="_blank">$class</a></td>\n);
 			$html_text .=qq(<td style = "text-align: center">$subclass</td>\n);
 			$html_text .=qq(</tr>);	
		}

# print rest of html form		
		$html_text .=qq(	
		</table>
		<!-- End Table Individual sequence matches -->	
	<p>	</p>
	<br />\n);
	

# offer opportunity to download tab-delimited file
	my $download_path = "$webpath/tmp2/$job_id/$job_id"."_subset_search_results.tab";
	$html_text .= qq(
	<b>Right-click to <a href="$download_path" target="_blank">DOWNLOAD</a></b> this table in tab-delimited format. 
	<br />);
	
	
$html_text .=qq(

	<a name = "further_analysis"></a>	
			
	<h2>Options</h2>
	
<!-- Start Table 2 analysis options -->		
	<table cellpadding = "5" cellspacing = "5" border = "0">);

	my $genomic_search_filename = "$data_dir/$job_id" ."_genomic_locations.html"; 	
	if (-s $genomic_search_filename && -s $genome_loc_tabfile)
	{
		#print "Found DNA line 1104<br />";
		$html_text .=qq(<tr>
		<td>&nbsp;</td>
		<td>
			<input type = "radio" name = "result_type" value = "genome_location" checked ="checked" ></input>
		</td>				
		<td>
			<b>View nucleotide coordinates </b>of all trimmed domain candidates.<br />
		</td>
	</tr>);	
	}	
	
	
	$html_text .=qq(
	<tr>
		<td>&nbsp;</td>
		<td>
			<input type = "radio" name = "result_type" value = "trimmed_nucs" checked ="checked" ></input>
		</td>				
		<td>
			<b>Output selected sequences </b>in fasta format<br />
		</td>
	</tr>);
	
# don't try to align if number is too high
my $search_results_page = "$webpath/tmp2/$job_id/$job_id"."subclass_tally.html"; 
if ($num_matched > 100)
{
	$html_text .=qq(Return to <a href = "$search_results_page" >previous page </a> 
	and choose a smaller subset (< 100 matches) 
	to access alignment and tree building functions.<br />	  
	<tr>
	<td colspan = "3"> 
	</td>
	</tr>
	</table>
	<!-- End Table 2 analysis options -->
	<table width = "90%"> 		
 	<tr>
 		<td>
 			<p>
				<button type="submit" >GET RESULTS</button>
			</p>			
		</td>
 	</tr>
	</form>		
	</table>);
	return $html_text;
}	
	
$html_text .=qq(	<tr>
		<td>&nbsp;</td>
		<td>
			<p><input type = "radio" name = "result_type" value = "msf"></input></p><br />
		</td> 
			<td><br /><b>Output Alignment</b> with closest database matches<br /> 
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


 $html_text .= qq(<table width = "90%"> 		
 	<tr>
 		<td>
 			<p>
				<button type="submit" >GET RESULTS</button>
			</p>			
		</td>
 	</tr>

	</table>);
		
# Assign pre-align group as all available KS or all available C domaiåns
if ($domain_type eq "KS")
{
	$html_text .= qq(
 <input type = "hidden" name = "align_group" value = "all_classified_KS" />
  );
 }
elsif ($domain_type eq "C")
{
	$html_text .= qq(
	<input type = "hidden" name = "align_group" value = "all_classified_C" </input>
	);
}
else
{
	&send_err("Couldn't find domain_type in domain_blast2.cgi");
}
  
  
$html_text .= qq( </form>		
	);
	return $html_text;
}
