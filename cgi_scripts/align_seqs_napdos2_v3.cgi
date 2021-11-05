#!/usr/bin/perl
# align_seqs_napdos2.cgi
# Sheila Podell
# September 27, 2020

#	Takes in CGI input from process_align_request_v2019.cgi
	# job_id, domain_type (KS or C)s
	# query_type (genome, aa)
	# script_source (allseq, subseq)
	# result_type (trim, align, tree), 
		# trim_type options (nucleic, amino)
		# align_format options (msf, fasta or clw)	
		# tree_type options (closest or align group)
# outputs results to HTML 
	
use warnings; 
use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use Benchmark;	
	
# global variables
	my $num_tree_matches = 5;
	my $max_domain_hits = 200;  # prevent server overload for blastx searches (domain coordinates)
	
# server-specific directories
	my $webpath = "https://npdomainseeker.sdsc.edu";	
	my $unixpath = "/misc/www/projects/npdomainseeker/napdos2/";	
	my $config_filename = "$unixpath/cgi-bin/napdos_dev_sdsc.cfg";
	my $muscle_app = "$unixpath/bin/muscle3.6/muscle";	
	my $tree_app = "$unixpath/bin/FastTree-2.1.11/FastTree";
	my $display_app = "$unixpath/bin/newick-utils-1.6/src/nw_display";
	my $pks_script_dir = "/misc/www/projects/npdomainseeker/scripts/";
	my $subseq_app = "$pks_script_dir/get_subsequence.pl";
	my $dereplicate_app = "$pks_script_dir/remove_replicates.pl";
	my $select_tab_instances_app = "$pks_script_dir/select_tab_instances.pl";
	
# initialize CGI 	
	my $cgi = CGI->new;
	my $source_file = "database_filename";
	my @input_blast_lines_list =  $cgi->param('list'); # contains blast match lines, with both query & match info
	my $queryfile = ""; #$cgi->param('queryfile'); #obsolete: originally came from hmm matches
	my $blastout = $cgi->param('blast_output'); #full, original blast to get closest db matches
	my $query_type = $cgi->param('query_type'); #possible values: genome aa 
	my $result_type = $cgi->param('result_type');	#"trimmed"  "msf" or "tree"
	my $align_format = $cgi->param('align_format');	#"msf"  "fasta" or "clw"
	my $script_source = $cgi->param('script_source'); # "allseq" or "subset"
	my $num_matched = $cgi->param('num_matched'); #probably not important - just what's selected 	
	# debug
		if ($result_type eq 'trimmed_nucs') # where is this coming from?
		{
			$result_type = 'trimmed';
		}
	my $trim_type = $cgi->param('trim_type');	#possble values: "nucleic" or "amino"	
	my $align_group = $cgi->param('align_group'); #not currently implemented
	my $job_id = $cgi->param('job_id');

# relative paths
	my $temp_dir = "/projects/napdos/tmp2";
	my $data_dir = "$temp_dir/$job_id";	
	my $site_log = "$temp_dir/pksdb.log";	
	my $custom_page_header = "$unixpath/napdos_templates/napdos2_header";	
	my $custom_page_footer = "$unixpath/napdos_templates/napdos2_page_footer";
	my $KS_ornament_mapfile ="$unixpath/napdos_templates/KS_refs_ornament.map";	
	my $C_ornament_mapfile ="$unixpath/napdos_templates/C_refs_ornament.map";
	my $pre_aln_file = "";
	if (defined $align_group && !($align_group =~ /blast_hits/))
	{
		$pre_aln_file = "$unixpath/alignments/$align_group"."_2019.afa";

	}	
# Prepare to write HTML response (site-wide templates)
	my $header = &get_file_text($custom_page_header);	
	$header =~ s/\$webpath/$webpath/g;	
	my $footer = &get_file_text($custom_page_footer);	
	
# Start writing HTML response
	print $cgi->header;
	print $cgi->start_html;
	my $html_out = "";	
	
# check to make sure all required accessory programs are available 
	&check_programs;
			
# check to make sure all required accessory programs are available 
	&check_programs;
		
# open logfiles to append messages
	open (SITE_LOG, ">>$site_log") or warn "couldn't open logfile $site_log for appending\n $!\n";
	
# Get parameter info from job_id file, if available, instead of CGI input
	my $param_info_file = "$data_dir/$job_id.params";
	my %saved_params = &read_params_file("$param_info_file");
		$source_file = $saved_params{ref_seq_file};		
		$query_type = $saved_params{query_type}; #genome or aa
		my $domain_type = $saved_params{domain_type}; #KS or C
		my $input_filesize= $saved_params{input_filesize};
		my $num_input_seqs= $saved_params{num_input_seqs};
		my $first_seq_len = $saved_params{first_seq_len} || 0 ;
		my $min_matchlength = $saved_params{min_matchlength};
		my $compare = $saved_params{compare};
		
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
		my $trim_nuc_coords_filename = "$data_dir/$job_id"."_trim_cands_diamond.coords"; #2879427045_trim_cands_diamond.coords # has both orig & orfname
			# $trim_nuc_coords_filename will not exist for large genomic data sets, >2000 input sequences 
		my $selected_nuc_coords_filename = "$data_dir/$job_id"."_selected_trim_nuc_coords"; #will not exist for large data sets, match_num >2000
		my $nuc_blastx_filename	= "$data_dir/$job_id"."_trim_cands_diamond.blastx"; #will not exist for large data sets, match_num >2000
		my $trim_nuc_domains_filename = "$data_dir/$job_id"."_trimmed_nucleotide_"."$domain_type"."domain_cands.fna"; #2879427045_trimmed_nucleotide_KSdomain_cands.fna
		my $selected_trim_nuc_domains_filename = "$data_dir/$job_id"."_selected_nucleotide_domain_cands.fna"; # does not pre-exist	
		my $selected_id_filename = "$data_dir/$job_id"."_selected_ids"; # for large pcr files - no coordinates
		my $selected_domain_id_filename = "$data_dir/$job_id"."_selected_domain_ids";

# get rid of previous versions of selected coords and trimmed files within same job, if exist 
		 if (-e "$selected_id_filename") #$data_dir/$job_id"."_selected_ids; for large pcr files - no coordinates
		{
			unlink $selected_id_filename;
		}
		if (-e "$selected_domain_id_filename") #$data_dir/$job_id"."_selected_ids; for large pcr files - no coordinates
		{
			unlink $selected_domain_id_filename;
		}
		
		 if (-e "$trimmed_aa_domains_fasta_filename") #"$data_dir/$job_id"."_align_select_ids";
 		{
 			unlink $trimmed_aa_domains_fasta_filename;
 		}			
		if (-e "$selected_trim_nuc_domains_filename") #"$data_dir/$job_id"."_align_select_ids";
		{
			unlink $selected_trim_nuc_domains_filename;
		}
				
# get list of currently selected ids, in order
# use hash to ensure uniqueness (duplications kill align and tree steps)
	my %selected_input_ids = (); # key = orf_name  value = number of domains in the orf
	# need orf_ids for short sequences/large PCR data sets where no trimming takes place
	# don't look for multidomain orfs
	my %multi_domain_orfs =(); # key = DOMAIN_name  value = number of domains in the orf
	
	my $input_count = 0; # keeps track of original input order
	foreach my $line (@input_blast_lines_list) # @input_blast_lines_list is a CGI parameter from domain_blast2_diamond.cgi
	{	
		# debug
		#print "line 161 input blast line: $line<br />";
		#gi|14794889|gb|AF357202.1|_3_24 amphotericin_KS08_cisAT 100.0 426 0 0 7401 7826 1 426 9.6e-248 853.6
		#line 161 input blast line: M02610:17:000000000-AACFL:1:2103:3924:14253-1:N:0:30_1_0_2_51 AliivibrioAPE_KS02_arylpolyeneKSa 56 50 3.9e-12 Aliivibrio fischeri aryl polyene type II aryl polyene KSa
		
		$line =~ s/ /\t/gs;
		my @tmp = split "\t", $line;

	# general procedure: for amino acids and genome queries
		my ($orf,$start,$end,$domain_id);
		$input_count++;	
		if ($script_source eq "allseq") 
		{
			$orf = $tmp[0];
			$start = $tmp[6];
			$end = $tmp[7];
			$domain_id = join "_", ($orf,$start,$end);
			#print "line 175 allseq orf=$orf<br /> ";
			#print "line 175 allseq domain=$domain_id<br /> ";
		}
	# selections from grouped subsets already have start and end coordinates appended 
		elsif ($script_source eq "subseq" && $query_type eq "genome")
		{
			$domain_id = $tmp[1];
		}
		elsif ($script_source eq "subseq" && $query_type eq "aa") 
		{
			$domain_id = $tmp[0];
		}
		if ($script_source eq "subseq" && $query_type eq "genome" && ($num_input_seqs > $max_domain_hits)) 
		{
			$domain_id = $tmp[0];
			#print "line 191 subseq large_query domain_id=$domain_id<br />";
		}
			
	# Do I need $selected_input_ids{$orf} the following in case of PCR queries 
	# that are going to be selected directly from orig file, not trimmed?
	# yes, but get them later - first need to deal with extra column (domain id) for input that comes
	# from the "subseq" cgi script
	
		$selected_input_ids{$domain_id} = $input_count;	
		
		if ($query_type eq "genome" && $script_source eq "allseq")
		{
			$selected_input_ids{$domain_id} = $input_count;
			#print "line 186 genome domain=$domain_id <br />";
		}	
		elsif ($script_source eq "subseq" && $query_type eq "genome") 
		{
		 # these lines start with an extra column for KS_#, and 1st column is domain name
			$domain_id = $tmp[1];
		# but this isn't true for genome query for large PCR data sets,  
		# so need to remove the aa coordinates from the orfname	
			if ($script_source eq "subseq" && $num_input_seqs > $max_domain_hits)
 			{		
 				$domain_id = $tmp[0];
 				my @splitname = split "_", $domain_id;
 				pop @splitname;
 				pop @splitname;
 				my $orfname = join "_", @splitname;
 			}	
			$selected_input_ids{$domain_id} = $input_count;			
		}
		elsif ($script_source eq "subseq" && $query_type eq "aa") 
		{
			$selected_input_ids{$domain_id} = $input_count;
			$multi_domain_orfs{$domain_id} = $input_count;	
		}
		
		$multi_domain_orfs{$domain_id} = $input_count;	
	}
	my $num_selected = $input_count;
# check to be sure subsequence program isn't messing up id numbers	
	my @sorted_input_id_list_array = (sort {$selected_input_ids{$a} <=> $selected_input_ids{$b}} keys %selected_input_ids);
	my @sorted_domain_id_list_array = (sort {$multi_domain_orfs{$a} <=> $multi_domain_orfs{$b}} keys %multi_domain_orfs);		
	
# both lists are now the same: contain both trimmed (domain) and untrimmed (orf) sequence id numbers	
	my $selected_input_id_list= join "\n", @sorted_domain_id_list_array;
	my $selected_input_domain_list	= join "\n", @sorted_domain_id_list_array;
			
	&text_to_file($selected_input_id_list,$selected_id_filename);
	&text_to_file($selected_input_domain_list,$selected_domain_id_filename);

#############################
# Trimming starts here
#############################
# select desired subsets of pre-trimmed aa and nucleotide sequences
	my ($trimmed_aa_seq_txt, $trimmed_nt_seq_txt, $unique_trimmed_nt_seq_txt);
	my @selected_domains = keys %selected_input_ids;
	
	my $num_selected_nt_coords = 0;
	my $selected_nucleotide_ids = "";
	if ($query_type eq "genome")
	{
	# this subroutine trims AND writes aa results to file
		$trimmed_aa_seq_txt = &trim_aa_domains(\%selected_input_ids, $aa_blastp_tophits_tab_filename,
		$translated_orf_filename, $trimmed_aa_domains_fasta_filename);
				
	# nucleotides have two cases:
		# case 1: need to select trimmed sequences from all domain candidates using nucleotide coordinates
			# domain candidates are named by protein coordinates, not nucleotide coordinates
			# so need to use selected blast lines to get them
			# can't just use orf id numbers, because there may be >1 domain per orf 
		 # case 2: $trim_nuc_coords_filename file does not exist (because match number is very high)
		 	# can select nucleotides using just id numbers, because nothing has been trimmed
	
	#previously generated by genome_search2_diamond.cgi
		 my $genomic_locations_tabfile = "$data_dir/$job_id"."_genomic_locations.tab";
		 my $parent_match_seqs_fasta = "$data_dir/$job_id"."_parent_match_seqs.fna";	 
		 
	# case 1: retrieving the selected nucleotide domains from a previously trimmed 
	# nucleotide file using domain ids (created by genome_search2_diamond.cgi)
	# this works - don't change it
		if ($num_input_seqs < $max_domain_hits) # || $min_matchlength > 150)
		{			
			$trimmed_nt_seq_txt = `$pks_script_dir/getseq_multiple.pl $selected_id_filename $trim_nuc_domains_filename`;
			$unique_trimmed_nt_seq_txt = &remove_duplicate_seqs($trimmed_nt_seq_txt);
				
		# text_to_file subroutine requred because CGI doesn't allow direct creation of files from unix command line
			&text_to_file($unique_trimmed_nt_seq_txt, $selected_trim_nuc_domains_filename);				 		
		}	
	# if >2000 nucleotide domain matches originally found, just select untrimmed fastas from parent_seqs
		elsif ($num_input_seqs > $max_domain_hits) #  || ($min_matchlength < 150)) 
		{	
			# 	my $selected_id_filename = "$data_dir/$job_id"."_align_select_ids";
				# M02610:17:000000000-AACFL:1:2108:17857:20192-1:N:0:30_1_0_1_85
				# these ids have frame, orf number, and start/stop appended
				# need to pop these off before retrieving untrimmed sequences from original input
	 		$selected_nucleotide_ids = "";
			foreach (@sorted_input_id_list_array)
			{
			chomp;
			 	my @tmp = split "_", $_;
			 	my $aa_end = pop @tmp;
			 	my $aa_start = pop @tmp;
			 	my $frame_num = pop @tmp;
			 	my $orf_num = pop @tmp;
			 	my $parent_id = join "_", @tmp;
			 $selected_nucleotide_ids .= "$parent_id\n";			 
		}
		&text_to_file($selected_nucleotide_ids, $selected_id_filename);
			if ($script_source eq "subseq")
			{
				$orig_nuc_filename = "$data_dir/$job_id".".fasta";
				# to find parent ids without frame info"
			}	
			$trimmed_nt_seq_txt = `$pks_script_dir/getseq_multiple.pl $selected_id_filename $orig_nuc_filename`;
			&text_to_file($trimmed_nt_seq_txt,$selected_trim_nuc_domains_filename);
			$unique_trimmed_nt_seq_txt = &remove_duplicate_seqs($trimmed_nt_seq_txt); 	
			&text_to_file("$unique_trimmed_nt_seq_txt", "$selected_trim_nuc_domains_filename");
		}
	}
	# for aa query, don't need to deal with nucleotide sequences
	elsif ($query_type eq "aa")
	{
		$trimmed_aa_seq_txt = &trim_aa_domains(\%selected_input_ids, $aa_blastp_tophits_tab_filename,
		$orig_aa_fasta_filename, $trimmed_aa_domains_fasta_filename);			
	}
	else 
	{
		&send_err ("Couldn't determine sequence input type for $query_type");
	}	

	my $num_selected_input_ids = scalar (keys %selected_input_ids);	
		
# if result_type = trimmed, write outstring to html <pre></pre>, end script
	if ($result_type eq "trimmed")
	{
		if ($trim_type eq "nucleic")
		{
			print "<pre>$unique_trimmed_nt_seq_txt<pre>";
		}
		elsif ($trim_type eq "amino")
		{
			print "<pre>$trimmed_aa_seq_txt<pre>";
			unless (length $trimmed_aa_seq_txt > 10)
			{
				print "<pre>Sorry, no amino acid sequences available for selection. Try nucleic acid instead. <pre>";
			}
		}
		else
		{
			&send_err ("print trim_type=$trim_type :unrecognized trim result type");
		}
	}
	
# check all sequence id numbers for illegal characters
# previous version looked only at first sequence - problem might be later
	if ($result_type eq "msf") #|| $result_type eq "tree")
	{		
		my $queryfile = "$trimmed_aa_domains_fasta_filename";
		my $query_id_cmd ="grep '>' $queryfile";
 		my $test_out = `$query_id_cmd`;
 		$test_out =~ s/>/\n/g;
 		my @ids = split "\n", "$test_out";
 		foreach my $test_id (@ids)
 		{
 			next unless (defined $test_id && length $test_id > 2);
 			my @tmp = split " ", $test_id;
			
 			if ($tmp[0] =~ /([:,])/)
			{
				&send_err("Query id number contains illegal character(s) ($1) <br /> $tmp[0]");
			}
			if (length $tmp[0] > 79)
			{
				&send_err("At least one query id number exceeds maximum length (79 chars) <br /> $tmp[0]");
			}		
 		}
	}	
	
# check to be sure < 100 matches are selected and result_type is align or tree,
	# else, end html and exit
 # pass selected trimmed aa file (or orf file) to aligment function (run muscle)
	# don't do alignments or trees on big files(might want to send error message, for safety)

#############################
# Alignment starts here
# don't do any of these steps if ($num_input_seqs < $max_domain_hits)
#############################

# run muscle for msf, output to file, keep track of time
	my $t_muscle_start = new Benchmark;
	my $aln_filename = "$data_dir/$job_id"."_align_out.fasta";
	my $tree_aln_filename = "$data_dir/$job_id"."_tree_aln.afa";
	my $align_out = "";
	my $muscle_output_type = "fasta";
	my $extended_db_fastas = "";		
		
# retrieve extended set (top $num_tree_matches) database matches (default 5)
	my $rawblast_file = "$data_dir/$job_id".".blastp".".raw";
	unless (-s $rawblast_file)
	{
		&send_err ("Couldn't find blastp file $rawblast_file<br />");
	}
	
# retrieve unique (tophit) database match fastas into string, with newlines instead of html breaks
	# can't use alignment file because it only contains a subset of potential blast matches
	my $source_file_path = "$unixpath/fasta/$source_file";	# fasta sequences of matches ($saved_params{ref_seq_file};	)
			#ref_seq_file=all_KS_191020_1877.faa	
	unless (-s "$source_file_path")
	{
		&send_err ("Couldn't find matched ref seq fasta file $source_file_path");
	}

# extract database matches to the currently selected sequences only from the master alignment file (.afa)	
	my $select_instances_filename = "$job_id"."_selected_matches.txt"; # a list of 2 columns: id number and closest match
	my $select_instances_txt = `$select_tab_instances_app $rawblast_file $num_tree_matches |cut -f1,2 >$data_dir/$select_instances_filename` ;
	
	
	#my $debug = "$select_tab_instances_app $rawblast_file $num_tree_matches |cut -f1,2 >$data_dir/$select_instances_filename" ;
	#print "$debug<br />";
		
	my %extended_db_id_matchlist = ();	# key = matched db sequence id # value = user query protein id
	# key = db match because don't want to duplicate matches
	
	 open (SELECT_DB_MATCHES, "$data_dir/$select_instances_filename") or warn "can't open select db matches file $select_instances_filename\n<br />$!\n";
 	 while (<SELECT_DB_MATCHES>)
 	 {
 	 	chomp;
 	 	my @tmp = split "\t", $_;
 		my $query = $tmp[0];
 		my $db_match = $tmp[1];
 		my $qstart = $tmp[6];
 		my $qend = $tmp[7];
 		my $domain_id = join "_", ($query,$qstart,$qend);
 		#$extended_db_id_matchlist{$db_match}=$query;
 		$extended_db_id_matchlist{$db_match}=$domain_id;
	 }
	 
	 close SELECT_DB_MATCHES;

	my @extended_db_matchlist = keys %extended_db_id_matchlist;
	my $count_extended_matchlist = scalar @extended_db_matchlist;
	#print "found $count_extended_matchlist extended matchlist, line 356<br />";	
	
	# source file is reference database should have domain ids, not orf ids
	$extended_db_fastas = &get_fastas(\@extended_db_matchlist, $source_file_path);
	
	# add $extended_db_fastas to $trimmed_aa_domains_fasta_filename and remove duplicates
		my $extended_fastas = "$extended_db_fastas"."\n"."$trimmed_aa_seq_txt";	
		my $unique_extended_fastas = &remove_duplicate_seqs($extended_db_fastas);
	# debug
		#print "$extended_db_fastas<br />";
		my $count_extended_fastas = ($extended_fastas =~ tr/>/>/);
		my $count_unique_extended_fastas = ($unique_extended_fastas =~ tr/>/>/);
		#print "found $count_extended_fastas extended fastas, $count_unique_extended_fastas unique.<br />";

# run muscle on extended protein matchfile 
		#$align_out = &run_muscle($unique_extended_fastas, $trimmed_aa_seq_txt, $aln_filename, $align_format);
		#sub run_muscle ($db_seqs, $query_seqs, $aln_filename, $muscle_out) = @_;
		
	# need to do extra alignments, if user wants msf or clw format, because these can't be used for trees
		my $aln_filename_msf = "$data_dir/$job_id"."_align_out.msf";
		my $aln_filename_clw = "$data_dir/$job_id"."_align_out.clw";
		
	# get rid of old versions of alignment files, in case users have chosen a different sequence subset
		if (-e $aln_filename_msf)
		{
			unlink $aln_filename_msf;
		}
		if (-e $aln_filename_clw)
		{
			unlink $aln_filename_clw;
		}
		if (-e $aln_filename)
		{
			unlink $aln_filename;
		}
	
# run muscle with user-selected format			
	if ($num_selected < $max_domain_hits)
	{ 	
		if ($align_format eq "msf")	
		{			
			$align_out = &run_muscle($unique_extended_fastas, $trimmed_aa_seq_txt, $aln_filename_msf, $align_format);			
		}
		elsif ($align_format eq "clw")
		{
			$align_out = &run_muscle($unique_extended_fastas, $trimmed_aa_seq_txt, $aln_filename_clw, $align_format);							
		}
		elsif ($align_format eq "fasta")
		{
			$align_out = &run_muscle($unique_extended_fastas, $trimmed_aa_seq_txt, $aln_filename, $align_format);
		}
	}
	else
	{
		#print "too many input sequences to run muscle: num_input_seqs=$num_input_seqs <br />";
		
		#print "too many input sequences to run muscle: num_selected_seqs=$num_selected <br />";
		unless ($result_type eq "trimmed") 
		{
		print "too many input sequences to run muscle: num_selected_seqs=$num_selected <br />";			
		}				
	}	
		
	my $t_muscle_finish = new Benchmark;
	my $muscle_time = timediff ($t_muscle_finish, $t_muscle_start);
	print SITE_LOG "\tMUSCLE ALIGNMENT:", timestr($muscle_time, 'nop'), "\n";
	
# get web paths for users to download alignment files (html users can't access unix paths)
	if ($result_type eq "msf" && $num_selected < $max_domain_hits)
	{	
		my $web_access_aln_filename = "$webpath/tmp2/$job_id/$job_id".".$align_format";
		if ($align_format eq "fasta") 
		{
			$web_access_aln_filename = "$webpath/tmp2/$job_id/$job_id"."_align_out.fasta";
		}
		elsif ($align_format eq "msf") 
		{
			$web_access_aln_filename = "$webpath/tmp2/$job_id/$job_id"."_align_out.msf";
		}
		elsif ($align_format eq "clw") # $job_id".fasta gets original fasta input, not selected alignment
		{
			$web_access_aln_filename = "$webpath/tmp2/$job_id/$job_id"."_align_out.fasta";
		}
		$html_out = &display_msf($web_access_aln_filename, $align_out);
	
	}
	# $html_out is written to html screen at end of this script

#############################
# Tree starts here
#############################

# run FastTree
	my $t_tree_start = new Benchmark;
	my $tree_filename = "$data_dir/$job_id.tre";
	my $tree_graphic_prettyfile = "$data_dir/$job_id.tre.svg";
	my $tree_graph_ascii = "$data_dir/$job_id.tre.ascii";
	my $newick = "";
	my @select_names = (keys %selected_input_ids);
	# need to make sure get ornaments for multiple independent domains within a single ORF
	# fasta amino acid sequences were re-named with coordinates
		
	open (TRIM_SEQS, $trimmed_aa_domains_fasta_filename) or die "unable to open trimmed aa file $trimmed_aa_domains_fasta_filename.<br />$!<br />";
	my $num_trim_seqs = 0;
	while (<TRIM_SEQS>)
	{
		if ($_ =~ /^>(.+)\s/)
		{
			push @select_names, $1;
			$num_trim_seqs++;
		}
	}	
	close TRIM_SEQS;
					
	#create a nw_display SVG mapfile for user selected sequences (give them a red dot)
	# make sure the marked names include both original ids & ids with coordinates appended			
		#my @select_names = (keys %selected_queries, keys %selected_queries_uniq, @modified_names);
		my $select_ornament_mapfile = "$data_dir/$job_id"."_ornaments.map";
		open (ORNAMENTS, ">$select_ornament_mapfile") or warn "couldn't open mapfile for selected seqs $!\n";
		 foreach my $name (@select_names)
 		{
 			print ORNAMENTS qq("<circle style='fill:red;stroke:black' r='5'/>" I $name \n);
 		}		
		close ORNAMENTS;
	
# build one of two tree types, based on closest neighbors or group reference	
	if ($result_type eq "newick" || $result_type eq "tree")
	{				
	# make sure there is an alignment in the correct format for FastTree	
		unless (-s "$aln_filename")
		{
			#print "Couldn't find or no data in alignment input file for trees, $aln_filename";
			$align_format = "fasta";
			$align_out = &run_muscle($unique_extended_fastas, $trimmed_aa_seq_txt, $aln_filename, $align_format);									
		}
		
		if ($align_group eq "closest")
		{
			$newick = &runFastTree($aln_filename, $tree_filename);
		}
		elsif ($align_group eq "all_classified_KS" || $align_group eq "all_classified_C")
		{
		# need to re-run muscle to combine $aln_filename (blast hits) with reference file	
			my $combo_afa_filename = "$data_dir/$job_id"."_ref_plus_queries.afa";
			my $dereplicted_combo_afa_filename = "$data_dir/$job_id"."_ref_plus_queries_uniq.afa";
	
			my $profile_cmd = "$muscle_app -profile -in1 $pre_aln_file -in2 $aln_filename -out $combo_afa_filename";
			`$profile_cmd`;
			#print "<br />$profile_cmd<br />";
			
			# need to remove duplicates for FastTree to work ($aln_filename  contains orig + reference matches )
			my $dereplicate_cmd = "$dereplicate_app $combo_afa_filename >$dereplicted_combo_afa_filename";
			#print "<br />$dereplicate_cmd<br />";
			`$dereplicate_cmd`;
					
			$aln_filename = "$dereplicted_combo_afa_filename";
	
			#print "align group = $align_group <br />";
			#print "align file = $aln_filename <br />";
			#print "pre-align file = $pre_aln_file <br />";
			
			$newick = &runFastTree($aln_filename, $tree_filename);
					
		}
		else
		{
			&send_err ("could't find align group name $align_group line 393 align_seqs2_v2.cgi");
		}
			
		unless (-s $tree_filename)
		{
			&send_err("can't find treefile line 399");
		}
		my $ugly = &create_graphic($tree_filename, $tree_graph_ascii);
		my $pretty = &create_graphic_svg($tree_filename, $tree_graphic_prettyfile,$select_ornament_mapfile );
				
		$html_out = &display_tree($tree_filename, $newick, $ugly, $tree_graphic_prettyfile);
	}
	my $t_tree_finish = new Benchmark;
	my $tree_time = timediff ($t_tree_finish, $t_tree_start);
	if ($result_type eq "newick")
	{
		print SITE_LOG "\tBUILD TREE:", timestr($tree_time, 'nop'), "\n";
	}
	# debug
		#print qq(<pre>$newick</pre>);
		#my $debug = &display_file($tree_filename);
		#print qq(Tree filename = $tree_filename<br />
		#<pre>$debug</pre>);
	
# send results to html, finish response
	print "$html_out";
	#print $footer;
	
	close SITE_LOG;



############################################
# SUBROUTINES
############################################

sub send_err
{	
	my ($err_msg) = @_;
	print qq(<p>Sorry, your request could not be processed.</p>);
	print $err_msg;
	
	 print JOB_LOG "\nERROR: $err_msg\n";
	
	print $cgi->end_html;
	close SITE_LOG;
	#close JOB_LOG;
	exit (0);
}

sub check_programs
{	
	my @program_list = (
		$muscle_app,
		$tree_app,
		$display_app,
		$dereplicate_app,
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
	
	unless (scalar keys %param_input > 1)
	{
		&send_err("No saved parameters recognized for job $job_id");
	}
    return %param_input;
 }
 
sub trim_aa_domains
{	
	my ($selected_ids_hashref, $blastp_filename, $source_fasta_filename, $destination_fasta_filename) = @_;	
	my %selected_ids_hash = %{$selected_ids_hashref};
	my $trim_coord_filename =  "$data_dir/$job_id"."_aa_match_trim_coords";
	if (-s "$trim_coord_filename")
	{
		unlink "$trim_coord_filename";
	}

# select coord lines for selected id numbers, write to file along with domain name
# DON'T select using egrep- problem with pipe characters in gbk id numbers!
	my $selected_coords = "";
	my $domain_num = 0;
	my %orf_domain_tally = ();  # key = orf_id   # value = # domains per orf
	open (BLASTOUT, "<$blastp_filename") or warn "couldn't open blastp file $blastp_filename for reading\n $!\n";
	while (<BLASTOUT>)
	{
		chomp;
		#print "$_<br />";
		my @tmp = split "\t", $_;
		my $orig_id = $tmp[0];
		my $id = $tmp[0];
		my $start = $tmp[6];
		my $end = $tmp[7];
	
		my @id_tmp  = split "_", $id;
		$id = join "_", @id_tmp;
		my $domain_id = join "_", ($orig_id,$start,$end);
		
	# create reference file that links orig_id (aa id) to domain id
		if (exists $selected_ids_hash{$domain_id})   
		{
			$orf_domain_tally{$id}++;  # don't need this for aa matches
			$domain_num = $orf_domain_tally{$id}; # don't need this for aa matches
			#debug
				#print "line 746 found domain $domain_num for $id<br />";
			my @coords = ($id, $start, $end);
			my $new_name = join "_", @coords; 
			my $coords_line = join "\t", (@coords, $domain_id);
			$selected_coords .= "$coords_line\n";			
		}
		else
		{
			#debug
				#print "line 759 domain $domain_id not found in selected ids hash <br />";
		}
	}
	close (BLASTOUT);
	
	&text_to_file($selected_coords, "$trim_coord_filename");

# get fastas for selected ids from source, write to destination fasta file		
	my $trimmed_fasta_txt =`$pks_script_dir/get_subsequence.pl $trim_coord_filename $source_fasta_filename`;

# make sure there are no orf id duplicates	
	my $unique_trimmed_fasta_txt = &remove_duplicate_seqs("$trimmed_fasta_txt");
	&text_to_file($unique_trimmed_fasta_txt, "$destination_fasta_filename"); 			
	return $trimmed_fasta_txt;
}

sub get_fastas
{
	my ($array_ref, $filename) = @_;
	my @namelist = @$array_ref;
	my $fastas = "";
	my %namehash = ();
	my $result = "";
		
# be sure there are sequences to get
	unless (scalar @namelist >0)
	{
		&send_err("Get fastas subroutine - no sequences selected.");
	}	
	
#  get list of unique id's into hash for sequence retrieval
	foreach (@namelist)
	{
		$namehash{$_} = " ";
	}
	
# retrieve sequences
	my $num_requests = scalar (keys %namehash);	
	my $select = 0;		# 0=false, 1= true
	my $count_hits = 0;
	my $sequence = "";
	open (INFILE, "<$filename") or send_err("Could not open database file $filename<br>");
	my $count = 0;
	while (my $line = <INFILE>)
	{
	# get header, (only) if id # is on selection list 
		if ($line =~ /^\>(.+)/)
		{
			last if ($count_hits == $num_requests);
			my $header = $1;			
			chomp $header;
			my @seq_words = split (" ", $header);
			my $seq_id = $seq_words[0];	
			if (exists $namehash{$seq_id})
			{
				$select = 1;
				$count++;
				$result .= "$line";
				$count_hits++;
			}
			else {$select = 0};
			next;
		}				
	# get body of sequence, (only) if id# was on selection list
		if ($select > 0)
		{
			$result .= "$line";
		}		
	}	
	close INFILE;

	return $result;
}

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
  # if (exists $selected_queries{$self->{id}})
# 	{
# 		$self->{selected} = 1;
# 	}
# 	else
# 	{
# 		$self->{selected} = 0;
# 	}
  return($self)
}


sub display_msf
{
	my ($filename, $display_txt) = @_;
	my $weblink = $filename; 
	my $webpath = "https://npdomainseeker.sdsc.edu";
	$weblink =~ s/$unixpath/$webpath/;	
	#print "subroutine looking for $filename<br />";
# page header
	my $html_out = $header;
	
	 $html_out .= qq(<link rel="stylesheet" href="/css/main.css" type="text/css" ></link>);
	 
	 $html_out .= qq(
	 
		<h2>Alignment results</h2>
	<p></p>);
	

# buttons for file retrieval	
 	$html_out .= qq(
 	Right click to <b><a href="$filename" target="_blank">DOWNLOAD</a></b> file
 	</form><br />);


# alignment file
	$html_out .= qq(<hr /><font face = "monospace"><pre>$display_txt</pre></font>\n<br />);

return $html_out;
}

sub display_tree
{
	my ($filename, $display_txt, $tree_ascii, $tree_pretty_file) = @_;
	#$html_out = &display_tree($tree_filename, $newick, $ugly, $tree_graphic_prettyfile);
	my $weblink = $filename;
	my $webpath = "https://npdomainseeker.sdsc.edu";
	my $user_accessible_path = "$webpath/tmp2/$job_id";
	 
	$weblink =~ s/$data_dir/$user_accessible_path/;
	$tree_pretty_file =~ s/$data_dir/$user_accessible_path/;
	
	my $tree_prettyfile_webpath = $tree_pretty_file;
	$tree_prettyfile_webpath =~ s/$unixpath/$webpath/;
	unless (defined $tree_ascii)
	{
		$tree_ascii = "no ascii data found";
	}
		
# page header
	my $html_out = $header;
	
	 $html_out .= qq(<link rel="stylesheet" href="/css/main.css" type="text/css" ></link>);
	 
	 $html_out .= qq(
	 
		<h2>Tree results</h2>
	<p> Additional biosynthetic gene cluster information for reference database matches in these trees, 
	is provided in the <br />
	<a href="$webpath/napdos2/pathways_v2.html"><b>BGCs</b></a> reference tab.</p>);
	#http://npdomainseeker.sdsc.edu//pathways_overview_pksdb13.html
	
$html_out .= qq(<table cellpadding = "10" cellspacing = "10">
<tr>
	<td>
		<h3>Tree image (svg format file)</h3>		
	</td>
	<td>	
		Right click to <b><u><a href="$tree_prettyfile_webpath" target="_blank" > DOWNLOAD </a></u></b>
	</td>

</tr>
<tr>
	<td>
		<h3>Newick format (plain text file)</h3>		
	</td>
	<td>		
		Right click to <b><u><a href="$weblink" target="_blank" >DOWNLOAD</a></u></b>
	</td>	
</tr>
</table>
	);

# alignment file
	#$html_out .= qq(<hr /><font face = "monospace"><pre>$display_txt</pre></font>\n<br />);

return $html_out;
}

sub pretty_seq
{
	my ($seq) = @_;
	my $pretty_seq = "";
	my $count = 0;
	my $line_length = 60;
	#$seq =~ s/ //sg;
	while ($count < length $seq)
	{
		$pretty_seq .= substr ($seq, $count, $line_length);
		$pretty_seq .= "\n";
		$count += $line_length;
	}
	chomp $pretty_seq;
	return $pretty_seq;
}
sub run_muscle
{
	my ($db_seqs, $query_seqs, $aln_filename, $muscle_out) = @_;
	# print "run_muscle_input_params:<br /> 	
 	#aln_filename=$aln_filename<br />
 	#muscle_out=$muscle_out<br /> <br />";
 	#print "db_seqs=$db_seqs<br />";
 	#print "query_seqs=$query_seqs<br />";
	
# check for muscle program
	unless (-s $muscle_app && -e $muscle_app)
	{
		&send_err("Couldn't find muscle alignment program $muscle_app.");
	}
	
# get blast-match sequences into temp file - this already exists
	#my $match_filename = "$data_dir/$job_id.matches";
	my $match_filename = "$data_dir/$select_instances_filename";
	unless (-s $match_filename)
	{
		print "run_muscle can't find file $match_filename <br />";
	}
	
	#2685528328_selected_matches
	
	#substitute newlines for html line breaks
	$db_seqs =~ s/\<br\>/\n/g;
	
# get user input sequences into temp file
	my $query_fasta_filename = "$data_dir/$job_id.trim";	
	#substitute newlines for html line breaks
	$query_fasta_filename =~ s/\<br\>/\n/g;
	open(OUTFILE, ">$query_fasta_filename") or  send_err qq(<br />can't open temp file $query_fasta_filename for writing.\n $!);
	print OUTFILE "$query_seqs";	
	close OUTFILE;
	
# combine both query and match sequences into a single file
	my $combo_filename = "$data_dir/$job_id.combo.faa";
	open(OUTFILE, ">$combo_filename") or  send_err qq(<br />can't open temp file $combo_filename for writing.\n $!);
	print OUTFILE "$query_seqs";
	print OUTFILE "$db_seqs";
	close OUTFILE;
	
	unless (-s $combo_filename)
	{
		print "run_muscle can't find file $combo_filename <br />";
	}
		
# de-novo alignment	
	my $cmd = "$muscle_app -in $combo_filename -$muscle_out";
	#print "run_muscle command <br />$cmd<br />";
	my $output = `$cmd`;
	
# remove duplicate fastas in alignment - muscle doesn't care, but FastTree chokes on them
	my $revised_output = $output;
	if ($align_format eq "fasta")	
	{
		$revised_output = &remove_duplicate_seqs($output);
		#print "align_seqs2 LINE 544<br />";
	}
	
# send program output to temp file	
	open(ALIGN, ">$aln_filename")
		or &send_err("<br />can't open temp file $aln_filename for writing.\n $!");
	print ALIGN "$revised_output";	
	close ALIGN;
			
# check to be sure 1st muscle output was successful
	unless (-s "$aln_filename")
	{
	
		&send_err("ERROR: run_muscle line 887 can't find Muscle program alignment output file: $aln_filename.");
	}	

# alternative option: used pre-aligned group as fixed "profile", query sequences as new input
	
	if (defined "$pre_aln_file"  && "$muscle_out" eq "fasta" && -s "$pre_aln_file")
	{
		 unless (-s "$pre_aln_file")
 		{
 			&send_err("ERROR: Can't find pre-aligned reference file $pre_aln_file line 886");
 		}
		
		# send program output to temp file	
		# open(ALIGN, ">$aln_filename")
# 			or &send_err("<br />can't open temp file $aln_filename for writing.\n $!");
# 		print ALIGN "$revised_output";	
# 		close ALIGN;
	}
	
# clean up temporary files
	#unlink $match_filename;
	#unlink $query_fasta_filename;
	#unlink $combo_filename;
	
	return $revised_output;
}

sub runFastTree
{
	my ($align_file, $treefile) = @_;
	
# double-check alignment file for duplicate or weird names
	my @namelist = &get_fasta_ids_from_file($align_file);
	my %namecount = ();	#key = id, value = # times appears
	
	foreach my $name (@namelist)
	{
		if ($name =~ /(\|:)/ || $name =~ /([,)(:])/)
		{
			&send_err("Sorry, can't build tree: illegal character(s) in sequence name:<br /> $name");
		}
		
		if (length $name > 61)
		{
			&send_err("Sorry, can't build tree: query id length exceeds maximum for tree software (55 chars) <br /> $name");
		}
		
		if (exists $namecount{$name})
		{
			&send_err("Duplicate name in muscle alignment $name");
		}		
		else
		{
			$namecount{$name}++;
		}
	}
	
	
	my $cmd = "$tree_app $align_file";
	my $tree_output = `$cmd`;
	#print "$cmd<br />";
	
	&create_outfile($tree_output, $treefile);
	unless (-s "$treefile")
	{
		&send_err("ERROR: Can't find Tree program alignment output file $treefile");
	}
	
	return "$tree_output";
}

sub create_graphic
{
	 my ($treefile, $graphic_file) = @_;
	 	
 	my $cmd = "$display_app $treefile";
 	my $graphic_output = `$cmd`;
 	
 	&create_outfile($graphic_output, $graphic_file);
 	unless (-s "$graphic_file")
 	{
 		&send_err("ERROR: Couldn't create graphic file $graphic_file \n $!");
 	}
 	
 	return "$graphic_output";
}
sub create_graphic_svg
{
	 my ($treefile, $graphic_file, $ornament_file) = @_;
	 
	#$ornament_file = "KS_refs_ornament.map";
	
# nw_display -s -S -v 40 -w 500 -b 'opacity:0' -i 'opacity:0' -l 'font-family:serif' \
 #	 -o test_ornament.map CNT088_KS_domains_big.tre > test2.svg	
	 	
 	#my $cmd = "$display_app -s $treefile";
 	
 	my $cmd = qq($display_app  -s -S -v 40 -w 500 -b 'opacity:0' -i 'opacity:0' -l 'font-family:serif' -o $ornament_file $treefile);
 	
 	my $graphic_output = `$cmd`;
 	
 	&create_outfile($graphic_output, $graphic_file);
 	unless (-s "$graphic_file")
 	{
 		&send_err("ERROR: Couldn't create graphic svg file $graphic_file \n $!");
 	}
 	
 	return "$graphic_output";
}

sub display_file
{
	my ($filename) = @_;
	my $output_text = "";
	open (INFILE, "$filename") or send_err("Could not open file for display $filename <br>$!<br>");
	while (<INFILE>)
	{
		$output_text .= "$_";
	}
	close INFILE;
	return $output_text;
}

sub create_outfile
{
	my ($input, $out_filename) = @_;
	open(OUTFILE, ">$out_filename") or  send_err qq(<br />can't open temp file $out_filename for writing.\n $!);
	print OUTFILE "$input";	
	close OUTFILE;
}

sub remove_duplicate_seqs
{
	my ($fasta_text) = @_;
	my @lines = split '\n', $fasta_text;
	my $num_lines = scalar @lines;
	# debug
		#print "found $num_lines lines before de-replication <br />";		

	my %object_list =();	#key = ID number value = sequence object
	my %id_tally = ();
	my @seq_list = ();	# list of id numbers - preserves original order
	my $current = "";

	foreach my $line (@lines)
	{
		next if $line =~ /^\s+$/;
		chomp $line;
	# get rid of multiple spaces	
		$line =~ tr / / /s;
		if ($line =~ />(.+)/) 
		{								
			my @tmp = split " ", $line;
			$current = &new_seq("record",\@tmp);
			$object_list{$current->{id}} = $current;
			$id_tally{$current->{id}}++;
			#print "$current->{id} = $id_tally{$current->{id}}<br />";
			
	# note: this filter may fail to detect duplicates for genbank sequence ID's with pipe characters (perl hash exists bug)	
			unless ($id_tally{$current->{id}} > 1)
			{
				push (@seq_list, $current->{id});
			}			
		}													
		else
		{			
			$current->{sequence} .= uc"$line";			
		}	
	}
	
# write to output string
	my $output = "";
	my $count = scalar @seq_list;
	#print "found $count lines after de-replication <br />";
	foreach my $next (@seq_list)
	{
		my $object = $object_list{$next};
		if (exists $object_list{$next} && defined $object->{sequence})
		{
			
			my $pretty_seq = &pretty_seq($object->{sequence});
			$output .= "$object->{header}\n$pretty_seq\n";
			#print "$object->{header}\n$pretty_seq<br />";
		}
		else
		{
			print "no sequence data for $next<br />"
		}
	}
	
	unless (defined $output && length $output > 0)
	{		
		#&send_err("ERROR removing duplicate sequences\n");
	}
	
	
	return $output;	
}
sub get_fasta_ids_from_file
{
	my ($filename) = @_;
	my @id_nums = ();
	open (FILE, "$filename") or warn "unable to open align file $!";
	while (<FILE>)
	{
		if ($_ =~ /^\>(.+)/)
		{						
			my $header = $1;			
			chomp $header;
			my @seq_words = split (" ", $header);
			my $seq_id = $seq_words[0];	
			push @id_nums, "$seq_id";
		}	
	}
	close FILE;

return @id_nums;
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
	
	print TEXTOUT "$text\n";
	close TEXTOUT;
}
sub get_parent_ids
{
	my ($arrayref) = @_;
	my @orf_id_list = @{$arrayref};
	my @parent_id_list = ();
	
	foreach (@orf_id_list)
	{
		chomp;
		 my @tmp = split "_", $_;
		 pop @tmp;  # remove terminal reading frame
		 pop @tmp; 	# remove terminal orf_id_num
		 my $next_id = join "_" , @tmp;
		push @parent_id_list, $next_id;		 
	}
	return @parent_id_list;
}	


__END__
	
		