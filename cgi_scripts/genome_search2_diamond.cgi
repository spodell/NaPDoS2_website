#!/usr/bin/perl
# genome_search2_diamond.cgi
# Sheila Podell
# January 13, 2020

# find and tally candidate domains from genomes and metagenomes
# revised December 2019 eliminate use of HMMER

use warnings; 
use strict;
use File::Basename;
use File::Spec;
use Cwd;

use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use Benchmark;

# global variables	
	my $num_processors = 2;
	my $max_domain_hits = 2000;  # prevent server overload for blastx searches (domain coordinates)
	
# server-specific directories
	my $webpath = "https://npdomainseeker.sdsc.edu";
	my $unixpath = "/misc/www/projects/npdomainseeker";
	my $diamond_dir = "$unixpath/bin/diamond-2.0.15";	
#	my $diamond_dir = "$unixpath/bin/diamond-0.9.29";
	my $transeq_dir = "$unixpath/bin/EMBOSS-6.6.0/emboss";	
	
# initialize CGI 
 	my $cgi = CGI->new;
 	my $job_id = $cgi->param('job_id');
 	my $database = $cgi->param('config_filename');
 	my $config_filename = "$unixpath/cgi-bin/$database";
	
# relative paths
	my $pks_script_dir = "$unixpath/scripts/";
	#my $hmm_pattern_dir = "$unixpath/hmm_models/";
	my $temp_dir = "$unixpath/tmp2";	
	my $data_dir = "$temp_dir/$job_id";	
	my $template_dir = "$webpath/templates/";		
	my $site_log = "$temp_dir/pksdb.log";
	my $INFILE = "$data_dir/$job_id.fasta";
	my $custom_page_header = "$unixpath/napdos_templates/napdos2_header";	
	my $custom_page_footer = "$unixpath/napdos_templates/napdos2_page_footer";
	
 # check for cgi input - otherwise might get these from a param file	
 	my $ref_seq_file = $cgi->param('ref_seq_file');
 	my $evalue = $cgi->param('evalue') || "1e-8"; 	
	my $min_matchlength = $cgi->param('min_matchlength') || "200";
	my $path_blast_evalue = $cgi->param('path_blast_evalue') || "1e-8";
	my $domain_type = $cgi->param('domain_type') || "KS";
	my $num_input_seqs = $cgi->param('num_input_seqs') || "0";
	my $compare = $cgi->param('compare') || "DNA";
	my $first_seq_len = $cgi->param('first_seq_len') || "0";
# probably can be removed (also from CGI sent to next scri0t
	my $max_hits =  $cgi->param('max_hits');
 	my $tree_hits =  $cgi->param('tree_hits');
# Prepare to write HTML response (site-wide templates)
	my $header = &get_file_text($custom_page_header);	
	$header =~ s/\$webpath/$webpath/g;	
	my $footer = &get_file_text($custom_page_footer);			
				
# Start writing HTML response to screen
	print $cgi->header;
	print $cgi->start_html(-title=>'Search Results',
							-style=>{'src'=>"../css/main.css"});
	print $header;	

# Start writing HTML response to file
	my $genomic_locations_page = "$job_id"."_genomic_locations.html"; # send delayed output here
	open (HTML_OUTPUT, ">$data_dir/$genomic_locations_page") or warn "couldn't open $genomic_locations_page\n $!\n";
	my $HTML_OUTPUT_text = ""; #&start_HTML_OUTPUT(); # this is an extra copy
	print HTML_OUTPUT "$HTML_OUTPUT_text";		

# also, write a tab-delimited table version
	my $genomic_locations_tabfile = "$job_id"."_genomic_locations.tab";
	open (TAB_OUTPUT, ">$data_dir/$genomic_locations_tabfile") or warn "couldn't open $genomic_locations_tabfile\n $!\n";
	
# open logfiles to append messages
	open (SITE_LOG, ">>$site_log") or warn "couldn't open logfile $site_log for appending\n $!\n";
	
# check to make sure all required accessory programs are available 
	&check_programs;
	
# Get parameter info from job_id file, if available, instead of CGI input
	my $param_info_file = "$data_dir/$job_id.params";
	my %saved_params = ();
	if (-s "$param_info_file")
	{
		%saved_params = &read_params_file("$param_info_file");
		$ref_seq_file = $saved_params{ref_seq_file};
		$domain_type = $saved_params{domain_type};
		$config_filename = $saved_params{config_filename};
		$max_hits = $saved_params{max_hits} || 1;
		$tree_hits = $saved_params{tree_hits} || 5;
		$evalue = $saved_params{evalue} || "1e-8";
		$job_id = $saved_params{job_id};
		$min_matchlength = $saved_params{min_matchlength} || "200";
		$path_blast_evalue = $saved_params{path_blast_evalue} || "1e-8";
		$num_input_seqs = $saved_params{num_input_seqs};
		$first_seq_len = $saved_params{first_seq_len};	
	}
	else
	{
		#&send_err("Can't find input data for JOB $job_id $param_info_file");
	}
	
# check for blast reference files	
	my $database_fasta = "$unixpath/fasta/$ref_seq_file";
	unless (-s $database_fasta)
	{
		&send_err("Couldn't find blast database file $database_fasta.");
	}
					
# Validate/parse input format (@sequences is list of sequence object references)
# check_filename exits if finds illegal characters (|?`> : ; < ..)
# convert dos and mac line endings to unix (not very efficiently)

	my $t_start_validate = new Benchmark;
	my @input_lines = ();
	if (length $INFILE > 2)
	{
		&check_filename($INFILE);
		@input_lines = &get_seqs_from_file($INFILE);
	}
	else
	{
		&send_err("Can't find temp infile $INFILE - line 150 genome_search.cgi");
	}
		
	my $seqref = \@input_lines;
	my @sequences = &get_seqs($seqref);
	unless (scalar @sequences > 0)
	{
		&send_err("empty input sequence array genome_search line 160");
	}
	my $sequence_type = &check_sequence_type(\@sequences);	
	my $t_end_validate = new Benchmark;
	
	my $validate_time = timediff ($t_end_validate, $t_start_validate);
	print SITE_LOG "\tVALIDATE INPUT:", timestr($validate_time, 'nop'), "\n";

# write temp file with unique name containing user input fasta sequences
	my $input_fasta_filename;
	my $input_fna_filename = "$data_dir/$job_id"."_orig.fna";
	
# this section makes sure can find original untrimmed input versus trimmed input 
	if ($sequence_type eq "nucleic")
	{
		$input_fasta_filename = $input_fna_filename;
	}	
	else
	{
		&send_err ("Problem with input sequence data (sequence_type=$sequence_type. \ngenomic search requires nucleic acid fasta format.");
	}

# set up container to hold query hit sequences for each input sequence
# make sure that USER hasn't entered sequences with duplicated id numbers
	my %query_hit_tally = ();	# key = query sequence_id  # value = num matches
	open(INPUT_FASTA, ">$input_fasta_filename") or warn "can't open temp file $input_fasta_filename for writing.\n $!";
	
	foreach (@sequences)
	{
		$query_hit_tally{$_->{ID}} = 0;
		print INPUT_FASTA "$_->{annotation}\n";
		print INPUT_FASTA "$_->{sequence}\n";
	}	
	close INPUT_FASTA;	
	
# translate nucleic acid to amino acid using EMBOSS TRANSEQ
# note: transeq flag "-sformat pearson" ensures that input sequence id numbers are not truncated
	my $t_start_transeq = new Benchmark;
	my $transeq_outname = "$data_dir/$job_id"."_orig.transeq";
	my $multifasta_outname = "$data_dir/$job_id"."_transeq.multifasta";
	my $min_transeq_matchlength = $min_matchlength;	
		
	if ($compare eq "DNA")
	{		
		my $cmd = "$transeq_dir/transeq -sequence $input_fna_filename -frame 6 -table 11 -sformat pearson -outseq $transeq_outname";
			my $response = `$cmd`;
		unless (-s $transeq_outname)
		{
			print SITE_LOG "No transeq file written:\n<br/>debug<br/>$cmd";
		}		
		
		$cmd = qq($pks_script_dir/transeq_to_multifasta.pl $transeq_outname $min_transeq_matchlength > $multifasta_outname);
			$response = `$cmd`;	
			
	}
	my $t_end_transeq = new Benchmark;
	
	my $transeq_time = timediff ($t_end_transeq, $t_start_transeq);
	print SITE_LOG "\tTRANSEQ TRANSLATION:", timestr($transeq_time, 'nop'), "\n";

# run diamond blast search (on orfs) to create candidate orflist	
	my $diamond_db = $database_fasta;
	$diamond_db =~ s/.faa/.dmnd/;

	my $blastp_result_filename = "$job_id".".blastp";
	my $program = "blastp";
	my $blastp_filtered_filename = "$job_id"."_filtered.blastp";
	
my $t_start_blastp = new Benchmark;	
# more matches in raw blastp to support later tree building
# diamond parameter "--length $min_matchlength" should filter out short matches, but doesn't seem to work
	my $blastp_raw_filename = "$blastp_result_filename".".raw";	
	my $blastp_command2 = "$diamond_dir/diamond $program -d $diamond_db -q $multifasta_outname -e $evalue --max-target-seqs 0 -p $num_processors --quiet -o $data_dir/$blastp_raw_filename";
		#print "$blastp_command2 <br />";
	my $blast_output = `$blastp_command2`;
	
# select top non-overlapping blastp matches for table
# DO trim metagenomes with large number of sequences (just not PCRs)
	if ($num_input_seqs	< $max_domain_hits || $min_matchlength > 150) # full pathway if input was reasonable size
	{	
		my $overlap_distance = 10; # num of nucleotides for alignment "slop" to detect a single domain
		# on short PCR fragments, just want the top match
		if ($min_matchlength < 100)
		 {
		 	$overlap_distance = 50
		 }
 		my $tophits_cmd = "$pks_script_dir/select_non-overlapping_tophits.pl $data_dir/$blastp_raw_filename $overlap_distance >$data_dir/$blastp_result_filename";			
  		$blast_output = `$tophits_cmd`;
	}
	else # for large num of input sequences (e.g. large PCR), just get tophit 
	{
		my $tophits_cmd = "$pks_script_dir/select_tab_first.pl $data_dir/$blastp_raw_filename  >$data_dir/$blastp_result_filename";
		#print "$tophits_cmd <br />";
		$blast_output = `$tophits_cmd`;	
	}
	
	my $blank_line_removal = `sed -i '/^$/d' $data_dir/$blastp_result_filename`;	

my $t_end_blastp = new Benchmark;
my $diamond_blastp_time = timediff ($t_end_blastp, $t_start_blastp);
	
# Filter blastp results for minimum align length (eliminate short matches)
  my $filter_blastp_cmd = "$pks_script_dir/select_tab_quant_greater_than_v2.pl $data_dir/$blastp_result_filename 3 $min_matchlength > $data_dir/$blastp_filtered_filename";
	`$filter_blastp_cmd`;
	`mv -f $data_dir/$blastp_filtered_filename $data_dir/$blastp_result_filename`;
# debug
	#print "$filter_blastp_cmd <br />";
	
my $t_end_blastp = new Benchmark;
my $diamond_blastp_time = timediff ($t_end_blastp, $t_start_blastp);
	print SITE_LOG "\tDIAMOND BLASTP:", timestr($diamond_blastp_time, 'nop'), "\n";	
	# debug
		#print "$blastp_command <br />";
	
	my $t_start_BLASTP_COORDS = new Benchmark;
	my $num_lines_blast_output = `wc -l $data_dir/$blastp_result_filename` || 0;
	chomp $num_lines_blast_output;	
					
# check to see if can exit early:
	if (length $num_lines_blast_output < 1)
	{		
		&send_err ( "Zero candidate domains found in preliminary screen.");
	}
	
# tally num of cand orfs/sequence and number of cand_domains/ORF from blastp results table
# parse blastp results to get protein coordinates into file
	my $blast_result_file = "$data_dir/$blastp_result_filename";
	my $blastp_cand_coords_filename = "$data_dir/$job_id"."_blastp_cand_coords";
	my %blastp_domain_cands_v2 = (); # keeps track of amino acid coordinates
		# key = coordinated-appended orf_id (parent_frame_orf#_start_end)
		# value = num hits. - should I change this to an object reference?
	my %blastp_cand_orf_ids = (); # key = 	 orfname  #value = num hits
		#keep track of which ones have hits and which do not
	my %blastp_domain_cands = ();	#key = orfname value = object reference 
		# DON'T USE THIS - doesn't work if there are multiple domains per orf
		#fields: cand_id, orf_id, orf_num, frame_num, parent_id, type, 
		# trim_aa_start, trim_aa_end, trim_nuc_start, trim_nuc_end
	my $current = "";	#placeholder for new objects	
	
	open(BLASTP_OUT, "<$blast_result_file") or die "can't open temp file $blast_result_file for reading.\n $!";
	open(BLASTP_COORDS, ">$blastp_cand_coords_filename") or die "can't open temp file $blastp_cand_coords_filename for writing.\n $!";
		
	while(<BLASTP_OUT>)
	{
		chomp;
		#push @blastp_lines, "$_";
		my @tmp = split '\t', $_;
		
		my $queryname = $tmp[0];
		my $db_match_name = $tmp[1];
		my $start = $tmp[6] || 0 ;
		my $end = $tmp[7] || 0;
		my $new_id = join "_", ("$queryname", "$start", "$end"); # query plus start and end coords		
	
	# get parent id from translated (transeq) blastp match name
		my ($parent_id, $orf_num, $frame_num);
		my @tmp2 = split "_", $queryname;
		$orf_num = pop @tmp2;
		$frame_num = pop @tmp2;
		$parent_id = join "_", @tmp2;
		
		#$blastp_domain_cands_v2{$new_id}++;
		$blastp_cand_orf_ids{$tmp[0]}++;	# keep track of number of candidate  domains per ORF
		push @tmp, $new_id;
		my $outstring = join "\t", @tmp;
		print BLASTP_OUT "$outstring\n";
		
		$query_hit_tally{$tmp[0]}++;
		$current = &new_domain("domain", $new_id);
			
		# object constructor fields
			$current->{domain_id} =$new_id; # includes amino acid coordinates	
			$current->{cand_id} = $queryname;
			$current->{type} = $domain_type;
			$current->{parent_id} = $parent_id;
			$current->{frame_num} = $frame_num;
			$current->{trim_aa_start} = $start;
			$current->{trim_aa_end} = $end;
			$current->{db_match_name} = $db_match_name;
			
		#$blastp_domain_cands{$queryname} = $current;
		#unless (exists $blastp_domain_cands_v2{$new_id}) # this isn't necessary
		#{
			$blastp_domain_cands_v2{$new_id} = $current;
		#}

		# create list for getting trimmed fastas (last column is new name, appending coordinates to frame# and orf#			
			#my @coords_info = ("$current->{orf_id}","$current->{trim_aa_start}","$current->{trim_aa_end}"); #, "$current->{cand_id}");
			my @coords_info = ("$current->{orf_id}","$current->{trim_aa_start}","$current->{trim_aa_end}", "$current->{domain_id}"); 

			my $coord_string = join "\t", @coords_info;
			print BLASTP_COORDS "$coord_string\n";
	}
	close BLASTP_OUT;
	close BLASTP_COORDS;
	my $t_end_BLASTP_COORDS = new Benchmark;	
	
	my $get_BLASTP_COORDS_time = timediff ($t_end_BLASTP_COORDS, $t_start_BLASTP_COORDS);
	print SITE_LOG "\tGET BLASTP_COORDS:", timestr($get_BLASTP_COORDS_time, 'nop'), "\n";

# debug  - how many lines in each hash?
#	my $num_blastp_domain_cands_v2 =scalar (keys %blastp_domain_cands_v2);
#	my $num_blastp_cand_orf_ids = scalar (keys %blastp_cand_orf_ids);
#	my $num_blastp_domain_cands = scalar(keys %blastp_domain_cands);
#	
#	print "num_blastp_domain_cands_v2=$num_blastp_domain_cands_v2 <br />
#	num_blastp_cand_orf_ids=$num_blastp_cand_orf_ids <br />
#	num_blastp_domain_cands=$num_blastp_domain_cands <br />";
	
# create fasta amino acid file containing all ORF candidates (untrimmed) that have at least one blast hit 
	my $t_start_getseq_cands = new Benchmark;
	my $blastp_cand_orf_ids_filename = "$data_dir/$job_id"."_blastp_cand_orf_ids";
	
	# create a smaller file than the entire input transeq to make subsequence retrieval faster
	open(OUTPUT, ">$blastp_cand_orf_ids_filename") or die "can't open temp file $blastp_cand_orf_ids_filename for writing.\n $!";
		
	foreach my $key (keys %blastp_cand_orf_ids) 
	{
		print OUTPUT "$key\n";
	}
	close OUTPUT;	
	my $cand_orf_faa_filename =  "$data_dir/$job_id"."_cand_orfs.faa";	
	my $cmd = "$pks_script_dir/getseq_multiple.pl $blastp_cand_orf_ids_filename $multifasta_outname > $cand_orf_faa_filename";

	`$cmd`;
	#debug
		#print "line 272 $cmd	<br /><br />";
	my $t_end_getseq_cands = new Benchmark;
	
	my $getseq_time = timediff ($t_end_getseq_cands, $t_start_getseq_cands);
	print SITE_LOG "\tGETSEQ CANDS:", timestr($getseq_time, 'nop'), "\n";
		
# get parent (nucleotide) id numbers for blast matches to ORFS (for faster subsequence selection on later blastx)
# put id numbers into file, select fasta sequences from these parents
# use hash so just get one copy of each parent id, even if it has multiple matches
# need these for getting candididate orfs (subsequences) and blastx later

	my $parent_match_ids_filename = "$data_dir/$job_id"."_parent_match_ids";
	my %parent_match_id_list =(); #key = transeq id  #value = original nucleotide id
	open(INPUT, "<$blastp_cand_orf_ids_filename") or die "can't open temp file $blastp_cand_orf_ids_filename for reading.\n $!";
	open(OUTPUT, ">$parent_match_ids_filename") or die "can't open temp file $parent_match_ids_filename for WRITING.\n $!";
		
	while (<INPUT>)
	{
		my @tmp = split "_", $_;
		my $num_fields = scalar @tmp;
		if ($num_input_seqs	< $max_domain_hits)
		{
			#my $end_num = pop @tmp;
			#my $start_num = pop @tmp;
		}
		my $orf_num = pop @tmp;
		my $frame_num = pop @tmp;
		my $orig_id = join "_", @tmp;
		$parent_match_id_list{$orig_id} = $orig_id;
	}
	close INPUT;
	
	#just get one line for each parent match
	foreach my $key (keys %parent_match_id_list)
	{
		print OUTPUT "$key\n";
	}	
	close OUTPUT;
	
 # use modified id numbers to select candidate parent fastas for getting coordinates
 # for multi-contig genomes/metagenomes, reduces number of source sequences for getting coordinates
	my $t_start_getseq = new Benchmark;
	my $parent_match_fasta_filename = "$data_dir/$job_id"."_parent_match_seqs.fna";
	 	
	my $select_cmd = "$pks_script_dir/getseq_multiple.pl $parent_match_ids_filename $data_dir/$job_id.fasta  > $parent_match_fasta_filename";
	#print "$select_cmd <br />";
	system($select_cmd);	
	
	my $t_end_getseq = new Benchmark;	
	my $getseq_time = timediff ($t_end_getseq, $t_start_getseq);
	#print SITE_LOG "\tGET MATCH SEQS:", timestr($getseq_time, 'nop'), "\n";	

# create file of trimmed orfs amino acid sequences from diamond BLASTX results
# but don't trim for pcr sequences with huge numbers of hits
# because blastx is too large and slow for server
	my $t_start_trim = new Benchmark;
	my $trimmed_blast_cands_faa_filename =  "$data_dir/$job_id"."_trimmed_blast_"."$domain_type"."domain_cands.faa";
	my $diamond_nuc_coords_filename = "$data_dir/$job_id"."_diamond_trimmed_nuc_coords";

	if ($num_lines_blast_output > 0)
	{					
		$cmd = "$pks_script_dir/get_subsequence.pl $blastp_cand_coords_filename $multifasta_outname > $trimmed_blast_cands_faa_filename";
		`$cmd`;
		#print "$cmd <br /><br />";				
	}
	
# Unless $num_matches < $max_domain_hits, get nuc coords of trimmed orfs by blastx vs orig.fna file, 
	my $getseq_filename = "$data_dir/$job_id"."_trimmed_nucleotide_getseq";
	my $trimmed_blast_cands_nuc_filename =  "$data_dir/$job_id"."_trimmed_nucleotide_"."$domain_type"."domain_cands.fna";
	my $parsed_diamond_blastx_filename = "$data_dir/$job_id"."_trim_cands_diamond.coords";
	#gi|14794889|gb|AF357202.1|	gi|14794889|gb|AF357202.1|_3_18_59_483	66102	67376
	my $num_matches = scalar (keys %blastp_domain_cands_v2);		
	unless ($num_matches > 0)
	{
			&send_err("No blast matches found for $domain_type domains\n");
	}

# get the coordinates only if number of matches isn't too high
	
	if ($num_matches > 0 && $num_matches < $max_domain_hits)
	{
		#update hashref for blastp_domain_cands
		my $updated_hashref = &get_nucleotide_coords_diamond ($diamond_nuc_coords_filename,\%blastp_domain_cands_v2, $parsed_diamond_blastx_filename);

		%blastp_domain_cands = %$updated_hashref;  
			#don't update blastp_domain_cands_v2 - otherwise gets rid of multiple domains per orf
		#my $printstring = &hashref_to_string(\%blastp_domain_cands);
		#print "line 418 <br >$printstring<br >";
		my $getseq_txt = `cut -f1,2,3,4 $parsed_diamond_blastx_filename`;
		# can't use unix command line to get columns in order 1,3, 4,2 
		#	but need column 2 to synchronize nt domain names with protein domain names
		# for trimmed output sequences - so revise lines below
		
		$parsed_diamond_blastx_filename = "$data_dir/$job_id"."_trim_cands_diamond.coords";
	#gi|14794889|gb|AF357202.1|	gi|14794889|gb|AF357202.1|_3_18_59_483	66102	67376
	# need to add a third column to getseq file with new (unique) subsequence names
	# appends coordinates to original source id so each new sequence has a unique id number
		my @getseq_lines = split "\n", $getseq_txt;
		my @revised_getseq_lines = ();
		foreach my $line (@getseq_lines)
		{
			my @tmp = split "\t", $line;
			#my $newname = join "_", ("$tmp[0]","$tmp[1]","$tmp[2]");
			my $parent_id = $tmp[0];
			my $domain_id = $tmp[1];
			my $start = $tmp[2];
			my $end = $tmp[3];
			#my $new_line = join "\t", ($line,$newname);
			my $new_line = join "\t", ($parent_id,$start,$end, $domain_id);
			push @revised_getseq_lines, ($new_line); 
		}
		$getseq_txt = join "\n", @revised_getseq_lines;
		&text_to_file("$getseq_txt",$getseq_filename);
		
		$cmd = "$pks_script_dir/get_subsequence.pl $getseq_filename $input_fna_filename > $trimmed_blast_cands_nuc_filename";
			#print "getting trimmed nucleotide sequences from genome_search2_diamond.cgi <br />";
			#print "$cmd <br /><br />"
		`$cmd`;
				
			
		# this would be a good place for an html spinner, while command is executing	
	}	
		else # don't try to trim if there are too many sequences
		{
			my $cmd = `cp $parent_match_fasta_filename $trimmed_blast_cands_nuc_filename`;
		}
					
	
	my $t_end_trim = new Benchmark;
	my $trim_time = timediff ($t_end_trim, $t_start_trim);
	print SITE_LOG "\tTRIM SEQS:", timestr($trim_time, 'nop'), "\n";
	
# pass results to display function
# file will get written in &display_domain_stats function 
	my $num_blast_matches = scalar (keys %blastp_domain_cands_v2);
	my $domain_stats = "";
	if ($num_blast_matches > 0)
	{
		if ($num_lines_blast_output <  $max_domain_hits)
		{
	# print full domain stats to HTML	
		$domain_stats = &display_domain_stats2 (\%blastp_domain_cands_v2, "Diamond BLASTX");
		print HTML_OUTPUT qq($domain_stats);
		
	# don't print stats to screen, but do print "continue" link
		$domain_stats = &display_pcr_stats(\%blastp_domain_cands_v2, "Diamond BLASTX");
		print qq($domain_stats);		
		
	# offer opportunity to download tab-delimited file of nucleotide locations
			my $download_path = "$webpath/tmp2/$job_id/$genomic_locations_tabfile";			
			print HTML_OUTPUT qq(
			<ul><li>Right-click to <a href="$download_path">DOWNLOAD</a> a table of nucleotide match 
			locations in tab-delimited format.</li></ul>);
			print qq(
			<ul><li>Right-click to <a href="$download_path">DOWNLOAD</a> a table of nucleotide match  
			locations in tab-delimited format.</li></ul>);
		}
		else  # too many matches to get nucleotide positions by blastx in a reasonable time 
		{
			$domain_stats = &display_pcr_stats(\%blastp_domain_cands_v2, "Diamond BLASTX");
			print HTML_OUTPUT qq($domain_stats);
			print qq($domain_stats);			
		}
		
	# offer immediate opportunity to search cand domains via next CGI
		my $further_anal_html = &call_blast_html($num_lines_blast_output,$trimmed_blast_cands_faa_filename);
		print $further_anal_html;
	
	# don't include the link to redo search on stored html 
		#print HTML_OUTPUT $further_anal_html;
	}
	else
	{
		print HTML_OUTPUT "No candidate domains found.";	
	}
						
# clean up temp files
	my @tmpfile_list = (
	#"$input_fna_filename",
	"$transeq_outname",
	"$multifasta_outname",	
	#"$blast_result_file",
	#"$blastp_cand_coords_filename",
	#"$blastp_cand_orf_ids_filename",
	#"$cand_orf_faa_filename",
	#"$trimmed_blast_cands_faa_filename",
	);
	
	foreach my $file (@tmpfile_list)
	{
		unlink $file;
	}
	
# Finish HTML response (custom ending)
	{
		print "$footer";
		#print HTML_OUTPUT "$footer";
	}
	
# close open files
	close SITE_LOG;
	close HTML_OUTPUT;
	close TAB_OUTPUT;

############################################
# SUBROUTINES
############################################

sub send_err
{
	
	my ($err_msg) = @_;
	print qq(<p>Sorry, your request could not be processed.</p>);
	print "$err_msg";	
	
	print SITE_LOG "\tERROR: $err_msg\n";
	print HTML_OUTPUT "\tERROR: $err_msg\n";
	
	print $cgi->end_html;
	print HTML_OUTPUT $cgi->end_html;
	
	close SITE_LOG;
	close HTML_OUTPUT;
	unlink "$data_dir/$genomic_locations_page";
	exit (0);
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

sub read_params_file
{
    # returns a hash, key = field name, value = value
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

sub get_seqs_from_file
{
	my ($INFILE) = @_;
	&check_filename($INFILE);

	my @input_lines = ();
	my @stats = stat($INFILE);
	my $filesize = $stats[7];
	my $filesize_string = $stats[7]/1000000;
	$filesize_string .=" MB";
	my $limit = 10000000; #10 MB	
#debug	
	$limit  = 1000000000; #1000 MB
	#$limit  = 2500000000; #2.5 GB
	
	my $limit_string = $limit/1000000;
	$limit_string .=" MB";
	
	open (FILEHANDLE, "$INFILE") or die "can't open debug filehandle\n";	
	while (my $line = <FILEHANDLE>)
	{
		chomp $line;		
		if ($filesize > $limit)
		{
			send_err("Input file size ($filesize_string) exceeds maximum ($limit_string).");
		}
		push (@input_lines, $line);
	}
	
# debug
	#my $num_input_lines = scalar @input_lines;
	#print "num_input_lines=$num_input_lines <br />";
	
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

sub get_seqs 
{
	my ($seqfile_name) = @_;
	my @seqlines = @$seqfile_name;

	my $counter = -1;
	my @seqs = ();	#a list of sequence object references
	my $limit = 500000;  # this should have already been set in process_request_napdos2.cgi
	
	if (scalar @seqlines <2)
	{
		send_err("no sequences entered");
	}
	my $current = "";	# sequence object reference		
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
		# these steps are really slow for large input sequences 
			#$seqs[$counter]->{sequence} =~ s/ //g; #remove any internal spaces
			#$seqs[$counter]->{sequence} =~ s/\r//g; #remove any internal linefeeds
			#$seqs[$counter]->{sequence} =~ s/\f//g; #remove any internal carriage returns						
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

sub new_domain 
{
	my ($className, $param) = @_;
 	my $self = {};
  	bless $self, $className;
  	 	
   $self->{cand_id} = $param;
 
 # If orig was nucleic, split into sections and re-assign values where necessary
	my ($parent_id, $orf_id, $orf_num, $frame_num, $frame_orfnum, $start, $end);
	my @tmp =split "_", $self->{cand_id};
	$end = pop @tmp;
	$start = pop @tmp;
	if ($sequence_type eq "nucleic")
	{			
		$orf_id = join "_", @tmp;
		if ($num_input_seqs	< $max_domain_hits) #for large # of input seqs, this never got aded
		{
			#$end = pop @tmp;
			#$start = pop @tmp;
		}
	# get parent id by popping of orf num and frame num	 
		$orf_num = pop @tmp;
		$frame_num = pop @tmp;
		$parent_id = join "_", @tmp;
		$self->{orf_id} = "$orf_id"; 
		 $self->{orf_num} = "$orf_num";
		 $self->{frame} = "$frame_num";
		 $self->{parent_id} = "$parent_id";
	}       

  	return($self);
}
sub check_sequence_type
{
	my ($seqlist_ref) = @_;
	my @seqlist = @$seqlist_ref;
	my ($seqstring, $seqtype);
	foreach (@seqlist)
	{
		my $seq_chars = $_->{sequence};
		$seq_chars =~ s/\s//g;
		$seqstring .= "$seq_chars";
	}
	my $seqlength = length $seqstring;
	
	my $num_ATCGN = (uc$seqstring =~ tr/ATUCGN//);
		my $num_aminos = (uc$seqstring =~ tr/EFILPQ//) || 0;
	
	if  ($num_ATCGN eq  $seqlength || $num_aminos < 1)
	{
		$seqtype = "nucleic";
		return $seqtype;
	}			
	elsif ($num_aminos > 1)
	{
		$seqtype = "protein";
		&send_err("Data input problem.<br> Fasta format nucleic acid sequence required for genomic analysis.");	
		#return $seqtype;
	}	
	else
	{
		$seqtype = "ambiguous";
		&send_err ("unable to determine sequence type genome_search.cgi line 756");
		return $seqtype;
	}
 
		if ($compare eq "protein")
		{
	 		unless ($seqtype eq "protein")
	 		{
	 			&send_err("Data input problem.<br> Fasta format amino acid sequence required for protein comparison.");
	 		}
	 	}
	 	if ($compare eq "DNA")
	 	{
	 		unless ($seqtype eq "nucleic")
	 		{
	 			&send_err("Data input problem.<br> Fasta format nucleic acid sequence required for DNA comparison.");
	 	 	}
	 	}
	
	return $seqtype;
}

sub parse_blastx_matches
{
	my ($blastx_filename) = @_;
	my $parsed_blastx = "";
	unless (-s $blastx_filename)
	{
		&send_err ("couldn't find blastx file $blastx_filename");
	}
	my $num_lines =0;
	my $num_hits = 0;
	open (BLASTX, "<$blastx_filename") or warn "couldn't open blastx file $blastx_filename for reading";
	{		
		while (<BLASTX>)
		{	
			chomp;
			$num_lines++;
			my @tmp = split "\t", $_;
			next unless (scalar @tmp > 7);
			my $query = $tmp[0];
			my $subject = $tmp[1];
			my $pct_id = $tmp[2];
			my $aln_len = $tmp[3];
			my $qstart = $tmp[7];
			my $quend =$tmp[8];
			next unless ($subject =~ /$query/);
			my @saved_fields = ($query, $subject,$qstart,$quend);
			my $saved_text = join "\t", @saved_fields;
			$parsed_blastx .= "$saved_text\n";
			$num_hits++;
		}
	}
	print "$num_hits saved from $num_lines in parse_blastx_matches<br />";
	close BLASTX;
	return $parsed_blastx;
}

sub get_nucleotide_coords_diamond
{
	my ($nuc_coords_filename, $domain_hashref, $parsed_diamond_blastx_filename) = @_;
	my %domain_objects = %$domain_hashref;  # domain objects should be the same as blastp_domain_cands
		
# Run blastx of orf candidates versus parent sequences
	# make sure this script waits for the command to be finshed using the "system" perl function
	# need to specify tmp file for formatdb, otherwise will fail on execution on web 
		my $diamond_formatdb_filename = "$data_dir/$job_id"."_trim_faa".".dmnd"; 
		my $diamond_format_cmd = "$diamond_dir/diamond makedb --in $trimmed_blast_cands_faa_filename -d $diamond_formatdb_filename -t $data_dir -p2";

		# use system() call to prevent race condition, to ensure that diamond makedb finishes before starting blastc
		my $test_output1 = system($diamond_format_cmd);					
				#print "$diamond_format_cmd <br /><br />";
			
		my $diamond_blastx_filename = "$data_dir/$job_id"."_trim_cands_diamond.blastx";
		my $diamond_blastx_cmd = "$diamond_dir/diamond blastx -d $diamond_formatdb_filename -q $parent_match_fasta_filename -p $num_processors -e $evalue --id 100 --max-target-seqs 0 --quiet -t $data_dir -o $diamond_blastx_filename";
			#print "$diamond_blastx_cmd <br />";
			system($diamond_blastx_cmd);
			
		# select lines from blastx file where identity = 100%  AND the parent seq matches the child
		# sort in decreasing order of bitscores, so don't get duplicate hits for same candidate orf
			my $parsed_diamond_blastx_filename = "$data_dir/$job_id"."_trim_cands_diamond.coords";
			my $parse_blastx_cmd =qq(sort -nr -k12 $diamond_blastx_filename |grep '100' |cut -f1,2,7,8); # > $parsed_diamond_blastx_filename);		
			#print "$parse_blastx_cmd <br />";		
			my $parsed_blastx_txt = `$parse_blastx_cmd`;
		
		# get only the tophit for each trimmed orf (some might match >1 db sequence at 100% identity	
			my %uniq_orf_blastx_hits = ();  #key = amino acid orfname (2nd col of trim_cands_diamond.coords # value = blastx line  
 			my @matchlines = split "\n", $parsed_blastx_txt;
 			foreach (@matchlines)
			{
 				chomp;
				my @tmp = split "\t", "$_";
 				my $orf_id = $tmp[1];
 				unless (exists $uniq_orf_blastx_hits{$orf_id})
 				{
 					$uniq_orf_blastx_hits{$orf_id} = "$_";
 				}	
			}
 			
 			$parsed_blastx_txt = join "\n", (values %uniq_orf_blastx_hits);		
 			&text_to_file("$parsed_blastx_txt", "$parsed_diamond_blastx_filename");
			
		unless (-s $parsed_diamond_blastx_filename)
		{
			&send_err("couldn't find diamond trim coords file $parsed_diamond_blastx_filename\n");
		}
	
# parse blastx output to get diamond trim coords file 						
	my $cand;
	my %uniq_orf_matches = (); # key = orfname value = num of matches at 100% identity
	
	open (INPUT, "$parsed_diamond_blastx_filename")  or (&send_err ("can't open file $parsed_diamond_blastx_filename\n$!\n"));
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s*$/);
		my @tmp = split "\t", $_;
		#print "parsed line: $_<br />";
			# NZ_CP017599.1	NZ_CP017599.1_5_147	7972651	7971260
			# parentid    orfid   aa_start   aa_end

		my $query = $tmp[1];
		$uniq_orf_matches{$query}++;
		unless ($uniq_orf_matches{$query} < 2)
		{
			#print "WARNING: more than one 100% identical domain for $query<br />";
			next;
		}
		
		#NZ_CP017599.1_5_147
		my $qstart = $tmp[2]; #nucleotide start
		my $qend = $tmp[3]; #nucleotide end
		unless (defined $qstart && $qstart >-1)
		{
			print "WARNING: no start coords for $query\n";
			next;
		}	
		
		my $full_domain_id = "$domain_objects{$query}->{domain_id}" || "not_found";
		#print "query=$query\n";
		#print "full_domain_id = $full_domain_id\n";
		
		unless (exists $domain_objects{$query})	
		{
			print "couldn't find domain object for domain_objects domain_id =$full_domain_id <br /> 
			query=$query, get_nucleotide_coords subroutine\n";
			next;
		}
	# debug: print candidate info	
		#print "line 856 cand =$cand domain_objects{query}=$domain_objects{$query} (should be a hashref) <br />";
			$cand = $domain_objects{$query};  # this is the full name, including protein coordinates
			$cand->{trim_nuc_start} = $qstart;
			$cand->{trim_nuc_end} = $qend;
		#my $printstring = &hashref_to_string(\%domain_objects);		
		#print "line 859 $printstring";
		#print "line 858 query=$query qstart=$qstart qend=$qend <br \>";
		#print "line cand->{trim_nuc_start} =  $cand->{trim_nuc_start} <br \>";
	}
	# debug 
		my $printstring = &hashref_to_string(\%domain_objects);		
		#print "line 864 $printstring";
	close INPUT;
	return $domain_hashref;			
}			


sub start_HTML_OUTPUT
{
my $webpath = "https://npdomainseeker.sdsc.edu";
# start html, general page  header
my $text =qq(<?xml version="1.0" encoding="iso-8859-1"?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US"><head><title>Search Results</title>
<link rel="stylesheet" type="text/css" href="/css/main.css" />
</head><body>);

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
<li><a href="$webpath/napdos2/pathways_v2.html">Pathways</a></li>
<li><a href="$webpath/napdos2/contact_v2.html">Contact Us</a></li>
</ul>
<div class="clear"></div>
</div> <!--End Topmenu-->
</div><!--End Header-->
<div id="main" class="centered shadowmain">
<div id="maininner">);
	
	return $text;
}

sub finish_HTML_OUTPUT
{
# general footer	
	my $text = qq(
	<p></p>
	</div>
	</div> <!--End Main-->
	
	<div id="footer"></div> <!--End Footer-->
	</body>
	</html>			
	);
	return $text;

}

sub call_blast_html
{
	my ($num_matches, $blast_cand_file) = @_;
	my $compare = "protein";	
	my $text = qq(<hr /><p>Click on the link below to: 
		<ul>
		<li> <p>Identify sequences assigned to different domain categories</p></li>
		<li> <p>Download matched sequences in fasta format</p></li> 
		<li> <p>Build trees comparing matched sequences with each other and database examples</p></li> 
		</ul>
	
		<form enctype="multipart/form-data" action="$webpath/cgi-bin/domain_blast2_v2.cgi" method="post">
			<input type="hidden" name="ref_seq_file"  value=$ref_seq_file />
			<input type="hidden" name="max_hits"  value=$max_hits />
			<input type="hidden" name="tree_hits"  value=$tree_hits />
			<input type="hidden" name="compare"  value= $compare />
			<input type="hidden" name="job_id"  value=$job_id />
			<input type="hidden" name="infile"  value="$trimmed_blast_cands_faa_filename" />
			<input type="hidden" name='num_input_seqs' value= $num_input_seqs />
			<input type="hidden" name="num_nuc_matches"  value=$num_matches />
			<button class = "button" type="submit" >CONTINUE ANALYSIS</button>		
		</form>	
		);
		return $text;
}

sub display_domain_stats2
{
	my ($domain_cands_hashref, $algorithm) = @_;
	my $webpath = "https://npdomainseeker.sdsc.edu";

	my %domain_cands = %$domain_cands_hashref;
			
	my %parent_seq_tally = ();
	my $num_domain_cands = scalar (keys %domain_cands);

	foreach my $key (keys %domain_cands)
	{
		my $domain = $domain_cands{$key};
		my $parent_id = $domain->{parent_id};
		$parent_seq_tally{$parent_id}++;
	}
	my $num_matches =	scalar keys %domain_cands;

# start output	
	my $tabfile_output = "";
	my $text_output = qq(
	<form enctype="multipart/form-data" action="$webpath/cgi-bin/domain_blast2_diamond.cgi" method="post">	

<h2>Nucleotide Sequence Search Results</h2>
	
<ul>
	<li>$num_matches $domain_type domains were detected from $num_input_seqs input sequences.</li>
	</ul>
<hr />

	<h3>Domain Locations</h3> 	
	);
						
	unless (scalar keys %domain_cands > 0)
	{
		$text_output .=  "\nZero matches found using $algorithm search\n";
		return($text_output);
	}

# start table			
	$text_output .= qq(
	<script src="../../scripts/sorttable.js" type="text/javascript"></script>
	<table class = "sortable2" width = "90%"  cellpadding = "5" cellspacing = "0" border = "1">\n);

	my @headers2 = ("domain id", "parent seq", "start", "end");
	if ($sequence_type eq "nucleic")
	{
		push @headers2, "read frame";
		push @headers2, "cand_id";
	}
# header row for table	
	$text_output .= qq(<tr>);
	foreach my $header (@headers2)
	{
		$text_output .= qq(
		<th>$header</th>);
		$tabfile_output .= "$header"."\t";
	
	}
	$text_output .=  qq(\n\t</tr>);
	$tabfile_output .= qq(\n);
		
	my %parent_tally = ();
	my %outlines = (); #key=text_line (parent_id, start,stop,strain,child_id  value=nt_start_pos

# sort keys by domain number assigned previously, e.g. KS1, KS2 etc
	my @cand_id_list = keys %domain_cands;  #complete line
	my $num_domain_cand_hashkeys = scalar keys %domain_cands;
	my $line_num = 0;
	foreach my $cand_id (sort keys %domain_cands)
	{
		my $cand = $domain_cands{$cand_id}; # hashtag
		$line_num++;
		my @elements = (
			"$cand->{type}",
			"$cand->{parent_id}",);
		# prevent overwrite of parent sequence with coords by parent sequences without hits
		unless (defined $cand->{trim_nuc_start})
			{
				#print "line 1255 undefined trim_nuc_start $line_num $cand->{parent_id} $cand->{parent_id} start=$cand->{trim_nuc_start} stop=$cand->{trim_nuc_end} frame=$cand->{frame}<br />"; 
				next;
			}
		if ($sequence_type eq "nucleic")
		{				
			push  @elements, "$cand->{trim_nuc_start}";
			push  @elements, "$cand->{trim_nuc_end}";			
			push  @elements, "$cand->{frame}";
			#print "codeline 1016,$line_num $cand->{parent_id} start=$cand->{trim_nuc_start} stop=$cand->{trim_nuc_end} frame=$cand->{frame}<br />";		
		}
		else
		{
			push  @elements, "$cand->{trim_aa_start}";
			push  @elements, "$cand->{trim_aa_end}";		
		}
		push  @elements, "$cand_id";
		
		my $printstring = join "\t", @elements;
	# keep track of start/stop coords so can sort output		
		$outlines{$printstring} = "$elements[2]";
		#print "sorting hash printstring=$printstring <br /> contents=$printstring=$elements[2]<br />";
	}	
	my $num_outlines = scalar keys %outlines;

# sort by coordinates, flag overlapping
	my ($current_start, $current_stop, $current_frame, $left_end, $right_end);
	my $prev_left = 0;
	my $prev_right = 0;
	my $prev_frame = 0;
	my $prev_parent = "none";
	my $frameshift_tally = 0;
	my $line_num = 0;
	my %ordered_cand_names = ();
	
	my $total_tally = 0;
	foreach my $next (sort {$outlines{$a}<=>$outlines{$b}} keys %outlines)
	{
		my @tmp = split "\t", $next;
		my $cand_type = $tmp[0];
		my $parent_id = $tmp[1];
		$current_start = $tmp[2];
		$current_stop = $tmp[3];
		$current_frame = $tmp[4];
		$total_tally++;  #for this particular key, but key here is a whole line
		#print "line  1049 $next _start=$current_start stop=$current_stop frame=$current_frame<br />";	
	#	deal with sequences on antisense strand
		if ($current_stop > $current_start)		
		{
			$left_end = $current_start;
			$right_end = $current_stop;
		}
		else
		{
			$left_end = $current_stop;
			$right_end = $current_start;
		}
	
	# for genomic seqs, count total num KS rather than than for a specific parent sequence
		#my $cand_type = $tmp[0]; #defined_above
		#my $parent_id = $tmp[1]; #defined above
	
		$line_num++;
		#$parent_tally{$parent_id}++;
		#my $current_tally = $parent_tally{$parent_id};
		$tmp[0] = "$cand_type"."$total_tally";
		#print "line 1073 cand_type=$cand_type total_tally=$total_tally parent_id=$parent_id start=$current_start stop=$current_stop frame=$current_frame<br />";	
	# keep track of tally number associated with this domain
		#$cand->{tally_num} = $line_num;
		
# create simpler candidate name		
		#my $alt_cand_name = "$parent_id"."_$tmp[0]";
		my $alt_cand_name = "$job_id"."_$tmp[0]";
		$ordered_cand_names{$tmp[5]}=$alt_cand_name;
		#$tmp[5] = $alt_cand_name; # if use this, have to change fastas & all other cgi scripts as well
		
		my $printstring = join "\t", @tmp; 
		
	# check for overlap with previous
	# this needs to be modified so it finds overlaps on the same parent only
		if ($left_end <= $prev_right && $current_frame != $prev_frame && $parent_id eq $prev_parent)
		{
			$printstring .= "*";
			$frameshift_tally++;
		}	

	#	print coords to saved output file
		my @cells = split "\t", $printstring;
		my @file_output_info = ($cells[1],$cells[2],$cells[3],$cells[5] ); # this is html order
			# debug 
				# print "line 1014 printstring=$printstring<br />";
# 				print "     parent_seq=$cells[1]<br />";
# 				print "     start=$cells[2]<br />";
# 				print "     end=$cells[3]<br />";
# 				print "     new_name=$cells[5]<br />";			

		my $file_output_text = join "\t",  @file_output_info ;
		# get rid of asterisk about possible frame-shift
		$file_output_text =~ s/\*//;
		print NUC_COORDS "$file_output_text\n";
		# debug 
			#print "line 1026 printstring=$file_output_text<br />";
	#	make html table row for each candidate, convert printstring to html row		
		$text_output .= qq(<tr>);
		
# wrap long names, so table doesn't get too wide	
		$text_output .= qq(<style>
		td {
				word-break: break-all;
		}
		</style>
		); 
		
		
		foreach my $cell (@cells)
		{
			$tabfile_output .= "\t$cell";
			if (length $cell < 3 || $cell =~ /^d+$/)
			{
				$text_output .= qq (<td align = "center" >$cell</td>\n)
			}
			else
			{
				$text_output .= qq (<td align = "left" >$cell</td>\n)
			}
		}
		$tabfile_output .= "\n";
		$text_output .= qq(</tr>);	
		$prev_left = $left_end;
		$prev_right = $right_end;
		$prev_frame = $current_frame;
		$prev_parent = $parent_id;
	}	
	
	$text_output .= qq(</table>);
	if ($frameshift_tally>0)
	{
		$text_output .= "\n<p>* possible frameshift</p>\n";
		$tabfile_output .= "\n* possible frameshift\n";
	}
	
	# explain cand_id abbreviations
		$text_output .= "\n<p><b>Note:</b> Each potential domain candidate has been given 
a unique candidate ID number, based on parent sequence id, reading frame number (1-6), 
gene number within the reading frame, and trimmed nucleotide start and stop coordinates 
within the reading frame.</p>";
	
	print TAB_OUTPUT "$tabfile_output";
	return $text_output;
}

sub display_pcr_stats 
{
	my ($domain_cands_hashref, $algorithm) = @_;
	my %domain_cands = %$domain_cands_hashref;
	my %parent_seq_tally = ();
	my $num_domain_cands = scalar (keys %domain_cands);
	foreach my $key (keys %domain_cands)
	{
		my $domain = $domain_cands{$key};
		my $parent_id = $domain->{parent_id};
		$parent_seq_tally{$parent_id}++;
	}
	my $num_matches =	scalar keys %domain_cands;
	
#my $tabfile_output = "";
	my $text_output = qq(
	<form enctype="multipart/form-data" action="$webpath/cgi-bin/domain_blast2_diamond.cgi" method="post">	

<h2>Nucleotide Sequence Search Results</h2>
	
<ul>
	<li>$num_matches preliminary $domain_type domain candidates were detected in $num_input_seqs input sequences.</li>
	
	</ul>
		);
	
	
	return $text_output;
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


sub re_order_faa
{
	my ($array_ref, $faa_file) = @_;
	my @order_list = @$array_ref;
	my $revised_fasta = "";
	# convert faa_file into hash of sequence objects: key = id num, value = fasta sequece)
	# for each item in ordered list, select fasta object from hash
		# add fasta object to new_faa_string
	# return faa string
		
	return $revised_fasta;	
}

sub hashref_to_string
{
	my ($hashref) = @_;
	my %domain_objects = %$hashref;
	my $outstring = "";	
foreach my $key (keys %domain_objects)
{
	my $cand =$domain_objects{$key};
$outstring .= qq (cand_id = $cand->{cand_id}<br />
orf_id = $cand->{orf_id}<br />
orf_num = $cand->{orf_num}<br />
frame_num = $cand->{frame_num}<br />
parent_id = $cand->{parent_id}<br />
type = $cand->{type}<br />
trim_nuc_start = $cand->{trim_nuc_start} <br />
trim_nuc_end = $cand->{trim_nuc_end} <br />
 <br /> <br />);
}		
return $outstring;	

}




__END__
# command-line debug
	# my $command_line = -1;
# 	if ($command_line > 0)
# 	{
# 		my $localdir = "/Users/spodell/Desktop/PKS_db_2010/";
# 		$temp_dir = "$localdir/pipeline_dev/tmp/";
# 		$hmm_pattern_dir = "$localdir/pipeline_dev/hmm_models/";
# 		$pks_script_dir = "$localdir/pipeline_dev/";
# 		my $max_hits = "10";
# 		my $compare = "nucleic"; #"protein" or "nucleic"
# 		my $evalue = "1e-5";
# 		my $out_format = "table";	# "table" or ??
# 		my $INFILE = $ARGV[0] || "$localdir/test_input_sequences/mupirocin.fna";
# 	}
# 	unless (-s $INFILE)
# 	{
# 		print STDERR "can't find input sequences in file $INFILE\n";
# 		exit (0);
# 	}
# 	# get filename for output
# 	my @suffixlist = (".fasta", ".txt", ".fa", ".faa", ".fna", ".fsa");
# 	my ($basename,$path,$suffix) = fileparse($INFILE,@suffixlist);			
#  			 $job_id = $basename;
