#!/usr/local/bin/perl
# load_pks_tables_v3.pl
# Sheila Podell
# November 25, 2019

use warnings;
use strict;
use Getopt::Long;
use DBI;

# creates, populates MySQL tables for PKS/NRPS database:
# organisms, pathways, ref_genes, and domains

# input files;

# configuration file needed with database info in following format:
# 	[db_name]=pks_01
# 	[db_program]=mysql
# 	[db_host]=localhost
# 	[db_user]=username
# 	[db_user_password]=your_password
# 	[max_lines_per_packet]=2000

# tab-delimited datafile with fields in the following order:
# 	field num	field name
# 	0	path name
# 	1	cluster type
# 	2	public-private status
# 	3	pubmed ref id
# 	4	mibig_id
# 	5	example name
# 	6	example struct filename
# 	7	Species name
# 	8	strain name
# 	9	source_location
# 	10	organism description
# 	11	gene name
# 	12	domain_name
# 	13	domain type
# 	14	domain sub-type
# 	15	Genbank prot ID
# 	16	trim start
# 	17	trim end
# 	18	curator
#	19	db_vers
# 	20	in_ref_tree
#	21	align_subtree
#	22	domain_subtype_abbrev

# revised 082520 tab-delimited datafile with fields in the following order:
# 	field num	field name
# 0	path_name
# 1	cluster type
# 2	public_status
# 3	pubmed_ref_id
# 4	mibig_id
# 5	example_product_name
# 6	example structure file-name
# 7	species_name
# 8	strain_name
# 9	source location
# 10	organism desc
# 11	gene_prod_name
# 12	domain name
# 13	domain_type
# 14	domain_class
# 15	domain_subclass
# 16	gbk_prot_id
# 17	gene_trim_start
# 18	gene_trim_end
# 19	curator
# 20	db_vers
# 21	in_ref_tree
# 22	align_subtree
# 23	domain_subclass_abbrev


# will not over-write old tables unless -o option specified on command line	
	
# get command line arguments 
	my ($db, $datafile, $sth, $sql, $overwrite, 
		$update, $result,  $lines_per_packet, $config_filename, $tablename);
	my $USAGE;
	my $message = qq(
	Usage: $0  -i source_datafile -c config filename <-o allow overwrites>)."\n\n";
	GetOptions( 
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    'o'   => \$overwrite,		#allow table overwrites
			    'u'   => \$update,
				);	     
	if($USAGE || !$config_filename || !$datafile) 
	{
		print STDERR $message;
		exit(0);
	} 
	unless (-s $datafile)
	{
		print STDERR "  Couldn't find source datafile $datafile \n";
		exit(0);
	}
	unless (-s $config_filename)
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		exit(0);
	}
	
# get SQL database parameters from config file
	my %config_params= &read_config_file($config_filename);
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_port = $config_params{db_port};
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";

# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	if (exists $config_params{db_port} && defined $db_port)
	{
		$db_path .= "$db_port";
	}
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";
	
# check to see if any of tables already exist
	my @tablenames = ("organisms", "pathways", "ref_genes", "domains");	
	my $found_table = 0;
	print STDERR "\n";

foreach my $name (@tablenames)
{
	$sql = "show tables;";
	$sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	my @row = ();
	my $element = ();
	while (@row = $sth->fetchrow_array)
	{
		foreach $element (@row)
		{
			if ($element =~/$name/)
			{
				$found_table = 1;
				unless (defined $overwrite || defined $update)
				{
					print STDERR "  Table $name already exists.\n";
					print STDERR "  Specify -o option if you wish to overwrite.\n";
					print STDERR "  Specify -u option to update without overwriting.\n";
					exit(0);
				}
				last;
			}
		}
	}
	
# if over-write option specified, drop table from database.
	if (defined $overwrite)
	{
		print STDERR "   Over-writing table $name.\n";
		$sql = qq(DROP TABLE IF EXISTS `$name`);
		$result = $dbh->do($sql);
		
	}
	if (defined $update)
	{
		print STDERR "updating table $name.\n";
	}
}

# define lines/packet to revent buffer over-runs
	if (defined $overwrite || $found_table == 0)
	{
		$lines_per_packet = 2000;
	}
	elsif (defined $update)
	{
		$lines_per_packet = 1;
	}
		print STDERR "\n";

# Create organisms table
	$tablename = "organisms";
	if (defined $overwrite || $found_table == 0)
	{
		$sql = qq(CREATE TABLE `$tablename` (
		`organism_id` int(11) NOT NULL,
		`species_name` varchar(100) NULL,
		`strain_name` varchar(100) NULL,
	   `source_location` varchar(255)  NULL,
	   `desc` varchar(255)  NULL,
	   PRIMARY KEY  (`organism_id`),
	   INDEX (`strain_name`),
	   INDEX (`species_name`) 
		) ;\n);
		$result = $dbh->do($sql);
		if ($result)
		{
			print STDERR "   New table $tablename created.\n";
		} 
		else
		{
			print STDERR "Unable to create new table $tablename";
			exit(0);
		}
    }
 
# Create pathways table
	$tablename = "pathways";
	if (defined $overwrite || $found_table == 0)
	{
		$sql = qq(CREATE TABLE `$tablename` (
		`path_id` int(11) NOT NULL,
		`path_name` varchar(50) NOT NULL,
		`cluster_type` enum('pks','nrps','pks-nrps', 'FAS') NULL,
		`public_status` enum('public','private') NULL,
	   `pubmed_ref_id` varchar(50)  NULL,
	   `mibig_id` varchar(50)  NULL,
	   `example_product_name` varchar(100)  NULL,
	   `example_struct_filename` varchar(100)  NULL,	   
	   `organism_id` int(11) NULL,
	   PRIMARY KEY  (`path_id`),
	   INDEX (`path_name`),
	   INDEX (`organism_id`)
		) ;\n);	

		$result = $dbh->do($sql);
		if ($result)
		{
			print STDERR "   New table $tablename created.\n";
		} 
		else
		{
			print STDERR "Unable to create new table $tablename";
			exit(0);
		}
    }
    
  # Create ref_genes table
	$tablename = "ref_genes";
	if (defined $overwrite || $found_table == 0)
	{
		$sql = qq(CREATE TABLE `$tablename` (
		`ref_gene_id` int(11) NOT NULL,
		`gbk_prot_id` varchar(50) NULL,
		`gene_name` varchar(100) NULL,
	     PRIMARY KEY  (`ref_gene_id`),
	   	 INDEX (`gene_name`) 
		) ;\n);	
		$result = $dbh->do($sql);
		if ($result)
		{
			print STDERR "   New table $tablename created.\n";
		} 
		else
		{
			print STDERR "Unable to create new table $tablename";
			exit(0);
		}
    }
          
# Create domains table
	$tablename = "domains";
	if (defined $overwrite || $found_table == 0)
	{		
		$sql = qq(CREATE TABLE `$tablename` (
		`domain_id` int(11) NOT NULL auto_increment,
		`domain_name` varchar(50) NOT NULL,
		`domain_type` enum('KS','C') NULL,
		`domain_class` varchar(100) NULL,
		`domain_subclass` varchar(100) NULL,
	   `path_id` int(11)  NULL,
	   `ref_gene_id` varchar(20)  NULL,
	   `gene_trim_start` int(11) NULL,
	   `gene_trim_end` int(11) NULL,	   
	   `curator` varchar(20) NULL,
	   `in_ref_tree` enum('yes','no') NULL,
	   `align_subtree` varchar(50) NULL,
	   `domain_subclass_abbrev` varchar(50) NULL,
	   PRIMARY KEY  (`domain_id`),
	   INDEX (`domain_name`)
		) ;\n);	
		
		$result = $dbh->do($sql);
		if ($result)
		{
			print STDERR "   New table $tablename created.\n";
		} 
		else
		{
			print STDERR "Unable to create new table $tablename";
			exit(0);
		}
    }    
        
# parse input file, storing only good lines 
	my @goodlines = ();
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s*$/);
		# skip hearder
		next if ($_ =~ /.*name\t/);
	# check for correct number of fields
		my @fields = split "\t", $_;
		unless (scalar @fields > 16)
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			next;
		}
		push @goodlines, $_;
	}
	close INPUT;
	my $num_goodlines = scalar @goodlines;
	my $datafile_wc = `wc -l $datafile`;
 	$datafile_wc  =~ /(\d+).+/;
 	my $total_num = $1;
 	print STDERR "\n  Found $num_goodlines/$total_num usable lines of input data in $datafile\n\n";

# populate organisms table. Include one entry for each unique combo of species_name, strain_name
	$tablename = "organisms";
	$sql  = " INSERT INTO $tablename VALUES ";
	my $count  = 0;
	my %uniq_organisms = ();
	foreach (@goodlines)
	{
		chomp;
		my $line_num = $.;
		my @fields = split "\t", $_;
			
		my $species_name = $fields[7] || "NULL";
		my $strain_name = $fields[8] || "NULL";
		my $source_location = $fields[9] || "NULL";
		my $desc = $fields[10]  || "NULL";
		my $combo_key = "$species_name"."$strain_name";
		
		unless (exists $uniq_organisms{$combo_key})
		{
			$count++;
			$uniq_organisms{$combo_key} = $count;
			my $organism_id = $uniq_organisms{$combo_key};
			$sql .= "('$organism_id','$species_name', '$strain_name', '$source_location', '$desc'),\n";
		}						
	}	
	$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
	$sql .= ";\n";
	$result = $dbh->do($sql);
	print STDERR "  Loaded $count line(s) into $tablename table.\n";
	
# populate pathways table - use only unique path names, cross referenced to organism_id
	my %uniq_pathnames = ();
	$tablename = "pathways";
	$sql  = " INSERT INTO $tablename VALUES ";
	$count = 0;
	
	foreach (@goodlines)
	{
		chomp;
		my @fields = split "\t", $_;
			
		my $path_name = $fields[0];
		unless (defined $path_name && length $path_name > 0)
		{
			print STDERR "problem processing $path_name in line $. Skipping line\n$_\n";
			next;
		}				
		my $cluster_type = $fields[1] || "NULL";
		my $public_status = $fields[2] || "NULL";
		my $pubmed_ref_id = $fields[3] || "NULL";
		my $mibig_id = $fields[4] || "NULL";
		my $example_product_name= $fields[5] || "NULL";
		my $example_struct_filename= $fields[6] || "NULL";
		my $organism_id= $fields[-1] || "NULL";
		my $species_name = $fields[7] || "NULL";
		my $strain_name = $fields[8] || "NULL";
		
		unless (exists $uniq_pathnames{$path_name})
		{
			$count++;
			$uniq_pathnames{$path_name} = $count;
			my $path_id = $uniq_pathnames{$path_name};
		
		# figure out organism_id used in other table	
			my $combo_key = "$species_name"."$strain_name";
			my $org_id = $uniq_organisms{$combo_key};	
			$sql .= "('$path_id','$path_name','$cluster_type', '$public_status', 
				'$pubmed_ref_id', '$mibig_id','$example_product_name', 
				'$example_struct_filename', '$org_id'),\n";
		}
	}	
	$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
	$sql .= ";\n";
	$result = $dbh->do($sql);
	print STDERR "  Loaded $count line(s) into $tablename table.\n";

# populate ref_genes table - use only unique ref_genes, defined by both name and gbk id #
	$tablename = "ref_genes";
	$sql  = " INSERT INTO $tablename VALUES ";
	$count  = 0;
	my %uniq_genes = ();
	foreach (@goodlines)
	{
		chomp;
		my @fields = split "\t", $_;
			
		my $gbk_prot_id = $fields[16] || "NULL";
		my $gene_name = $fields[11] || "NULL";
		my $gene_key = "$gbk_prot_id"."$gene_name";
		
		unless (exists $uniq_genes{$gene_key})
		{
			$count++;
			$uniq_genes{$gene_key} = $count;
			my $gene_count = $uniq_genes{$gene_key};
			$sql .= "('$gene_count','$gbk_prot_id', '$gene_name'),\n";
		}						
	}	
	$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
	$sql .= ";\n";
	$result = $dbh->do($sql);
	print STDERR "  Loaded $count line(s) into $tablename table.\n";

# populate domains table (cross-referenced to ref_gene_id and path_id)
	$tablename = "domains";
	$count = 0;
	$sql  = " INSERT INTO $tablename VALUES ";
	$count = 0;
	
	foreach (@goodlines)
	{
		chomp;
		my @fields = split "\t", $_;
			
		my $path_name = $fields[0];
		my $path_id = $uniq_pathnames{$path_name};

		my $gbk_prot_id = $fields[16] || "NULL";
		my $gene_name = $fields[11] || "NULL";
		my $combo_key = "$gbk_prot_id"."$gene_name";
		my $ref_gene_id = "$uniq_genes{$combo_key}"; 	
		my $domain_name  = $fields[12] || "NULL";
		
		my $domain_type = $fields[13] || "NULL";
		my $domain_class = $fields[14] || "nd";
		my $domain_subclass = $fields[15] || "nd";
		
		unless (defined $domain_subclass && length "$domain_subclass" >0)
		{
			$domain_subclass = "nd";		
		}

		my $gene_trim_start = $fields[17] || "NULL";
		my $gene_trim_end= $fields[18] || "NULL";
		my $curator= $fields[19] || "NULL";
		my $in_ref_tree = $fields[21] || "NULL";
		my $align_subtree = $fields[22] || "NULL";
		my $domain_subclass_abbrev = $fields[23] || "NULL";
		
		$count++;
			
		$sql .= "('','$domain_name','$domain_type', '$domain_class', '$domain_subclass','$path_id', '$ref_gene_id', 
				'$gene_trim_start', '$gene_trim_end','$curator','$in_ref_tree','$align_subtree',
				 '$domain_subclass_abbrev'),\n";
	}	
	$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
	$sql .= ";\n";
	$result = $dbh->do($sql);
	print STDERR "  Loaded $count line(s) into $tablename table.\n";

# cleanup
	$sth->finish();
	$dbh->disconnect();
	
###########################
# SUBROUTINES
###########################

sub read_config_file
{
    my ($filename) = @_;
    open (INPUT, $filename)  or die "can't open DarkHorse config file $filename\n$!\n";
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


__END__
