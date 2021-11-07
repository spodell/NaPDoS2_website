#!/usr/local/bin/perl
# load_group_tables.pl
# Sheila Podell
# March 5, 2010

use warnings;
use strict;
use Getopt::Long;
use DBI;

# takes as input two files:

	# directory containing curated alignment file(s) in
	# tab-delimited file describing group:
		# group_name
		# align_filename
		# example domain name
		# description

# creates, populates MySQL tables for PKS/NRPS database:
# alignment groups, memberships
# needs to have domains table already loaded

# id numbers for pre-curated aligned sequences need to 
# match names already entered in "domains table"

# will not over-write old tables unless -o option specified on command line

# get command line arguments 
	my ($db, $datafile, $sth, $sql, $overwrite, $align_file_directory,
		$update, $result,  $lines_per_packet, $config_filename, $tablename);
	my $USAGE;
	my $message = qq(
	Usage: $0  -i source_datafile -d align_file directory -c config filename <-o allow overwrites>)."\n\n";
	GetOptions( 
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    "d=s" => \$align_file_directory,
			    'o'   => \$overwrite,		#allow table overwrites
			    'u'   => \$update,
				);	     
	if($USAGE || !$config_filename || !$datafile || !$align_file_directory) 
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
	unless (-s $align_file_directory)
	{
		print STDERR "  Couldn't find align_file_directory $align_file_directory \n";
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
	
	
# check to see if table(s) already exist
	my @tablenames = ("alignment_groups", "memberships");
	my %requirements = ("domains" => "1");
	my $found_table = 0;
	my $found_requirement = -1;
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
			if ($element =~/$name/ && exists $requirements{$name})
			{
				$found_requirement =1;
				next;
			}			
			elsif ($element =~/$name/)
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
	
	foreach my $key (keys %requirements)
	{
		unless ($requirements{$key} > 0)
		{
			print STDERR "  ERROR: Could not find pre-requisite table $name in database $db_name.\n";
			exit(0);
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


# Create alignment groups table
	$tablename = "alignment_groups";
	if (defined $overwrite || $found_table == 0)
	{
		$sql = qq(CREATE TABLE `$tablename` (
		`group_name` varchar(50) NOT NULL,
		`align_filename` varchar(50) NULL,
	   `example_domain_name` varchar(255)  NULL,
	   `group_desc` varchar(255)  NULL,
	   PRIMARY KEY  (`group_name`) 
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
 
# Create memberships table
	$tablename = "memberships";
	if (defined $overwrite || $found_table == 0)
	{
		$sql = qq(CREATE TABLE `$tablename` (
		`membership_id` int(11) NOT NULL auto_increment,
		`group_name` varchar(50) NOT NULL,
	   `domain_name` varchar(50)  NOT NULL,
	   PRIMARY KEY  (`membership_id`),
	   INDEX (`group_name`),
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
    
# parse input files, storing only good lines
	my @align_groups =();	#list of object references
	my $num_goodlines = 0;
	my $current = "";
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s*$/);
		# skip hearder
		next if ($_ =~ /.*name\t/);
	# check for correct number of fields
		my @fields = split "\t", $_;
		unless (scalar @fields >3)
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			next;
		}
		$current = &new_group("group",\@fields);					
		next unless (defined $current);
		push @align_groups, $current;

		$num_goodlines++;
		
	}
	close INPUT;
	
	my $datafile_wc = `wc -l $datafile`;
 	$datafile_wc  =~ /(\d+).+/;
 	my $total_num = $1;
 	print STDERR "\n  Found $num_goodlines/$total_num usable lines of input data in $datafile\n\n";

# populate groups table	
	$tablename = "alignment_groups";
	$sql = " INSERT INTO $tablename VALUES \n";
	my $count = 0;
		
	foreach my $next (@align_groups)
	{
		my $group_name = $next->{group_name};
		my $align_filename = $next->{align_filename};
		my $example_domain_name = $next->{example_domain_name};
		my $desc = $next->{desc};
								
		$sql  .= "('$group_name','$align_filename','$example_domain_name', '$desc'),\n";
		$count++;
	}
	$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
	$sql .= ";\n";
	$result = $dbh->do($sql);
	#print STDERR "$sql\n";
	print STDERR "  Loaded $count line(s) into $tablename table.\n";

# populate memberships table
	$tablename = "memberships";
	$sql = " INSERT INTO $tablename VALUES \n";
	$count = 0;
	
	foreach my $next (@align_groups)
	{
		my $group_name = $next->{group_name};
		my @id_list = @{$next->{domain_names}};	
		foreach my $domain_name (@id_list)		
		{		
			$sql  .= "('','$group_name','$domain_name'),\n";
			$count++;
		}
		
	}
	$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
	$sql .= ";\n";
	$result = $dbh->do($sql);
	#print STDERR "$sql\n";
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

sub new_group
 {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  
   $self->{group_name} = $properties[0];
   $self->{align_filename} = $properties[1];
   $self->{example_domain_name} = $properties[2];
    $self->{desc} = $properties[3];
  
	my $alignpath = "$align_file_directory"."/"."$properties[1]";
	unless (defined $properties[1] && -s "$alignpath")
		{
			print STDERR "ERROR: unable to find align file $alignpath for group  $self->{group_name}\n";
			return(undef);
		}
   
  
# parse align input file for domain names
  my @domain_names = ();	
  $self->{domain_names} = \@domain_names; 
  
  open (FILE,"$alignpath") or die "can't open alignment file $alignpath\n $!\n";
  while (my $line = <FILE>)
  {
  	next if $line =~ /^\s+$/;
		chomp $line;		
		if ($line =~ />(.+)/) 
		{	
			my @tmp = split " ", $1;
			my $id_num = $tmp[0];
			if (defined $id_num && length $id_num >0)
			{
				push @domain_names, $id_num;
			}
			else 
			{
				print STDERR "skippling problem fasta sequence (no id number) $line\n";
			}
		}  	
  }

  close FILE;
   
   return($self); 
  }
  


__END__

