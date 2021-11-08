# NaPDoS2_website
This repository contains the code used to construct version 2 of the Natural Product Domain Seeker website:

	https://npdomainseeker.sdsc.edu/napdos2/
  
  It also contains file size management scripts for selecting data subsets, for example:
  
  	- Size Filtering (size_limit_seqs.pl)
	
	- File Splitting (serialize_sequences.pl)
	
	- Random Sampling (get_seq_info.pl, randomize_lines.pl, serialize_large_list.pl)
	
	- Sliding window subsequences (sequence_subdivider.pl)
		
What the NaPDoS2 website does:

	- Identifies ketosynthase (KS) and condensation (C) domains from user-submitted 
	genomic, metagenomic, and PCR amplicon data, based on sequence similarity to a 
	manually curated reference database of experimentally verified examples.
	
	- Uses a chemical phylogeny-based classification scheme to make functional
	predictions about the polyketide synthase (PKS) and non-ribosomal peptide
	synthetase (NRPS) genes in which these domains reside.
	
	- Does not require fully assembled biosynthetic gene clusters (BGCs), 
	instead relying on short KS and C domain sequence tags. This approach is particularly 
	useful for amplicons, poorly assembled genomes, or metagenomes.
	
	- The output can be used to assess the biosynthetic potential of large
	datasets and make predictions about the types of specialized metabolites that
	might be produced.
  
  Dependencies (pre-compiled binary programs called by CGI scripts)
	
	- EMBOSS, version 6.6.0, transeq module
	
	- DIAMOND, version 0.9.29
	
	- MUSCLE, version 3.8
	
	- FastTree, version 2.1.11
	
	- newick-utils, version 1.6
	
	- MySQL, version 8 (pre-loaded with classified reference database)
	
  MySQL database schema:
  
 ![NaPDoS2_schema_2021_small](https://user-images.githubusercontent.com/24737584/140801976-920fbe79-a962-4c23-af2a-347abc9313c2.png)

  Basic pipeline processing steps:
  
	1. 6-frame translation of nucleic acids to proteins(if necessary).
	
	2. Diamond BLASTP search against a curated reference database to retrieve protein domains.
	
	3. Retrieve reference database information on blast matches from a custom MySQL database.
	
	4. Tally category summary for matches in submitted data set.
	
	5. Retrieve trimmed protein sequences for selected domains via BLASTP search coordinates.
	
	6. Retrieve trimmed nucleic acid sequences via BLASTX search coordinates.
	
	7. Retrieve amino acid sequences for database neighbors identified by original BLASTP search.
	
	8. Align selected user domain protein sequences with database neighbors.
	
	9. Build tree from protein alignment.

  Script-calling order:
  
  ![napdos2_script_order_small](https://user-images.githubusercontent.com/24737584/140625623-f7516ab5-cbb5-4009-adcc-e26571246f92.png)

 
  Reference Citations:
  
	Klau LJ, Podell S, Creamer KE, Demko AM, Singh HW, Allen E, Moore BM, Ziemert N,
	Letzel AC, Jensen PR. The Natural Product Domain Seeker (NaPDoS) version 2:
	Relating ketosynthase phylogeny to biosynthetic function. (in preparation)

	Ziemert N, Podell S, Penn K, Badger JH, Allen E, Jensen PR The Natural Product
	Domain Seeker NaPDoS: A Phylogeny Based Bioinformatic Tool to Classify Secondary
	Metabolite Gene Diversity. PLoS One. 2012;7(3):e34064 Epub 2012 Mar 29. PubMed
	PMID: 22479523
    
