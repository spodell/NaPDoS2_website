# NaPDoS2_website
This repository contains the code used to construct version 2 of the Natural Product Domain Seeker website:

	https://npdomainseeker.sdsc.edu/napdos2/
  
  It also contains file size management scripts that users can run to select subsets 
  of their data to reduce computational time and resources, for example:
  
  	- Size Filtering
	
	- File Splitting
	
	- Sequence randomization
	
  	
What the NaPDoS2 website does:

	- NaPDoS2 detects and classifies ketosynthase (KS) and condensation (C) domains
	from genomic, metagenomic, and PCR amplicon sequence data.
	
	- A phylogeny-based classification scheme is used to make functional
	predictions about the larger polyketide synthase (PKS) and non-ribosomal peptide
	synthetase (NRPS) genes in which these domains reside.
	
	- NaPDoS2 provides a rapid method to assess biosynthetic potential without
	the need for fully assembled biosynthetic gene clusters (BGCs), instead relying
	on short KS and C domain sequence tags-- this approach is particularly useful
	for amplicons, poorly assembled genomes, or metagenomes.
	
	- The output can be used to assess the biosynthetic potential of large
	datasets and make predictions about the types of specialized metabolites that
	might be produced.
  
  Basic back-end pipeline steps:
  
	1. Translate nucleic acids to proteins (if necessary).
	
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
    
