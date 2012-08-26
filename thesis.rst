==============
Lior's Thesis
==============


Results
========

Description of the dataset
--------------------------

	The first step in training our classifier was obtaining a dataset of peptides that had their activity level with the respect to the HDAC8 enzyme verified experimentally. *Fierke et al* created a dataset composed of 361 6-mer peptides with the sequence GXK(Ac)YGC (where X,Y are all the amino acids except Cysteine). For each of these peptides, a level of activity with respect to HDAC8 was determined by measuring the percentage of deacetylation after 1 hour.(?) (**Add reference to the proper section in the supplementary material**)
	The dataset was divided to training and test sets by sorting the peptides by their activity , taking all the even rows to be the test set and all the odd rows to be the training set. That division assured even distribution of peptides with respect to their activity levels (avoiding a situation where one set holds a large number of high/low activity decoys).
	
Template selection
----------------------

	The three dimensional structure of HDAC8 was solved on numerous occasions and under different conditions in the last few years. [Add references] Below is a table that summarizes the HDAC8 structures that were tested as templates for our protocol:

======	=========	============================================================
PDB ID	Reference	Remarks
------	---------	------------------------------------------------------------
2v5w	[?]		This structure was solved as a dimer with a peptide
			substrate with 2 acetylated lysines, the peptide was taken 
			from p53 (a known substrate of HDAC8)
1vkg	[?]		First structure of HDAC ever solved, dimer with Na and Zn 
			Hydroxamate inhibitors
3ew8	[?]		Complexed as a monomer
3mz6	[?]		Fe co-factor, with D101L mutation. I mutated this 
			residue back to D, to create 3mz6_A.
======	=========	============================================================


Preparation of starting structure (+derivation of constraints)
-----------------------------------------------------------------
	#) taking the backbone from the native peptide of 2v5w, fixing the acetylated lysine , extending and mutating the residues to the desired peptide sequence.

Calibration of the protocol
------------------------------
	#) Explanation of scoring terms
	#) Different approaches for modeling (10 representatives)
		#) Naive scoring mechanism
		#) approach with different scoring 
		#) Clustering and how it influenced the results

Whole data set analysis
--------------------------
	#) measures of success
	#) determination of cutoff
	#) statistical tests

Phosphosite database
------------------------
