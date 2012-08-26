==============
Lior's Thesis
==============


Results
========

Description of the dataset
--------------------------

	The first step in training our classifier was obtaining a dataset of peptides that had their activity level with the respect to the HDAC8 enzyme verified experimentally. *Fierke et al* created a dataset composed of 361 6-mer peptides with the sequence GXK(Ac)YGC (where X,Y are all the amino acids except Cysteine). For each of these peptides, a level of activity with respect to HDAC8 was determined by measuring the precentage of deacetylation after 1 hour.(?) (**Add reference to the proper section in the supplementary material**)
	The dataset was divided to training and test sets by sorting the peptides by their activity , taking all the even rows to be the test set and all the odd rows to be the training set. That division assured even distribution of peptides with respect to their activity levels (avoiding a situation where one set holds a large number of high/low activity decoys).
	
	#) how it was prepared
	#) activity level
	
Template selection
----------------------
	#) overview on the different Zn templates

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
