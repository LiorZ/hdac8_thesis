.. role:: ref

.. role:: label

.. raw::  latex

  \newcommand*{\docutilsroleref}{\ref}
  \newcommand*{\docutilsrolelabel}{\label}
  \newcommand*{\docutilsrolecaption}{\caption}
  
.. role:: raw-math(raw)
    :format: latex html

==============
Lior's Thesis
==============

Introduction
=============

	1) Acetylation - what it is (a reversible ptm), why is it important?
	2) Histone Deacetylases (and HAT) - description of the HDAC family, its function
	3) HDAC8 - why is it *interesting* (the specificity switch , monomer, the syndrome that was published on nature, its role in different cancers, etc)
	
The Rosetta Framework
----------------------
	
	Rosetta is a well known framework that serves as a multi-purpose toolbox in a variety of scientific studies that involve the three dimensional modeling of a macro-molecule; From design of new enzymes [citation] and symmetric proteins to predicting the structure of an RNA molecule [citation]. In its early days, Rosetta started merely as a protocol for predicting the three-dimensional structure of a protein from sequence alone, *ab-initio* modeling, a heuristic to a difficult problem which is long known to be NP-complete [9]_ . Critical to all molecular modeling problems - from design to *ab-initio* structure prediction are a reasonably accurate free-energy function and a sampling method capable of locating the minima of this function for the biomolecular system under study. 
	
	**Rosetta's scoring function** attempts to capture several hallmark features that exists in all folded structures of macro-molecules, particularly in proteins. One of these features is the nearly void-free packing of non-polar groups burying them away from water, as well as the formation of intramolecular hydrogen bonds among all buried polar atoms [10]_ . This feature is a direct consequence of the hydrophobic effect discovered by Kauzmann and was shown to be the dominant driving force in the folding of proteins [11]_ . Another feature reflects the Van-der Waals interactions between buried atoms - particularly the strong size dependence between the free energy cost of forming a cavity in the solvant to accomodate the macro molecule, and third, the free energy cost of striping water molecules from polar residues, that has to be compensated by the formation of intramolecular network of hydrogen bonds. 
	
	These features are captured in Rosetta to some extent, atom-atom interactions are computed using a Lennard-Jones potential to describe packing, a solvation model in which interactions with water molecules aren't modeled explicitly, (an implicit solvent model), to describe the hydrophobic effect and the electrostatic desolvation cost associated with burial of polar atoms, and an explicit hydrogen-bonding potential to describe hydrogen bonding. The energy function employed by Rosetta, although was proved to be robust in a plethora of studies is only a rough approximation; For start, long range electrostatic interactions that were shown to be incredibly difficult to compute because of the *induced polarization effect* , are not handled in the classic implementation of the energy function of Rosetta (Lately, a rough approximation was proven useful in a number of cases, particularly in the modeling of Protein-DNA interactions [5]_ ). Rosetta's scoring function also does not compute the entropic change that is associated with the protein attaining an ordered structure, the underlying assumption behind this omission is that entropies of different well-packed proteins are similiar.
	
	With all that said, we must note that an accurate scoring function that captures all the physical properties that are associated with protein folding and interactions is not a necessesity for the success of most variants of structural modelnig problems such as structure prediction and protein docking, rather, the success stemms from the large free-energy gap between the native structure and all the other possible conformations. 
	
	**Rosetta employs several sampling strategies** that battle the ragged energy landscape that is generally associated with macro-molecular modeling. One such powerful approach that was initially developed in *ab-initio* structure prediction is smoothing the energy landscape by modeling a low-resolution version of the interaction with a corresponding low-resolution energy function; each residue is assigned with a *centroid sphere* that encompasses its chemical properties - such hydrophobicity , polarity , etc, leading to a smoother energy landscape in which local minima are easily identified. Another important tool that aids in the location of local minima is the incorporation of a library of fragments of amino acids with defined backbones to the simulations in its early stages. The library is constructed based on sequence similarity to the query seqeunce, usually a short peptide, and on the secondary structure predicted for the peptide by PSIPRED [12]_. Fragment libraries were used extensively in our study of flexible peptide protein interactions [13]_.

Specificity prediction of peptide protein interactions
-------------------------------------------------------

	*London et al* have previously developed a general pipeline for the prediction of binding specificity of flexible peptides to protein receptors. In this pipeline, termed FlexPepBind, he modeled the structure of a collection of peptides with variable sequences to a target receptor using a high resolution peptide docking protocol - FlexPepDock [15]_ and use the energy estimation given by this protocol to each of the peptide - receptor complexes to determine their relative binding affinities and subsequently train a classifier that will be able to distinguish binders from non-binders. 
	
	This protocol has proven itself in 2 distinct biological systems - the interaction between Bcl2-like proteins and BH3 domains [7]_ which is a key feature in the regulation of apoptosis and  the farnesyltransferase (FTase) enzyme [citation] that catalyzes the attachment of farnesyl group to a protein via a thioether bond to a cysteine at or near the carboxy terminus of the protein [1,2 citation from nir's article]. *London et al* modeled the interaction between a collection of helical BH3 domains and some proteins from the Bcl-2 family and was successful in recapitulating a significant part of their specificity profile, as well as unraveling novel interactions.
	
	Unlike Bcl2-BH3, FTase is a catalytic protein that interacts primarily with *substrates*. Since FlexPepBind only models the interface between a peptide and a receptor, *London et al* assumed that binding equals catalysis and showed that this assumption is valid for the vast majority of cases. 
	
	This study is yet another adaptation of this protocol to the intriguing enzyme HDAC8 to determine its binding specificity and potentially find novel substrates. In our study we also assume that binding equals catalysis, demonstrating that this assumption is valid across a wide range of peptides. The pipeline can be summarized as follows; First, we calibrate and test our protocol for the binding of peptides that were validated experimentally by *Fierke et al*. Then, we derive a classifier and show that it indeed possesses an ability to differentiate between experimentally validated low and high affinity peptides substrates. Last, we try to find novel substrates from a large database of lysine-acetylated proteins.

Methods
========

Overview
---------
	
	We adapted FlexPepBind to predict the substrate specificity of Histone Deacetylase 8. First, we prepared a coarse starting complex of the enzyme and an array of peptides that were experimentally tested for catalytic activity, then we calibrated our protocol on a small subset of that experimentally curated dataset and obtained an initial coarse set of parameters - such as perturbation size of backbone movement and weight of certain types of features in the scoring function, this coarse set of parameters was refined by applying the pipeline on the whole training set. The performance of each set of parameters was evaluated by Pearson's correlation and on the case of the whole training set - by Kolmogorov - Smirnov goodness of fit test.


Flexible peptide - protein interactions with FlexPepDock
---------------------------------------------------------
	
	We use the previously described FlexPepBind protocol in our substrate specificity prediction of Histone Deacetylase 8. One of the most important building blocks of this protocol is a high resolution flexible peptide - protein docking protocol, FlexPepDock [15]_ . This protocol was shown to robustly refine coarse models of peptide–protein complexes into high resolution models and was later extended to model *ab-initio* peptide - protein complexes in which only the binding site and the sequence of the peptide is known. The general problem of modeling peptide - receptor interactions can roughly be divided to these subsections; 
	
	1) Model the receptor structure
	2) Predict potential binding sites on the receptor structure
	3) Model the peptide backbone on the binding site
	4) Refine the complex to higher resolution
	
	In most cases including the one we're describing in this study, the last step is sufficient - several variants of receptor structures or even closely related homologs can be obtained from the PDB database, accompanied with proteins or peptides that are already located at the binding site and provide an approximate starting structure for the refinement process [16]_ [17]_.
	
	The first step of each FlexPepDock simulation is the prepacking of the input structure to provide better packing and remove internal clashes. Side chain conformations are optimized by determining the best rotamer combination for both the protein and the peptide separately [15]_ . The second step involves 10 outer cycles of optimization. In the first cycle, the weight of the repulsive van der Waals term is reduced to 2% of its normal magnitude, and the attractive van der Waals term is increased by 225%. This allows significant perturbations within the binding pocket, while preventing the peptide and protein to separate during energy minimization. During refinement, the repulsive and attractive terms are gradually ramped back towards their original values (so that in the last cycle the energy function corresponds to the standard Rosetta score). Within each outer cycle, The rigid body orientation between the protein and the peptide is optimized, and then the peptide backbone is optimized for the new orientation. 
	
	In each such outer cycle there are 8 inner cycles in which monte carlo search with energy minimization is performed, a rigid body perturbation that is sampled from a gaussian distribution is performed and followed by sidechain repacking and minimization (The default implementation of the minimization algorithm is DFP [18]_ ) of interface residues. The metropolis criterion is applied right after the energy minimization step.
	
	Additional 8 cycles involve the optimization of the peptide backbone by applying the same, monte-carlo search with energy minimization. The backbone perturbations alternate between 2 types of moves - small and shear moves [19]_ with perturbation size of 6 degrees by default.
	
	The following figure , taken from ref [15]_ outlines the FlexPepDock protocol

.. figure:: images/fpdock.png
	:scale: 35%

	:label:`fpdock` an outline of the FlexPepDock protocol . 
	
Preparation of starting structure
---------------------------------

	For each of the peptide sequences a coarse model of the complex was generated based on the selected template, that coarse model is the starting structure that serves as input to the FlexPepDock protocol. We tested 2 approaches to create the starting complex, one involved threading the peptide sequence onto the backbone configuration taken from solved structures, the other approach included the extension of the peptide to a complete linear polypeptide (all phi angles were set to -135.0 degrees, all psi angles to +135.0 degrees) and superimposing only the acetylated Lysine onto a position taken from the crystal structure. 

	The *no free lunch* theorem suggests that all search algorithms have the same average performance over all problems [4]_, and thus implies that to gain in performance on a certain application one must use a specialized algorithm that includes some prior knowledge about that problem. In previous studies we found that incorporating key interactions between the peptide and the receptor as constraints in FlexPepDock's search algorithm greatly improves the performance of the resulting predictor. 

	Unlike previous studies, where the key interactions from which the constraints were derived relied heavily on backbone atoms, we found that the dominant interactions in our case are mostly mediated through the acetylated lysine sidechain. Furthermore, Our computational results suggests that the sidechains adjacent to the acetylated lysine form stablilizing stacking interactions with the receptor. Indeed, experimental data shows that aromatic amino acids at these positions are over represented in highly active peptides. However, we still lack a crystal structure that validates our structural hypothesis.

.. figure:: images/figure_1.png
	:scale: 20%

	:label:`keyint` The key interactions from which the constraints were derived, taken from a solved crystal complex (PDB: 2v5w).

	The interaction between D101 in the receptor and the N atom in the acetylared Lysine is critically important, a mutation D101A resulted in a complete loss of enzyme activity on the peptidic substrate and also on purified histones. [1]_ Additional constraints were derived from the interaction between the acetyl group and the two His, Asp in the active site - mostly in the purpose of fixating the acetylated Lysine in the active site.
	
Calibration of the protocol
------------------------------
	
	*London et al* [citation] developed a general framework for the prediction of binding specificity of flexible peptides to protein receptors. In general, the scheme of this framework follows a pipeline in which a collection of peptides are modeled in complex with the receptor using a high resolution peptide docking protocol [citation], then the energy estimations (termed *score*) for the modeled complexes are used to determine the relative binding affinity of each peptide to the receptor. In case the receptor is actually an enzyme that catalyzes a chemical reaction, we assume that binding = catalysis, although there are examples in which this assumption fails.[citation] 
	Previous studies have shown that a calibration process of a FlexPepBind protocol results in a more accurate predictor than a predictor that's created using a default set of parameters [citation]. The calibration process usually involves the selection of a template, adapting the scoring function and finding the right amount of sampling needed to achieve time - performance balance. [citation to bcl]

Sampling
..........
	
	The term *Sampling* in the context of FlexPepDock takes 2 different meanings. Since the entire Rosetta framework is based on non-deterministic simulation pathways, the resulting output is different from one simulation to the next and in order to capture the conformation of a complex, several simulation runs should be made so that several will eventually find the global minimal energy conformation. The other meaning of *sampling* in the context of FlexPepDock is the perturbation size of small/shear moves of the peptide backbone. A large perturbation size increases the sampling space , causing the peptide to explore more conformations.
	
	Calibrating the amount of sampling in our FlexPepBind protocol in the context of number of simulations, requires us to find the trade-off between computation time (each simulation run is computationally intensive) and number of near-native output structures (in optimal cases, the more we sample, the larger our signal/noise ratio). In the sampling space context, we aim at finding the trade-off between sampling different peptide conformations and the size of the sample space. If the peptide native structure is relatively different than the starting structure of the simulation (in term of phi/psi angles) then larger perturbations are a necessity in order to find it. Increasing the perturbation size however, can pose a probelm as it also increases the space of possible conformations, potentially decreasing the signal/noise ratio.
	
	Threading a peptide onto an existing backbone conformation in our case proved to be problematic. As we've previously mentioned, the lack of proper substrate - receptor crystal structure didn't allow us to obtain a genuine peptide - receptor complex and as a result, we couldn't reuse a reliable backbone conformation. We tried to reuse the existing peptide backbone that was present in *2v5w*, this complex was far from optimal - the peptide was located right in the interface between the two HDAC8 dimers that formed in the crystalization process, and interacted heavily with both of them. Furthermore, it contains a fluorescent coumarin residue and two acetylated lysine residues - these facts prevented the backbone conformation of this peptide from being an optimal solution, and indeed - this approach didn't yield a better predictor than the one we got when we used an extended peptide as a starting structure for our simulations.
	

.. figure:: images/2v5w_complex.png
	:scale: 25 %

	The interface between the peptide substrate that was crystallized with *2v5w*. 
	
	This backbone of this peptide was found to be a poor starting structure since it interacts with both monomers in the dimer, contains a coumarin residue (which potentially has different backbone preferences than conventional amino acids ) and two acetylated lysines.

Template selection
...................

	As we've previously discussed, our protocol models the interaction between a peptide and its corresponding receptor. FlexPepDock takes as input a three dimensional structure of the receptor and a low resolution approximation of the peptide. In our case, the receptor is HDAC8, its three dimensional structure was solved on numerous occasions and under different conditions in the last few years. In this study we tested multiple structures as templates for the FlexPepBind protocol, summarized in the table below.

.. table:: Structures of HDAC8 that were tested as templates

	======	=========	============================================================
	PDB ID	Reference	Description
	------	---------	------------------------------------------------------------
	2v5w	[1]_		HDAC8 in complex with a p53-derived diacetylated peptide 
				with a Y306F catalysis abolishing mutation
	3f07	[2]_		HDAC8 complexed with APHA
	1t67	[3]_		HDAC8 complexed with hydroxamate inhibitor (MS-344), 
				residues 62-68 were discarded from the model
	======	=========	============================================================

..

	Choosing the right template is a formidable challenge - some structures were solved with inhibitors - a thing that could induce a different *bound* structure than the actual real substrates. Others were solved with mutations that abolished catalysis and/or binding. And most of all, most structures were solved as dimers that interacted with their highly flexible regions (even though the biological active form is a monomer [1]_ ) creating crystal contacts and potential interactions that might have altered the specificity profile of the enzyme.

	In order to select a template we applied a short FlexPepDock run on each of the above recetors, complexed with the top and bottom 5 binders and used Pearson's correlation to determine how well we could distinguish between the two classes. We note that *London et al* merely used a short minimization to the template structure to select a proper template in the case of Bcl2 [7]_ , In our case, the highly flexible interface of HDAC8 indicated that a more extensive approach is needed. This short pipeline suggested that 2v5w is the best candidate for the structural template, this structure was solved together with an actual peptide, not along with a small molecule or in its free form - a fact which probably contributed to its performance as a structural template.

	In comparison, the 3f07 structure contains 3 monomers, 2 of which interact with their flexible interfaces. The ligand that interacts with the receptor is a small molecule calls APHA (aroyl pyrrolyl hydroxamate) that functions as an inhibitor. 1t67 however was solved as a monomer - a form which is identitical to the biologically active one, but some of its residues were discarded from the model and it too, was solved with an hydroxamate inhibitor.
	
.. figure:: images/interface_allReceptors.png
	:scale: 30 %

	:label:`interreceptor` **A** - The interface of 2v5w with the lysine acetylated peptide and the coumarin residue up close. **B** - An alignment of the structures from Table 1, demonstrating the conformational flexibility of the interface of HDAC8.

Scoring function
.................

	The FlexPepDock simulations were performed using both the standard Rosetta scoring schema (*score12*) and a slightly modified *score12* that includes several minor adjustments that were shown to improve the resulting classifier. The most critical change was the introduction of a weak short range Coulombic electrostatic energy term (hack_elec) In this term, a simple, linearly increasing distance-dependent dielectric was used to model solvent screening effects, with all interactions truncated at 5.5 Å, thereby preserving the short-ranged nature of the all-atom potential. *Bradly et al* [5]_ demonstrated that the incorporation of the explicit electrostatics term in addition to Rosetta's orientation-dependent hydrogen bonding potential [6]_ helped to prevent unfavorable short-range electrostatic interactions, modulated the interaction strength of charged and polar hydrogen bonds and generally, improved the performance of their DNA-protein interaction specificity predictions. This slight modification was also used by *London et al* in their Bcl-2 - BH3 specificity predictions [7]_ and in our calibration process we validated some of these parameters, verifying that they indeed introduce an improvement to the resulting predictor.
	
	We've seen in several studies conducted in our lab that a slight *post-simulation* change to the scoring function might be beneficial in determining the relative binding affinity of the peptide to the receptor. In other words, the scoring function that is used for the modeling process might be slightly different than the scoring function used to evaluate the modeled complexes after the simulation has been completed. These changes are:

	#) **Peptide score** - includes an estimation of the internal energy of the peptide
	#) **Interface score** - includes an estimation of the interactions across the interface
	#) **Reweighted score** - the sum of peptide score, interface score and total score.


	It is yet to be determined if the modification of the scoring function in the following fashion in the simulation phase itself also results in better estimation of the relative binding affinity.

Rigid body movements
.....................
	
	FlexPepDock applies rigid body movements to the peptide relative to the receptor. The transformations that define these movements are calculated using an axis and the point of center of mass of the peptide. By default , the axis equals to the vector that connects the closest peptide CA atom to the center of mass the peptide , to the closest receptor atom. Since the interaction between HDAC8 and its acetylated peptidic substrate involves a deep pocket in which the acetylated Lysine lies, we tested several alternative axes (described in figure :ref:`mc` )

.. figure:: images/anchor_arrows.png
	:scale: 30 %
	
	:label:`mc` The main axes we tested in the calibration process. One, rotating the peptide around the Lysine residue, the other approx. around the vector that is formed by the linear conformation of the peptide, X4-Ca (X - a variable position), is the default choice of the protocol.

Constraints
............
	
	HDAC8 has the ability to catalyze a deacetylation reaction with several different substrate [citation]. We believe that its ability to maintain such a diverse specificity profile stems from the fact that its binding motif is encoded in the structure of its substrates. One of our most basic assumptions when applying the FlexPepBind protocol is that the ability to characterize the structural interaction motif properly correlates the capacity to reconstruct the entire specificity profile. To this date (10/2012) there is only one solved complex containing a peptidic substrate bound to HDAC8 (PDB *2v5w*) , so finding a motif in our case was somewhat a challenge. Figure :ref:`keyint` illustrates the conserved interactions we derived from the solved complexes.
	
	Once a structural motif is determined, the scoring function must be modified to favor conformations that include that particular strucural motif. This step subsequently directs the search algorithm to sample structures that satisfy this collection of constraints. The most common types of constraints that are available in Rosetta are summarized below:
	
.. table:: Types of constraint functions in Rosetta

	=================	==========	=======================================
	Type of function	Parameters			Formula
	-----------------	----------	---------------------------------------
	Harmonic		x0, sd		.. image:: images/harmonic.png
							:scale: 50%
	Circular Harmonic	x0, sd		.. image:: images/circular_harmonic.png
							:scale: 50%
	Gaussian		mean,sd		.. image:: images/gaussian.png
							:scale: 50%
	=================	==========	=======================================

..
	
	Since we didn't want to alow much flexibility in the particular conserved interactions we defined as *conserved*, we used the harmonic function as our constraint, testing several standard deviations in our calibrations.
	
	**TODO**: add a reference to supp for the constraint file

Diffrentiation between binders and non binders
------------------------------------------------

	We used several statistical tests to evaluate the performance of our protocol and its set of parameters. The short calibration runs were evaluated by Pearson's correlation coefficient.

	While Pearson's correlation functions well on the small data set used for calibration, In larger data sets such as the training set, Pearson's correlation was shown to function poorly and doesn't provide reliable evaluation of the potential predictor's performance. In the small calibration set the of zero-activity peptides and their corresponding scores could be somewhat correlated linearly among themselves, and so does the high activity peptides. But fot the larger training set that contains peptides with all ranges of activity, this isn't necessarily the case, as the energy estimations given to each of the peptides by our protocol aren't necessarily in a *linear* correlation with the level of activity. For the purpose of evaluating our ability to differentiate between binders and non binders in the whole training set we used the Kolmogorov Smirnov goodness-of-fit test. This test quantifies a distance between the empirical distributions of two samples - in our case - binders and non-binders. The resulting p-value is calculated under the null hypothesis that the samples are drawn from the same distribution.

Results
========


Description of the dataset
--------------------------

	*Fierke et al* created a dataset composed of 361 6-mer peptides with the sequence GXK(Ac)YGC (where X,Y are all the amino acids except Cysteine). For each of these peptides, a level of activity with respect to HDAC8 was determined by measuring the percentage of deacetylation after 1 hour.(?) (**Add reference to the proper section in the supplementary material**)
	We divided the dataset to training and test sets by sorting the peptides according to their experimental activity , taking all the even rows to be the test set and all the odd rows to be the training set. That division assured even distribution of peptides with respect to their activity levels (avoiding a situation where one set holds a large number of high/low activity decoys).
	

Calibration of the protocol
------------------------------
	
	Below we describe the results obtained in the calibration process. This process resulted in a coarse set of parameters, to be refined on the whole training set as part of the classifier learning process. Usually, Each step of the calibration process involved changing one degree of freedom of a certain feature (such as - amount of sampling, constraints, etc) while maintaining the others fixed.
	The performance of each simulation was evaluated by the Pearson correlation coefficient by averaging the score of the top 3 models with the lowest peptide , interface and reweighted score against. The tables that summarize the performance of each of these simulations can be found in the  `Calibration simulations and their performance` section, in the `Supplementary Material`_. Plots that show the distribution of score of each sequence against its experimental activity are available in section `Calibration`_ in the `Supplementary Material`_.
	
	The first calibration round was made by taking 5 best binders and 5 bad binders, trying to generate a coarse set of parameters to be refined later using the entire training set.

.. table:: A short version of the dataset used for coarse calibration of our protocol.
	
	+---------------+----------------------+------------------+
	|Sequence	|	% deacetylation|	annotation|
	+===============+======================+==================+
	|GYK(ac)FGC	|93		       |		  |
	+---------------+----------------------+		  |
	|GYK(ac)WGC	|80		       |		  |
	+---------------+----------------------+     Binders	  |
	|GLK(ac)FGC	|66		       |		  |
	+---------------+----------------------+		  |
	|GIK(ac)FGC	|64		       |		  |
	+---------------+----------------------+		  |
	|GRK(ac)YGC	|62		       |		  |
	+---------------+----------------------+------------------+
	|GQK(ac)YGC	|0		       |		  |
	+---------------+----------------------+		  |
	|GIK(ac)VGC	|0		       |		  |
	+---------------+----------------------+   Non Binders	  |
	|GMK(ac)VGC	|0		       |		  |
	+---------------+----------------------+		  |
	|GDK(ac)YGC	|0		       |		  |
	+---------------+----------------------+		  |
	|GMK(ac)YGC	|0		       |		  |
	+---------------+----------------------+------------------+
..


	This set of short simulations allowed us to quickly distinguish between sets of parameters.
	

Sampling
.........

	We inspected different amounts of sampling in which the number of decoys generated and the amount of perturbation size were modified together, since As we've previously mentioned, the larger the perturbation size - the larger the space of possible peptide conformations.
	
	Since the amount of sampling was the first feature we decided to calibrate, we initialized the other features with values that were found optimal in previous studies such as:
	
	#) Weight of *hackelec* (electrostatic term) = 0.5
	#) Standard deviation of constraints = 0.2
	#) Number of decoys generated per simulation = 200

	These features were of course, validated and perturbed in later phases.
	
	We also figured that the default anchor chosen in the FlexPepDock protocol will not be optimal in our case, so we started with a predefined anchor that we found to be suitable, and verified its optimality later on when other sets of parameters were calibrated. Furthermore, since it is unlikely that the amount of sampling will be different from one template to another, we selected 2v5w since it is the one that has the best chances to serve as a template, due to the properties we mentioned earlier (primarily since it was solved with an actual peptide and not a small molecule)

.. table:: Calibration of the amount of sampling.

	+---------------+--------------------------------+----------------------------------------------------+
	|		|	 **Sampling**        	 |       **Scoring scheme** (correlation coefficient) |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|No.		|Perturbation size |  No. decoys | Peptide score | Interface score | Reweighted score |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|1		|30		   |  200	 | -0.45	 | -0.69	   | -0.32	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|2		|60		   |  500	 | -0.38	 | -0.65	   | -0.26	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|3		|90		   |  900	 | -0.27	 | -0.58	   | 0.48	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|4		|30		   |  500	 | -0.46	 | -0.75	   | -0.21	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|5		|20		   |  200	 | -0.464	 | -0.76	   | -0.24	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|8		|6 (default value) |  200	 | -0.24	 | -0.72	   | -0.121	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|9		|15		   |  200	 | -0.41	 | -0.77	   | -0.24	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+
	|16		|15		   |		 |		 |		   |		      |
	|		|low resolution    |  		 |		 | 		   |		      |	
	|		|pre-optimization  |		 |		 |		   |		      |
	|		|(centroid mode)   |  200	 | -0.41	 | -0.77    	   | -0.24	      |
	+---------------+------------------+-------------+---------------+-----------------+------------------+


..


	Our findings above suggests that a modest amount of sampling (in the context of number of simulation runs) is sufficient to generate a reliable predictor. Our findings correlate with an earlier study conducted by *London et al* [8]_ , that found that 200 simulation rounds are indeed sufficient for this purpose, and that a larger number of simulation rounds doesn't necessarily yield significant improvements in the perdictor's performance. However, in terms of the perturbation size, we found that the default amount of sampling in FlexPepDock (simulation number 8) that was sufficient for all previous studies, wasn't optimal in our case, perhaps since our simulation started from an extended peptide conformation, while all other studies reused an existing backbone conformation as a template that all the sequences were threaded on. Furthermore, this short set of calibration runs suggests that the interface scoring scheme functions better than the rest in the task of diffrentiating between binders and non binders.
	
Template selection
...................

	We applied a short FlexPepDock run on each of the possible templates complexed with the top and bottom 5 binders and used Pearson's correlation to determine how well we could distinguish between the two classes. 
	
	+----------------------------------+----------------------------------------------------+
	|			 	   |       **Scoring scheme** (correlation coefficient) |
	+---------------+------------------+---------------+-----------------+------------------+
	|No.		|Template	   | Peptide score | Interface score | Reweighted score |
	+---------------+------------------+---------------+-----------------+------------------+
	|9		|2v5w		   | -0.41	   | -0.77	     | -0.24   		|
	+---------------+------------------+---------------+-----------------+------------------+
	|13		|3f07		   | 0.44	   | -0.51	     | -0.51   		|
	+---------------+------------------+---------------+-----------------+------------------+
	|15		|1t67		   | -0.11	   | -0.11	     | -0.6   		|
	+---------------+------------------+---------------+-----------------+------------------+	

	These short simulations validate our initial assumption that *2v5w* is the best candidate for a template. 
	
Scoring function
.................

	In our calibration of the scoring function we were interested to see whether our initial parameters - the use of the short electrostatic term (hackelec) and the lazaridis karplus modification should be refined or modified. For that, we tried to use Rosetta's default scoring function *score12* and decreased the weight of hackelec in the scoring function.
	
	+----------------------------------------------+----------------------------------------------------+
	|		                	       |       **Scoring scheme** (correlation coefficient) |
	+---------------+------------------------------+---------------+-----------------+------------------+
	|No.		|Scoring function  	       | Peptide score | Interface score | Reweighted score |
	+---------------+------------------------------+---------------+-----------------+------------------+
	|9		|weight of hackelec = 0.5      | -0.41         | -0.77	         | -0.24   	    |
	+---------------+------------------------------+---------------+-----------------+------------------+	
	|10		|weight of hackelec = 0.25     | -0.45         | -0.56	         | -0.31   	    |
	+---------------+------------------------------+---------------+-----------------+------------------+
	|7		|*score12* (hackelec=0)        | -0.48         | -0.7	         | -0.28   	    |
	+---------------+------------------------------+---------------+-----------------+------------------+
	
	Looking at the results, clearly, our initial assumption looks valid - the correlation coefficient is optimal in simulation 9 where the weight of hackelec is 0.5. 
	
Rigid body movements
.....................
	
	We've tested several approaches to the way we perform rigid body movements. As we've previously mentioned, the axis that determines the transformations of the peptide relative to the receptor equals to the vector that connects the closest peptide CA atom to the center of mass the peptide , to the closest receptor atom. We've tried to cleaverly select these two atoms so that different axes will be used by the protocol , so that consequently, different axes will be used for the rigid body transformations.
	
	+--------------------------------------------------------+----------------------------------------------------+
	|		                		         |       **Scoring scheme** (correlation coefficient) |
	+---------------+----------------------------------------+---------------+-----------------+------------------+
	|No.		|Anchor (residue) 	  	         | Peptide score | Interface score | Reweighted score |
	+---------------+----------------------------------------+---------------+-----------------+------------------+
	|9		| 366 (CA atom)		                 | -0.41         | -0.77	   | -0.24            |
	+---------------+----------------------------------------+---------------+-----------------+------------------+
	|6		| 367 (chosen automatically -		 |		 |		   |		      | 
	|		| center of mass of the peptide)         | -0.49         | -0.65	   | -0.51            |
	+---------------+----------------------------------------+---------------+-----------------+------------------+
	|12		| 366 (anchor atom was CH, instead of CA)| -0.45         | -0.77	   | -0.41            |
	+---------------+----------------------------------------+---------------+-----------------+------------------+
	|17		| 366 , receptor anchor was 		 |		 |		   |		      |
	|		| the CA atom of residue no. 289	 | -0.48	 | -0.74	   | -0.38            |
	+---------------+----------------------------------------+---------------+-----------------+------------------+		
	
	Looking at the results we see that either of the atoms in residue 366 can be selected as anchors, yielding similar ability to distinguish between binders and non binders.
	
	TODO: Insert a figure of all the axes.
	
Constraints
............

	We tested few different values for the standard deviations of the constraints that were introduced to the simulations. (see figure `keyint`) We note that a simulation with no constraints at all generated model structures in which the peptide didn't bind the active site at all and thus, weren't relevant for comparison.
	
	+------------------------------------------------+----------------------------------------------------+
	|		                		 |       **Scoring scheme** (correlation coefficient) |
	+---------------+--------------------------------+---------------+-----------------+------------------+
	|No.		|Constraints (standard deviation)| Peptide score | Interface score | Reweighted score |
	+---------------+--------------------------------+---------------+-----------------+------------------+
	|9		| 0.2 Å 	                 | -0.41         | -0.77	   | -0.24            |
	+---------------+--------------------------------+---------------+-----------------+------------------+
	|18		| 0.15 Å 	                 | -0.45         | -0.54	   | -0.38            |
	+---------------+--------------------------------+---------------+-----------------+------------------+
	|19		| 0.25 Å 	                 | -0.47         | -0.51	   | -0.28            |
	+---------------+--------------------------------+---------------+-----------------+------------------+

	Surprisingly, a slight modification to the standard deviation of the constraints yields drastic change in our ability to distinguish binders from non binders.
	
Threading the peptide
......................
	
	In the Methods section we've discussed the reasons that led us to use primarily extended conformations as the starting structure for the peptide. We verified this hypothesis in a simulation that incorporated the threading of peptides onto the existing starting structure from *2v5w* with a parameter-set that's identical to simulation 9 that achieved the best performance in terms of Pearson's correlation coefficient:
	
	* Pearson's Correlation coefficient for the following scoring schemes:
		* Interface score: -0.784
		* Peptide score: -0.64
		* Reweighted score: -0.003
		
	Comparing to simulation #9 and its set of parameters and in contrast to our initial assumption, this simulation achieved the best correlation with experimental data. 
	
Summary of calibration runs
............................
	
	This phase of calibration allowed us to select an initial set of parameters lately to be refined on the whole training set. With this calibration approach we could easily discard sets of parameters that failed to identify highly reactive substrates, and falsly identified zero activity substrates. We note simulation #11 and simulations #9 and its set of parameters, using the interface scoring scheme yielded the best performance in terms of Pearson's correlation coefficient. We also noticed that the interface scoring scheme achieved superior performance to the rest of the schemes for every parameter set we've tested. Moreover, the reweighted score scheme that demonstrated good ability to distinguish binders from non binders in previous studies, failed in our case.
	
	In the next phase , in which we run our peptide modeling protocol on the whole training set, we mainly use the set of parameters that exhibited superior performance in the short calibration phase.

Whole data set analysis
--------------------------
	
Training a classifier
.....................

	After an initial phase of calibration , we were set to examine the parameters learned from the brief simulations on the whole training set, this step allowed us to refine our initial, coarse set of parameters. Below is a table that summarizes the simulations we've performed on the whole training set.

	For each of these simulations and for each scoring scheme we calculated the Pearson's correlation coefficient to evaluate its fitness to experimental data. 
	Let us remember that our dataset contains sequences of lysine acetylated peptides that are ranked by their level activity as substrates. The peptide's level of activity is not represented in a binary fashion (binder / non-binder) , but rather as a continous value in [0,1]. In order to train a binary classifier, we needed to adapt our dataset accordingly, to a binary representation. To accomplish that, we selected an experimental level of activity to serve as a cutoff so that each sequence with activity that is lower than the cutoff is labeled as a non-binder and vice versa. We derived that cutoff by applying 2 samples KS test on all possible activity levels ([0,1], in resolution of 0.01), the activity level that was chosen as cutoff is the one that obtained the lowest p-value in the KS test, thus, the one that could best differentiate between the 2 distributions of *scores* - that of the binders and the score distribution of non binders.  (see figure :ref:`cutoff`)
	
.. figure:: plots/cutoff.png
	:scale: 50 %

	:label:`cutoff` log(p-value) of KS test when using the cutoff from the X axis (simulation 1). Clearly, the best cutoff we can choose in this case is 0.34.

..

TODO: Replace that figure with one that doesn't have a red underline in the word deacetylation.


	This table summarizes the simulations we performed on the whole training set, each of the columns describe a different aspect of the parameter set used.
	
	
.. table:: Summary of training set simulations

	======		================	===============================	===========	===================
	No.		Anchor (residue)	Sampling			Template	Scoring function
	======		================	===============================	===========	===================
	1		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 simulations per peptide.			* hack_elec = 0.5

	2		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 simulations per peptide.	(threaded)	* hack_elec = 0.5	

	3		366			* perturbation size = 15	3f07		* Lazaridis-Karplus
						* 200 simulations per peptide.			* hack_elec = 0.5

												  
	4		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
			anchor was CH		* 200 simulations per peptide.			* hack_elec = 0.5
			atom									

	5		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
			anchor was CH		* 200 simulations per peptide.			* hack_elec = 0.5
			atom			* low resoultion preopt.							

	6		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 simulations per peptide.			* hack_elec = 0.5
												* sd of constraints
												  is 0.15

	7		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 simulations per peptide.			* hack_elec = 0.5
												* sd of constraints
												  is 0.25
	======		================	===============================	===========	===================

..

	
	Simulations 6 and 7 achieved the best KS p-values on the training set, 1.51×10\ :sup:`-5` and 2.79×10\ :sup:`-5` respectively, using the peptide scoring scheme. However the cutoff that's responsible for these low p-values is 0.44 which is relatively high and isn't sensitive enough (there are only 11 out of 181 peptides with higher activity levels). Simulation #4 showed a potentially good ability to differentiate between binders and non-binders with cutoff of 0.35 and KS p-value of 4.63×10\ :sup:`-5`. 
	
	We applied a clustering step [citation] to the structures from each simulation and averaged the top 3 ranking decoys in the largest cluster to get a score for each peptide. In cotrast to previous findings in earlier studies [7]_ , [8]_, we found that clusering improves the differentiation between binders and non binders in several orders of magnitude. For example, Simulation #4, the one with the lower standard deviation for the constraints demonstrated the best performance with the interface scoring scheme and a KS p-value of 4.89×10\ :sup:`-7` which is two orders of magnitudes increment from the lowest p-values that we obtained without clustering. Another notable candidate was Simulation #2 , in this simulation we threaded the peptide onto the existing backbone conformation, using the peptide scoring scheme it showed a p-value of 4.03×10\ :sup:`-6` using a cutoff of 0 activity level. This parameter set indeed demonstrate both specificity and a very high sensitivity in differentiating between binders and non-binders.
	
	Interestingly, we saw the level of activity of 0.34 reccur as a cutoff for a number of well performing parameter sets that achieved low p-values after clustering under different scoring schemes. For example , simulation #1 that has the parameter set that was one of the best performing in the first initial calibration phase with the interface scoring scheme achieved a p-value of 4.4×10\ :sup:`-6` - three orders of magnitudes improvement comparing to its performance without clustering.

	The `Training set simulations and their performance`_ concentrates a summary of all simulations with and without a clustering step, including the statistical evaluation of their performance. 

	To visualize the comparison of our ability to distinguish binders from non binders with and without clustering, we plotted *score vs. activity* plots for all simulations. They are available in the `Supplementary Material`_ - `Training set analysis`_
	From the results above we were able to derive a modeling scheme that could serve us in our future predictions for additional substrates - the scheme we used in simulation #1 together with a clustering step achieved best AUC together with the 0.34 cutoff we obtained. (see figure :ref:`roc`)
	
Comparison to a minimization only based classifier
...................................................

	Previous studies have indicated that a minimization only scheme could yield suprisingly good predictors and as a result, posses a ability to distinguish binders and non binders in several biological systems [7]_ [8]_. The FlexPepDock protocol applies a minimization scheme in which only the corresponding peptide and the interface residues are minimized while the whole receptor structure stays fixed. We've applied several different minimization schemes to our training set to evaluate and compare the ability of both methods - the full optimization that uses the FlexPepDock modeling protocol and the a simple minimzation of the interface and peptide employed by FlexPepDock. We've tried several approaches:
	
	1) Minimization with *score12*, rest is similar to Simulation #1 applied to the whole training set
	2) Minimization with the same modification to the scoring function as Simulation #1 (hackelec, Lazaridis-Karplus) applied to the whole training set
	3) Minimization starting from threaded peptides, identical to simulation #2 applied to the whole training set
	
	Surprisingly , the 1st approach - the one that didn't require any changes to the scoring function was the one that best correlated with experimental data and showed the best ability so far to distinguish binders from non binders with a KS p-value of 5.95×10\ :sup:`-10` and a cutoff of 0.34 using the peptide scoring scheme - three orders of magnitude improvement to full optimization simulations. The 2nd approach also performed well with a KS p-value of 4.6×10\ :sup:`-8` and a cutoff of 0.34, using the peptide scoring scheme. The 3rd approach failed to improve any of the p-values obtained in the full simulation runs. Figure :ref:`roc` shows an ROC plot comparing the performance of possible predictors derived from both types of best performing simulations - minimization only and full optimization.

Test set analysis
..................

	With our insights from training a classifier on the training set, we applied it on the other part of the sequences - the test set. The simulation scheme used the set of parameters and constraints identical to that of simulation #1 in the training set runs, as its resulting predictor has the best ability to distinguish between binders and non binders (ROC plot AUC of 0.873).
	The below ROC plot summarizes the performance of our classifier on the test set, comparing to its performance on the training set and to a minimization only scheme.


	.. figure:: plots/ROCPlots/roc.png
		:scale: 50 %

		:label:`roc` Comparison of the minimization and full optimization schemes that included clustering on both training and test sets.
	
		The minimization step uses the *peptide scoring scheme*, while in the full optimization the inteface scoring scheme performed better on the training set and thus - served as the basis for the predictor on the test set.

Searching for novel, non-histone substrates
--------------------------------------------

	We used the minimization only version of our predictor - the one that performed best on the experimental dataset - to search for potential novel substrates of HDAC8.
	We've obtained a copy of the Phosphosite database from PhosphoSitePlus (PSP) - an online systems biology resource providing comprehensive information and tools for the study of protein post-translational modifications and queried it for lysine acetylated proteins. We've trimmed the sequences so they will be of the same size as the sequences that are present in the experimental dataset - **YYK(ac)YYY**. 

	To demonstrate the ability of our classifier to recognize potential substrates among the large database of acetylated sequences we plotted the distribution of scores of all the acetylated sequences from the database against a background distribution of random peptides that were sampled from the distribution of amino acids in the acetylated sequences (figure :ref:`phosphodist`) and under the null hypothesis that both sequences were originated from the same distribution, we used the Kolmogorov-Smirnov test to calculate a p-value of 5.07×10\ :sup:`-5`.
	It is important to note that surely, not all sequences in the Phosphosite database are substrates of HDAC8, but nevertheless, we were managed to diffrentiate between a collection of random sequences and a collection of acetylated sequences that some of them were putatively originated from potential substrates of HDAC8. This finding could suggest that there are quite a lot potential substrates of HDAC8 that are yet to be discovered.

	.. figure:: plots/PhosphositeDisr/plot.png
		:scale: 50 %

		:label:`phosphodist` Distribution of scores in both acetylated and random sequences
	
		The rightmost bar concentrates all the peptides that have a minimization score above 10. (a high score that suggests that these peptides were not modeled successfully)

HDAC8 and CdLS syndrome
........................
	
	A recent study [23]_ nominated the loss of function of HDAC8 as one of the causes to the Cornelia de Lange syndrome (CdLS) that occurs due to a malfunction in the cohesin acetylation cycle. In humans the cohesin is a multisubunit complex that is made up of SMC1A, SMC3, RAD21 and a STAG protein. These form a ring structure that is proposed to encircle sister chromatids to mediate sister chromatids cohesion [20]_ and also has key roles in gene regulation [21]_ . Using a monoclonal antibody specific for acetylated SMC3 the researchers found that the total levels of SMC3 is constant throughtout the cell cycle while SMC3-ac levels rapidly decline during mitosis, a finding that suggested a coordinated deacetylation. The researchers therefore used RNAi for each of the known histone deacetylases and sirtuins and identified HDAC8 as the primary SMC3 deacetylase. Indeed, SMC3 has 6 known acetylation sites [22]_ , 3 of them obtained low scores indicating them as HDAC8 deacetylation sites by our protocol: 
	
.. table:: SMC3 known acetylation sites with FlexPepBind scores
	
	=================	============	============
	Position
	of Deacetylation	Sequence	FPBind score
	-----------------	------------	------------
	106			AKK(ac)DQY 	672.779
	1190			GVK(ac)FRN 	125.366
	336			LEK(ac)IEE 	25.855
	215			YQK(ac)WDK 	-2.082
	105			GAK(ac)KDQ 	-4.027
	140			IVK(ac)QGK 	-6.222
	=================	============	============

..

	
	**Are there any more deactylation sites?** We were interested to see whether our protocol can capture additional deacetylation sites that aren't known yet. For that, we trimmed the SMC3 sequence to short peptides , 6 residues, wherever there was a lysine ( in format identical to the YYK(ac)YYY format, see Figure :ref:`smc3seq`).
	
.. figure:: images/peptide_collection_arrows.png
	:scale: 55%

	:label:`smc3seq` From each possible acetylation site (each lysine in SMC3 sequence) we created a peptide as input to our protocol to find putative deacetylation sites

..

	Results from the minimization version of our protocol that achieved superior results in earlier tests indicate that there are 13 additional possible deacetylation sites, assuming these sites undergo acetylation in the first place. see table in *HDAC8 and CdLS syndrome* in the supplementary material.
	
	Mutant SMC1A proteins account for ~ 5% of the cases of CdLS and is shown to have several mutations in a number of patients and number of sites [24]_. We tested whether any of these mutations are known acetylation sites and whether these acetylation sites are recognized by our protocol as HDAC8 deacetylation positions.
	
.. figure:: images/SMC1A_mutations.png
	:scale: 40%

	:label:`smc1amut` Known acetylation sites and observed mutations in SMC1A, see summary on the table below
	
	**A** - SMC1A sequence annotated with known acetylation sites and mutations, as well as peptides trimmed from the protein that predicted to bind when tested as potential acetylated peptides. (peptides > 6 residues indicate overlapping) **B** - Reproduced from [24]_ , A schema of SMC1A structure annotated with mutations that were discovered in different patients
	

.. table:: Lysine acetylation positions

	+--------+
	|Position|
	+--------+
	|282	 |
	+--------+
	|437	 |
	+--------+
	|536	 |
	+--------+	
	|648	 |
	+--------+	
	|713	 |
	+--------+
	
..
	
	
.. table:: Mutations that were observed in different patients in the SMC1A protein

	=========	==================
	Position	Mutation Type
	---------	------------------
	58-62		deletion: V58-R62
	133		F133V
	196		R196H
	493		E493A
	496		R496C, R496H
	711		R711W
	790		R790Q
	832		D831_Q832delinsE
	1122		R1122L
	=========	==================
	
..

	
	Worth noting is the mutation **R711W** that is located right close to a known acetylation site in the coiled coil region and was predicted by our classifier as a binder. A mutated version of the peptide - **WLKYSQ** was predicted as a non-binder. The authors of the study in ref [24]_ used the Coils program [25]_ , that predicts the probability of protein to form a coiled coil and concluded that the R711W mutation has a low likelihood of disrupting the coiled coil. However, the authors speculate that the alterations caused by this mutation may affect the angulation of the coiled-coil resulting in impaired intra or intermolecular approximation of the SMC head domains, or disrupt binding of accessory proteins to the cohesin ring. Our findings suggests yet another possibility - the R711W mutation might disrupt the acetylation or deacetylation of SMC1A at position 713, and that might contribute to the protein inability to bind accessory proteins or attain a non-functioning structure.
	
	Position K437 is also a known acetylation site according to ref [22]_ and the peptide **IEKLEE**  that overlaps this position is predicted by our protocol to undergo deacetylation by HDAC8. 
	
	
Summary
--------

	We have previously used structure-based prediction of binding specificity to successfully identify both known and novel protein farnesyltransferase (FTase) substrate peptides and BH3 peptides to Bcl-2-like proteins. The HDAC8 system presents additional challenges to systems we studied previously - the extremely flexible loops in the interface has the ability to move and accomodate different substrates for each conformation, the lack of solved crystals that incorporated a genuine substrate and the acetylated lysine - a post translational modification that was poorly addressed in previous computational studies.
	In this study, We've applied the FlexPepBind modeling scheme to a series of peptide sequences in order to train a predictor that will have the ability to distinguish between peptides that serve as substrates of HDAC8 and peptides that are doesn't. Since FlexPepDock only models the interface between the two , and not the catalytic process, we've assumed that peptides that bind the receptor are necessarily deacetylated and going through the whole catalytic process. 

	We calibrated a set of parameters that included the amount of sampling and movement, degree of constraints and some other energy terms in the scoring function and compared the resulting predictor to a predictor that was obtained by applying much simpler and less computationally intensive approach - the FlexPepDock minimization scheme. The minimization only predictor performed better in the task of separating between binders and non binders in the experimental dataset we used. Its ability, in addition to the fact that this scheme is much less computationally intensive, lead us to utilize it to find new potential substrates to HDAC8 in a large database of acetylated proteins.

Supplementary Material
=======================

Calibration
------------

Calibration simulations and their performance
.............................................

Summary of calibration runs
````````````````````````````

.. table:: Description and summary of calibration simulations.

	======		================	===============================	===========	===================
	No.		Anchor (residue)	Sampling			Template	Scoring function
	------		----------------	-------------------------------	-----------	-------------------
	1		366			* perturbation size = 30	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5
	
	2		366			* perturbation size = 60	2v5w		* Lazaridis-Karplus
						* 500 decoys per peptide.			* hack_elec = 0.5
						
	3		366			* perturbation size = 90	2v5w		* Lazaridis-Karplus
						* 900 decoys per peptide.			* hack_elec = 0.5

	4		366			* perturbation size = 30	2v5w		* Lazaridis-Karplus
						* 500 decoys per peptide.			* hack_elec = 0.5
	
	5		366			* perturbation size = 20	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5

	6		367 (chosen		* perturbation size = 20	2v5w		* Lazaridis-Karplus
			automatically		* 200 decoys per peptide.			* hack_elec = 0.5
			since its the 
			center of mass)	
			
	7		366			* perturbation size = 20	2v5w		* Rosetta's default
						* 200 decoys per peptide.			  score function
												  (score12)
	8		366			* perturbation size = 6 
						  (default)			2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5

	9		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5

	10		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.25
	
	11		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.	(threaded)	* hack_elec = 0.5
										[*]_	
														
	12		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
			(anchor was CH		* 200 decoys per peptide.			* hack_elec = 0.5
			atom, instead of
			CA)	
	
	13		366			* perturbation size = 15	3f07		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5
	
	14		366			* perturbation size = 15	3f07		* Lazaridis-Karplus
			(anchor was CH		* 200 decoys per peptide.			* hack_elec = 0.5
			atom instead of
			CA)								
	
	15		366			* perturbation size = 15	1t67		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5

	16		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5
						* low resolution step 
						  (centroid mode)						
	
	17		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
			receptor anchor		* 200 decoys per peptide.			* hack_elec = 0.5
			was 289 
			(manually)
			[*]_
	
	18		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5
												* sd of constraints
												  is 0.15
												  
	19		366			* perturbation size = 15	2v5w		* Lazaridis-Karplus
						* 200 decoys per peptide.			* hack_elec = 0.5
												* sd of constraints
												  is 0.25		
	======		================	===============================	===========	===================
	
..

.. [*] The sequence was threaded on the peptidic substrate backbone in the 2v5w crystal. Since this peptidic substrate was only 4 amino acid long (the train/test sequences were 6 residues long), the 2 extra amino acids backbone conformation attained an extended conformation.

.. [*] Setting the receptor anchor to be the 289 residue , creating an axis that aligns with the Lysine residue side-chain. This axis is directed inside the pocket , and allowed the peptide to rotate while the Lysine residue stays fixed (see figure :ref:`mc`)

Peptide Score
``````````````

.. table:: Results for short calibration runs, by peptide score.

	=====	==========================================
	No.	Pearson correlation coefficient
	-----	------------------------------------------
	1	* R: -0.45
		* p-Value: 0.18
		
	2	* R: -0.38
		* p-Value: 0.27

	3	* R: -0.27
		* p-Value: 0.44

	4	* R: -0.46
		* p-Value: 0.18

	5	* R: -0.464
		* p-Value: 0.176
		
	6	* R: -0.493
		* p-Value: 0.146
		
	7	* R: -0.48
		* p-Value: 0.152
		
	8	* R: -0.24
		* p-Value: 0.498
		
	9	* R: -0.41
		* p-Value: 0.230

	10	* R: -0.45
		* p-Value: 0.185

	11	* R: -0.64
		* p-Value: 0.043
		
	12	* R: -0.45
		* p-Value: 0.202
		
	13	* R: 0.44
		* p-Value: 0.185

	14	* R: 0.79
		* p-Value: 0.006
		
	15	* R: -0.11
		* p-Value: 0.75
		
	16	* R: -0.3
		* p-Value: 0.39
		
	17	* R: -0.48
		* p-Value: 0.153
		
	18	* R: -0.45
		* p-value: 0.15

	19	* R: -0.47
		* p-value: 0.16

	=====	==========================================


Interface Score
`````````````````

.. table:: Results for short calibration runs, by interface score.

	=====	==========================================
	No.	Pearson correlation coefficient
	-----	------------------------------------------
	1	* R: -0.69
		* p-Value: 0.02
		
	2	* R: -0.65
		* p-Value: 0.04

	3	* R: -0.58
		* p-Value: 0.07

	4	* R: -0.75
		* p-Value: 0.012

	5	* R: -0.76
		* p-Value: 0.01
		
	6	* R: -0.65
		* p-Value: 0.04
		
	7	* R: -0.7
		* p-Value: 0.02
		
	8	* R: -0.72
		* p-Value: 0.018
		
	9	* R: -0.77
		* p-Value: 0.008

	10	* R: -0.56
		* p-Value: 0.085

	11	* R: -0.784
		* p-Value: 0.007
		
	12	* R: -0.77
		* p-Value: 0.009
		
	13	* R: -0.51
		* p-Value: 0.130

	14	* R: -0.174
		* p-Value: 0.62
		
	15	* R: -0.11
		* p-Value: 0.75
		
	16	* R: -0.542
		* p-Value: 0.1
		
	17	* R: -0.74
		* p-Value: 0.013
		
	18	* R: -0.54
		* p-Value: 0.1

	19	* R: -0.51
		* p-value: 0.13
	=====	==========================================


Reweighted Score
`````````````````

.. table:: Results for short calibration runs, by reweighted score.

	=====	==========================================
	No.	Pearson correlation coefficient
	-----	------------------------------------------
	1	* R: -0.32
		* p-Value: 0.35
		
	2	* R: -0.26
		* p-Value: 0.46

	3	* R: 0.48
		* p-Value: 0.156

	4	* R: -0.21
		* p-Value: 0.54

	5	* R: -0.24
		* p-Value: 0.49
		
	6	* R: -0.51
		* p-Value: 0.13
		
	7	* R: -0.28
		* p-Value: 0.42
		
	8	* R: -0.121
		* p-Value: 0.738
		
	9	* R: -0.24
		* p-Value: 0.496

	10	* R: -0.31
		* p-Value: 0.382

	11	* R: -0.003
		* p-Value: 0.99
		
	12	* R: -0.41
		* p-Value: 0.23
		
	13	* R: -0.51
		* p-Value: 0.130

	14	* R: -0.6
		* p-Value: 0.06
		
	15	* R: -0.19
		* p-Value: 0.59
		
	16	* R: -0.008
		* p-Value: 0.98
		
	17	* R: -0.38
		* p-Value: 0.27
		
	18	* R: -0.28
		* p-value: 0.08

	19	* R: -0.09
		* p-value: 0.2
	=====	==========================================

Score vs. Activity plots
.........................
.. list-table:: Training set - score vs. activity plots for the short calibration phase
   :widths: 5 30 30 30
   :header-rows: 1

   * - No.
     - Reweighted Score
     - Peptide Score
     - Interface Score
   * - 1
     - .. image:: plots/ShortCalibration/calibration2_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration2_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration2_I_sc_activity_score.png
     	:scale: 20%
   * - 2
     - .. image:: plots/ShortCalibration/calibration3_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration3_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration3_I_sc_activity_score.png
     	:scale: 20%
   * - 3
     - .. image:: plots/ShortCalibration/calibration4_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration4_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration4_I_sc_activity_score.png
     	:scale: 20%
   * - 4
     - .. image:: plots/ShortCalibration/calibration5_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration5_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration5_I_sc_activity_score.png
     	:scale: 20%
   * - 5
     - .. image:: plots/ShortCalibration/calibration6_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration6_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration6_I_sc_activity_score.png
     	:scale: 20%
   * - 6
     - .. image:: plots/ShortCalibration/calibration7_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration7_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration7_I_sc_activity_score.png
     	:scale: 20%
   * - 7
     - .. image:: plots/ShortCalibration/calibration8_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration8_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration8_I_sc_activity_score.png
     	:scale: 20%
   * - 8
     - .. image:: plots/ShortCalibration/calibration9_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration9_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration9_I_sc_activity_score.png
     	:scale: 20%
   * - 9
     - .. image:: plots/ShortCalibration/calibration10_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration10_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration10_I_sc_activity_score.png
     	:scale: 20%
   * - 10
     - .. image:: plots/ShortCalibration/calibration12_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration12_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration12_I_sc_activity_score.png
     	:scale: 20%
   * - 11
     - .. image:: plots/ShortCalibration/calibration13_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration13_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration13_I_sc_activity_score.png
     	:scale: 20%
   * - 12
     - .. image:: plots/ShortCalibration/calibration14_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration14_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration14_I_sc_activity_score.png
     	:scale: 20%
   * - 13
     - .. image:: plots/ShortCalibration/calibration33_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration33_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration33_I_sc_activity_score.png
     	:scale: 20%
   * - 14
     - .. image:: plots/ShortCalibration/calibration32_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration32_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration32_I_sc_activity_score.png
     	:scale: 20%
   * - 15
     - .. image:: plots/ShortCalibration/calibration34_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration34_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration34_I_sc_activity_score.png
     	:scale: 20%
   * - 16
     - .. image:: plots/ShortCalibration/calibration36_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration36_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration36_I_sc_activity_score.png
     	:scale: 20%
   * - 17
     - .. image:: plots/ShortCalibration/calibration45_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration45_pep_sc_activity_score.png
     	:scale: 20%
     - .. image:: plots/ShortCalibration/calibration45_I_sc_activity_score.png
     	:scale: 20%

Training set analysis
----------------------

Training set simulations and their performance
...............................................

.. list-table:: Pearson's correlation coefficient for training set simulations (Interface score)
   :widths: 5 20 20
   :header-rows: 1

   * - No.
     - Pearson correlation
     - KS Test
   * - 1
     - * R: -0.22
       * p-value: 0.002
     - * Cutoff: 0.35
       * p-value: 0.008
   * - 2
     - * R: -0.168
       * p-value: 0.020
     - * Cutoff: 0.35
       * p-value: 0.02
   * - 3
     - * R: 0.003
       * p-value: 0.96
     - * Cutoff: 0.35
       * p-value: 0.001
   * - 4
     - * R: -0.21
       * p-value: 0.004
     - * Cutoff: 0.28
       * p-value: 0.0004
   * - 5
     - * R: -0.08
       * p-value: 0.27
     - * Cutoff: 0.22
       * p-value: 0.13
   * - 6
     - * R: -0.22
       * p-value: 0.002
     - * Cutoff: 0.35
       * p-value: 0.0005
   * - 7
     - * R: -0.27
       * p-value: 0.0002
     - * Cutoff: 0.35
       * p-value: 0.007

.. list-table:: Pearson's correlation coefficient for training set simulations (Peptide score)
   :widths: 5 20 20
   :header-rows: 1

   * - No.
     - Pearson correlation
     - KS Test
   * - 1
     - * R: -0.15
       * p-value: 0.04
     - * Cutoff: 0.44
       * p-value: 0.0001
   * - 2
     - * R: -0.13
       * p-value: 0.06
     - * Cutoff: 0.53
       * p-value: 0.0003
   * - 3
     - * R: -0.1
       * p-value: 0.14
     - * Cutoff: 0.03
       * p-value: 0.02
   * - 4
     - * R: -0.14
       * p-value: 0.04
     - * Cutoff: 0.35
       * p-value: :raw-math:`$$ 4.63 \times 10^{-5} $$`
   * - 5
     - * R: -0.21
       * p-value: 0.004
     - * Cutoff: 0.63
       * p-value: 0.002
   * - 6
     - * R: -0.15
       * p-value: 0.03
     - * Cutoff: 0.44
       * p-value: :raw-math:`$$ 1.51 \times 10^{-5} $$`
   * - 7
     - * R: -0.15
       * p-value: 0.03
     - * Cutoff: 0.44
       * p-value: :raw-math:`$$ 2.79 \times 10^{-5} $$`

.. list-table:: Pearson's correlation coefficient for training set simulations (Reweighted score)
   :widths: 5 20 20
   :header-rows: 1

   * - No.
     - Pearson correlation
     - KS Test
   * - 1
     - * R: -0.09
       * p-value: 0.2
     - * Cutoff: 0.31
       * p-value: 0.0005
   * - 2
     - * R: -0.03
       * p-value: 0.68
     - * Cutoff: 0.09
       * p-value: 0.04
   * - 3
     - * R: 0.004
       * p-value: 0.95
     - * Cutoff: 0.52
       * p-value: 0.15
   * - 4
     - * R: -0.08
       * p-value: 0.04
     - * Cutoff: 0.31
       * p-value: 0.003
   * - 5
     - * R: -0.02
       * p-value: 0.7
     - * Cutoff: 0.31
       * p-value: 0.017
   * - 6
     - * R: -0.07
       * p-value: 0.28
     - * Cutoff: 0.31
       * p-value: 0.0015
   * - 7
     - * R: -0.09
       * p-value: 0.19
     - * Cutoff: 0.31
       * p-value: 0.0005
       
--------------------------------------



 .. list-table:: Pearson's correlation coefficient and KS-test values for training set simulations after a clustering step (Interface score)
   :widths: 5 20 20
   :header-rows: 1

   * - No.
     - Pearson correlation
     - KS Test
   * - 1
     - * R: -0.25
       * p-value: 0.002
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 4.4 \times 10^{-6} $$`
   * - 2
     - * R: -0.187
       * p-value: 0.012
     - * Cutoff: 0
       * p-value: 0.005
   * - 3
     - * R: 0.005
       * p-value: 0.84
     - * Cutoff: 0.363
       * p-value: 0.02
   * - 4
     - * R: -0.24
       * p-value: 0.0007
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 4.48 \times 10^{-7} $$`
   * - 5
     - * R: -0.04
       * p-value: 0.55
     - * Cutoff: 0.09
       * p-value: 0.14
   * - 6
     - * R: -0.28
       * p-value: 0.0001
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 2.64 \times 10^{-6} $$`
   * - 7
     - * R: -0.27
       * p-value: 0.00017
     - * Cutoff: 0.31
       * p-value: :raw-math:`$$ 1.53 \times 10^{-6} $$`

.. list-table:: Pearson's correlation coefficient and KS-test values for training set simulations after a clustering step (Peptide score)
   :widths: 5 20 20
   :header-rows: 1

   * - No.
     - Pearson correlation
     - KS Test
   * - 1
     - * R: -0.22
       * p-value: 0.003
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 2.64 \times 10^{-6} $$`
   * - 2
     - * R: -0.17
       * p-value: 0.02
     - * Cutoff: 0
       * p-value: :raw-math:`$$ 4.03 \times 10^{-6} $$`
   * - 3
     - * R: -0.1
       * p-value: 0.167
     - * Cutoff: 0.11
       * p-value: 0.05
   * - 4
     - * R: -0.214
       * p-value: 0.003
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 5.89 \times 10^{-7} $$`
   * - 5
     - * R: -0.126
       * p-value: 0.09
     - * Cutoff: 0.18
       * p-value: :raw-math:`$$ 1.82 \times 10^{-5} $$`
   * - 6
     - * R: -0.24
       * p-value: 0.001
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 2.64 \times 10^{-6} $$`
   * - 7
     - * R: -0.23
       * p-value: 0.001/
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 4.4 \times 10^{-6} $$`

.. list-table:: Pearson's correlation coefficient and KS-test values for training set simulations after a clustering step (Reweighted score)
   :widths: 5 20 20
   :header-rows: 1

   * - No.
     - Pearson correlation
     - KS Test
   * - 1
     - * R: -0.2
       * p-value: 0.007
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 4.4 \times 10^{-6} $$`
   * - 2
     - * R: 0.09
       * p-value: 0.18
     - * Cutoff: 0
       * p-value: 0.01
   * - 3
     - * R: 0.005
       * p-value: 0.938
     - * Cutoff: 0.44
       * p-value: 0.14
   * - 4
     - * R: -0.215
       * p-value: 0.003
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 5.9 \times 10^{-7} $$`
   * - 5
     - * R: -0.08
       * p-value: 0.24
     - * Cutoff: 0.31
       * p-value: 0.006
   * - 6
     - * R: -0.234
       * p-value: 0.001
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 4.81 \times 10^{-6} $$`
   * - 7
     - * R: -0.217
       * p-value: 0.003
     - * Cutoff: 0.34
       * p-value: :raw-math:`$$ 7.27 \times 10^{-6} $$`

Score vs. Activity plots
.........................


.. list-table:: Training set - score vs. activity plots
   :widths: 5 30 30 30
   :header-rows: 1

   * - No.
     - Reweighted Score
     - Peptide Score
     - Interface Score
   * - 1
     - .. image:: plots/TrainingSetAnalysis/calibration16_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration16_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration16_I_sc_activity_score.png
     	:scale: 21%     
   * - 2
     - .. image:: plots/TrainingSetAnalysis/calibration18_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration18_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration18_I_sc_activity_score.png
     	:scale: 21%    
   * - 3
     - .. image:: plots/TrainingSetAnalysis/calibration33_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration33_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration33_I_sc_activity_score.png
     	:scale: 21%     
   * - 4
     - .. image:: plots/TrainingSetAnalysis/calibration38_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration38_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration38_I_sc_activity_score.png
     	:scale: 21%     
   * - 5
     - .. image:: plots/TrainingSetAnalysis/calibration39_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration39_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration39_I_sc_activity_score.png
     	:scale: 21%   
   * - 6
     - .. image:: plots/TrainingSetAnalysis/calibration42_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration42_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration42_I_sc_activity_score.png
     	:scale: 21%     
   * - 7
     - .. image:: plots/TrainingSetAnalysis/calibration43_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration43_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/calibration43_I_sc_activity_score.png
     	:scale: 21%     
     	

.. list-table:: Training set - score vs. activity plots after clustering
   :widths: 5 30 30 30
   :header-rows: 1

   * - No.
     - Reweighted Score
     - Peptide Score
     - Interface Score
   * - 1
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration16_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration16_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration16_I_sc_activity_score.png
     	:scale: 21%     
   * - 2
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration18_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration18_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration18_I_sc_activity_score.png
     	:scale: 21%    
   * - 3
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration33_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration33_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration33_I_sc_activity_score.png
     	:scale: 21%     
   * - 4
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration38_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration38_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration38_I_sc_activity_score.png
     	:scale: 21%     
   * - 5
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration39_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration39_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration39_I_sc_activity_score.png
     	:scale: 21%   
   * - 6
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration42_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration42_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration42_I_sc_activity_score.png
     	:scale: 21%     
   * - 7
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration43_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration43_pep_sc_activity_score.png
     	:scale: 21%
     - .. image:: plots/TrainingSetAnalysis/Clustering/calibration43_I_sc_activity_score.png
     	:scale: 21%

HDAC8 and CdLS syndrome
------------------------

.. table:: Additional putative deacetylation sites for SMC3 suggested by our protocol.

	========================	===========	=============	
	Position of K(ac)		Sequence	FPBind score		
	------------------------	-----------	-------------
		157			RLK(ac)LLR	-1.665
		215			YQK(ac)WDK	-2.082
		304			RTK(ac)LEL	-3.588
		1046			FQK(ac)LVP	-3.957
		105			GAK(ac)KDQ	-4.027
		621			FDK(ac)AFK	-4.050
		400			ELK(ac)SLD	-4.140
		1012			GYK(ac)SIM	-4.619
		388			TSK(ac)EER	-4.747
		493			EKK(ac)QQL	-4.976
		984			VNK(ac)KAL	-5.243
		745			KEK(ac)RQQ	-6.122
		138			IVK(ac)QGK	-6.222
		695			EAK(ac)LNE	-6.646
		1105			TGK(ac)QGE	-6.986
		1052			GGK(ac)ATL	-7.044
	========================	===========	=============

..

References
===========

.. [1] Vannini A, Volpari C, Gallinari P, et al. Substrate binding to histone deacetylases as shown by the crystal structure of the HDAC8-substrate complex. EMBO Rep. 2007;8(9):879-84.
.. [2] Dowling DP, Gantt SL, Gattis SG, Fierke CA, Christianson DW. Structural studies of human histone deacetylase 8 and its site-specific variants complexed with substrate and inhibitors. Biochemistry. 2008;47(51):13554-63.
.. [3] Somoza JR, Skene RJ, Katz BA, et al. Structural snapshots of human HDAC8 provide insights into the class I histone deacetylases. Structure. 2004;12(7):1325-34.
.. [4] English, T. (2004) No More Lunch: Analysis of Sequential Search, Proceedings of the 2004 IEEE Congress on Evolutionary Computation, pp. 227–234.
.. [5] Yanover C, Bradley P. Extensive protein and DNA backbone sampling improves structure-based specificity prediction for C2H2 zinc fingers. Nucleic Acids Res. 2011;39(11):4564-76.
.. [6] Kortemme T, Morozov AV, Baker D. An orientation-dependent hydrogen bonding potential improves prediction of specificity and structure for proteins and protein-protein complexes. J. Mol. Biol. 2003;326:1239-1259.
.. [7] London N, Gullá S, Keating AE, Schueler-furman O. In silico and in vitro elucidation of BH3 binding specificity toward Bcl-2. Biochemistry. 2012;51(29):5841-50.
.. [8] London N, Lamphear CL, Hougland JL, Fierke CA, Schueler-furman O. Identification of a novel class of farnesylation targets by structure-based modeling of binding specificity. PLoS Comput Biol. 2011;7(10):e1002170.
.. [9] Berger B, Leighton T. Protein folding in the hydrophobic-hydrophilic (HP) model is NP-complete. J Comput Biol. 1998;5(1):27-40.
.. [10] Baldwin RL. Energetics of protein folding. J Mol Biol. 2007;371(2):283-301.
.. [11] Kauzmann W. Some factors in the interpretation of protein denaturation. Adv Protein Chem. 1959;14:1-63.
.. [12] Gront D, Kulp DW, Vernon RM, Strauss CE, Baker D. Generalized fragment picking in Rosetta: design, protocols and applications. PLoS ONE. 2011;6(8):e23294.
.. [13] Raveh B, London N, Zimmerman L, Schueler-furman O. Rosetta FlexPepDock ab-initio: simultaneous folding, docking and refinement of peptides onto their receptors. PLoS ONE. 2011;6(4):e18934.
.. [14] Schueler-furman O, Wang C, Bradley P, Misura K, Baker D. Progress in modeling of protein structures and interactions. Science. 2005;310(5748):638-42.
.. [15] Raveh B, London N, Schueler-furman O. Sub-angstrom modeling of complexes between flexible peptides and globular proteins. Proteins. 2010;78(9):2029-40.
.. [16] Cesareni G, Panni S, Nardelli G, Castagnoli L. Can we infer peptide recognition specificity mediated by SH3 domains?. FEBS Lett. 2002;513(1):38-44.
.. [17] Niv MY, Weinstein H. A flexible docking procedure for the exploration of peptide binding selectivity to known structures and homology models of PDZ domains. J Am Chem Soc 2005;127:14072– 14079.
.. [18] Davidon WC. Variable metric method for minimization. SIAM Journal on Optim 1991;1:1–17.
.. [19] Rohl CA, Strauss CE, Misura KM, Baker D. Protein structure pre- diction using Rosetta. Methods Enzymol 2004;383:66–93.
.. [20] Nasmyth K, Haering CH. Cohesin: its roles and mechanisms. Annu Rev Genet. 2009;43:525-58.
.. [21] Dorsett D. Cohesin: genomic insights into controlling gene transcription and development. Curr Opin Genet Dev. 2011;21(2):199-206.
.. [22] Choudhary C, Kumar C, Gnad F, et al. Lysine acetylation targets protein complexes and co-regulates major cellular functions. Science. 2009;325(5942):834-40.
.. [23] Deardorff MA, Bando M, Nakato R, et al. HDAC8 mutations in Cornelia de Lange syndrome affect the cohesin acetylation cycle. Nature. 2012;489(7415):313-7.
.. [24] Deardorff MA, Kaur M, Yaeger D, et al. Mutations in cohesin complex members SMC3 and SMC1A cause a mild variant of cornelia de Lange syndrome with predominant mental retardation. Am J Hum Genet. 2007;80(3):485-94.
.. [25] Lupas A, Van dyke M, Stock J. Predicting coiled coils from protein sequences. Science. 1991;252(5009):1162-4.

.. footer::
	Page ###Page### of ###Total###
