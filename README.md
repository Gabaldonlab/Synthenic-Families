Pipeline to classify lncRNA from different species into families according to synthenic relationships.

This pipeline was developed to classify genes from 4 nematode species: C. elegans, C. briggsae, C. remanei and C. brenneri, but it could be easily modified to be used in any other species.

REQUIREMENTS: ftp://ftp.wormbase.org/pub/wormbase/releases/WS248/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS248.orthologs.txt


###### 1st: for each species, to obtain a file with a sorted list of gene IDs. In this list all chromosomes are concatenated. 
	#ex:
	#WBGene00189949 
	#WBGene00021681 
	#WBGene00004274 
	#XLOC_001892 
	#WBGene00199899 
	#WBGene00199597 
NOTE: this list can contain all genes in each genome, only protein coding genes, etc. 


###### 2nd: to create a file including all orthology relationships for the four nematode species
#usage: python wormbase_orthoParalogsGH.py elegans_genes.list c_elegans.PRJNA13758.WS248.orthologs.txt celegans_ortho.out
NOTE: protein coding genes are more prone to have orthology relationships annotated


###### 3rd: to compare gene order between species.The script script produces a file including all lncRNA having conserved syntheny with lncRNA from any of the three other species 
#usage: synteny_nematodesv4GH.py sp1_geneOrder.list sp2_geneOrder.list sp3_geneOrder.list sp4_geneOrder.list celegans_ortho.out 4spv4.out cluster overlap minSideGenes noHomology

	- cluster (integer): number of genes considered at each side of a given lncRNA; ex: 3. NOTE: if stated 3, the considered size of the cluster is 3+3=6
	- overlap (integer): minimum number of shared genes for each pairwise comparisson betwen species; ex: 3
	- minSideGenes (integer): minimum number of shared genes at each side of a given lncRNA to be considered as members of the same family; ex: 1
	- noHomology (yes/no): yes: to consider genes lacking orthology relationships; no: do not consider genes lacking homology	
	ex: python synteny_nematodesv4GH.py celegans_geneOrder.list cbriggsae_geneOrder.list cbrenneri_geneOrder.list cremanei_geneOrder.list celegans_ortho.out out 3 3 1 no


###### 4th: to classify lncRNA from 4 different species into families
#usage: python classifyFamiliesv5_VennGH.py cele cbrig cbren crem 4spv4.out 4spv4_Families.fam 4spv4_Families.txt 4spv4_FAMvenn.R 4spv4_GENESvenn.R >4spv4_Families.counts

	- ex INPUT file: 4spv4.out
	cele	cbren	XLOC_014828	XLOC_024666
	cele	cbren	XLOC_014828	XLOC_024657
	- ex INPUT file: cele -> #number of C. elegans genes to classify
	- ex INPUT file: cbrig -> #number of C. briggsae genes to classify
	- ex INPUT file: cbren -> #number of C. brenneri genes to classify
	- ex INPUT file: crem -> #number of C. remanei genes to classify

	- ex OUTPUT file: 4spv4_Families.fam
	fam1    XLOC_010119cbrig
	fam1	XLOC_019544crem
	fam1	XLOC_024997cbren
	- ex OUTPUT file: 4spv4_Families.txt
    	fam1	  bren	brig	 rem
    	fam2   	bren   rem
	- ex OUTPUT file: 4spv4_Families.counts
    	EleBrenBrig 56
    	EleBrenRem 51
	- ex OUTPUT file: 4spv4_GENESvenn.R -> R script to draw a venn diagram with the number of overlapped genes using R
	- ex OUTPUT file: 4spv4_FAMvenn.R -> R script to draw a venn diagram with the number of overlapped families using R

###### 5th: to refine families classification using blast matchs
the script calculates a score indicating the blast suport for a given family baseed on blast hits
#usage: python refineFamiliesWithBlastGH.py 4spv4_Families.fam blastOut 4spGeneTransID.txt out1 out2
	- ex INPUT file: 4spv4_Families.fam
	fam1    XLOC_010119cbrig
	fam1	XLOC_019544crem
	fam1	XLOC_024997cbren
	- ex INPUT file: blastOut.txt
   	TCONS_00000606cbren	TCONS_00000606cbren	100.00	486	0	0	1486	1	486	0.0	 898
	- ex INPUT file: 4spGeneTransID.txt
   	XLOC_000040cele	TCONS_00000067cele
	- ex OUTPUT file: OUT1
    	fam408	XLOC_009936crem	XLOC_010893cbrig	0
    	fam408	XLOC_009936crem	XLOC_012709cbrig	0
	- ex OUTPUT file: OUT2
	famID	  blastHits	possiblePairwiseComparissons	score	#sp
    	fam408	0.0	10.0	0.0	3	rem	brig	ele



