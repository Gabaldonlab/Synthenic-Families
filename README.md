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


###### 2nd: to create a file including all orthology relationships for the 4 nemateode species
	ex: python wormbase_orthoParalogsGH.py elegans_genes.list c_elegans.PRJNA13758.WS248.orthologs.txt celegans_ortho.out
NOTE: protein coding genes are more prone to have orthology relationships annotated


###### 3rd: to compare gene order between species.The script script produces a file including all lncRNA having conserved syntheny with lncRNA from any of the three other species 
usage: synteny_nematodesv4_allGenes_LeftRigth.py sp1_geneOrder.list sp2_geneOrder.list sp3_geneOrder.list sp4_geneOrder.list celegans_ortho.out 4spv4.out cluster overlap minSideGenes noHomology
	- cluster (integer): number of genes considered at each side of a given lncRNA; ex: 3. NOTE: if stated 3, the considered size of the cluster is 3+3=6
	- overlap (integer): minimum number of shared genes for each pairwise comparisson betwen species; ex: 3
	- minSideGenes (integer): minimum number of shared genes at each side of a given lncRNA to be considered as members of the same family; ex: 1
	- noHomology (yes/no): yes: to consider genes lacking orthology relationships; no: do not consider genes lacking homology	
	ex: python synteny_nematodesv4GH.py celegans_geneOrder.list cbriggsae_geneOrder.list cbrenneri_geneOrder.list cremanei_geneOrder.list celegans_ortho.out out 3 3 1 no


###### 4th: to classify lncRNA from 4 different species into families
#usage: python classifyFamiliesv3_VennGH.py 4spv4.out 4spv4_Families.fam 4spv4_Families.txt 4spv4_venn.R >4spv4_Families.counts
	- ex INPUT file: 4spv4.out
	cele	cbren	XLOC_014828	XLOC_024666
	cele	cbren	XLOC_014828	XLOC_024657
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
	- ex OUTPUT file: 4spv4_venn.R -> R script to draw a venn diagram using R



