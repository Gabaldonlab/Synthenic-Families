# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 13:51:53 2016

@author: cpegueroles
"""

#script to calculate a score indicating the blast suport for a given family
#ex 4spv4_6cluster3minLeftRight_Families.fam:
	#fam1	XLOC_009958crem
#ex 4sp_blastOut.txt:
   #TCONS_00000606cbren	TCONS_00000606cbren	100.00	486	0	0	1486	1	486	0.0	 898
#ex 4spGeneTransID.txt
   #XLOC_000040cele	TCONS_00000067cele
#ex OUT1:
    #fam408	XLOC_009936crem	XLOC_010893cbrig	0
    #fam408	XLOC_009936crem	XLOC_012709cbrig	0
#ex OUT2:
    #famID	  blastHits	possiblePairwiseComparissons	score	#sp
    #fam408	0.0	10.0	0.0	3	rem	brig	ele

import sys
import itertools
#a = [1,2,3,4]
#print list(itertools.combinations(a,2))

in1= sys.argv[1] #ex: 4spv4_6cluster3minLeftRight_Families.fam
in2= sys.argv[2] #ex: 4sp_blastOut.txt
in3= sys.argv[3] #ex: 4spGeneTransID.txt 
temp= open(sys.argv[4], 'w') #
temp.write('famID\tsp1\tsp2\tblastHit\n')
out= open(sys.argv[5], 'w') #
out.write('famID\tblastHits\tpossiblePairwiseComparissons\tscore\t#sp\n')

#create a dictionary of all families
dictFam={}
for row in open(in1, 'r').readlines():
    row = row.rstrip().split('\t')
    key= row[0]
    val=[]
    if not key in dictFam:
        dictFam[key]=val
    dictFam[key].append(row[1])   
#to retrieve all possible pairs within a given family. 
#Direction is not important, pairs are not repeated; ex: [A B] and [B A] are consirered once as [A B]
dictFP={}
for key, val in dictFam.iteritems():
    keyFP=key
    valFP=list(itertools.combinations(val,2))
    dictFP[keyFP]=valFP
    
#create a dictionary with blast matchs; dictionary of geneIds instead of transcriptID
dictBlastMatch={}
for row1 in open(in2, 'r').readlines():#4sp_blastOut.txt
    row1=row1.rstrip().split('\t')
    first=""
    second=""
    for row2 in open(in3, 'r').readlines():#4spGeneTransID.txt
        i=0
        row2=row2.rstrip().split('\t')
        if row1[0] ==row2[1]:
            first=row2[0]
            i=i+1
        if row1[1] ==row2[1]:
            second=row2[0]  
            i=i+1
        if i==2:
            break
    if first != second:
        key=first
        val=second
        if not first in dictBlastMatch.keys() and not first in dictBlastMatch.values() : 
            dictBlastMatch[key]=val

finalDict={}
for keyFP, valFP in dictFP.iteritems():
    finalKey=keyFP
    finalval=[] 
    finalDict[finalKey]=finalval
    for x in valFP:
        found=0
        for key, val in dictBlastMatch.iteritems(): 
            if key in x and val in x:
                #print keyFP,x,1
                finalval.append([x[0], x[1], 1])
                found=1
        if found==0:
            #print keyFP,x,0
            finalval.append([x[0], x[1], 0])

for key, val in finalDict.iteritems(): 
    mySp=[]
    for x in val:
        #print key, x[0], x[1], x[2]
        temp.write('%s\t%s\t%s\t%s\n' % (key, x[0], x[1], x[2]))
        sp1=x[0].split('c')
        sp1=sp1[1]
        sp2=x[1].split('c')
        sp2=sp2[1]
        mySp.append(sp1)
        mySp.append(sp2)
    score=float(sum(x[2] for x in val))/float(len(val))
    mySp= list(set(mySp))
    #print key, score
    out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (key,sum(x[2] for x in val),len(val), score,len(mySp),'\t'.join(list(set(mySp)))))


    