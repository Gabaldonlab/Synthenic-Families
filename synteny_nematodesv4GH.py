# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:52:45 2016
@author: cinta
"""
#script to obtain a file including all lncRNA having conserved syntheny with lncRNA from any of the three other species 
#to run: celegans_geneOrder.list celegans_ortho.out out
#celegans_geneOrder.list: contains the gene order of the whole genome concatenating chr 
#head celegans_ortho.out: contains othology relationships (orthologs and paralogs)
    #elegans	remanei	briggsae	brenneri  myID
    #WBGene00022277	WBGene00057179	WBGene00026679	WBGene00151451	myID1
#ex out:
    #cele   	cbren	XLOC_014828	XLOC_024666
    #cele	  cbren	XLOC_014828	XLOC_024657

#time python synteny_nematodesv4_allGenes_LeftRigth.py 
#celegans_geneOrder.list cbriggsae_geneOrder.list cbrenneri_geneOrder.list cremanei_geneOrder.list 
#celegans_ortho.out 4spv4_6cluster2minFull.out  3 3 1 no
#	3: number of genes considered at each side of a given lncRNA
#	3: minimum number of shared genes 
#	1: minimum number of shared genes at each side
#	yes -> to consider those genes without homology; no -> only genes with homology

import sys,re

def recodeList(file1, file2, string):#gene order list for spX, orthology file, ex: 'cele'
    #returns a list with the gene order of a given sp were geneIDs were renamed to myID 
    #to allow direct comparisson between species;
    #genes with no homology are also included!
    #ex: ['myID1', 'WBGene000xxx' 'myID2', 'myID3', 'myID4', 'myID5','XLOC_001892', 'myID6', 'myID7', ...]  
    mylist = []
    matched=0; nonmatched=0; found=0; lncRNA=0
    for row1 in open(file1, 'r'):  
        row1 = row1.rstrip()
        if re.search(r'XLOC_', row1):
            mylist.append(row1)
            found=1
            lncRNA=lncRNA+1
        for row2 in open(file2 , 'r').readlines():  
            row2= row2.rstrip().split("\t")
            if string=='cele':
                if row1== row2[0]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='cbrig':
                if row1== row2[2]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='cbren':
                if row1== row2[3]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='crem':
                if row1== row2[1]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
                    
        if found==0:
            nonmatched=nonmatched+1
            #TO IMPROVE: If found==0 I could append row1 (a gene with no orthology relationship) to mylist
            if nonMatch=='yes':
                mylist.append(row1)
        found=0
        #print mylist
    print "number of matched IDs for",string,": " ,matched
    print "number of NON matched IDs for",string,": " ,nonmatched
    print "number of candidate lncRNA",string,": " ,lncRNA,"\n" 
    return mylist

def dictionaryOfClusters(myidx,mylist):#positions in myList for each lncRNA; myList=renamed gene order
    #to create a dict; for each lncRNA (key) stores the number of nearby renamed geneID indicated by the user
    #TO IMPROVE: I could not count lncRNAID...    
    mydict={}
    for idx in myidx:
        key=mylist[idx]
        val={'left':[], 'right':[], 'all':[]}
        if not key in mydict:
            mydict[key]=val
        i=1
        while i <=genesNearby:
            try:
                mydict[key]['left'].append(mylist[idx-i])
                mydict[key]['right'].append(mylist[idx+i])
                mydict[key]['all'].append(mylist[idx-i])
                mydict[key]['all'].append(mylist[idx+i])
            except IndexError:
                pass
            i=i+1
        #print key, mydict[key]
            
    return mydict
    
def comparingDict(sp1, sp2):
    #compares dictionaries from dictionaryOfClusters() for two species; 
    #if the number of shared genes is >= minOverlap, lncRNA are stored in myHomologs list
    if sp1=='cele':
        dict1= dictionaryOfClusters(cele_idx,celeList)
    if sp1=='cbrig':
        dict1= dictionaryOfClusters(cbrig_idx,cbrigList)
    if sp1=='cbren':
        dict1= dictionaryOfClusters(cbren_idx,cbrenList)
    if sp1=='crem':
        dict1= dictionaryOfClusters(crem_idx,cremList)
    if sp2=='cele':
        dict2= dictionaryOfClusters(cele_idx,celeList)
    if sp2=='cbrig':
        dict2= dictionaryOfClusters(cbrig_idx,cbrigList)
    if sp2=='cbren':
        dict2= dictionaryOfClusters(cbren_idx,cbrenList)
    if sp2=='crem':
        dict2= dictionaryOfClusters(crem_idx,cremList)
    myHomologs=[]
    homologFound='false'
    for key1, val1 in dict1.iteritems():
        print key1, val1
        for key2, val2 in dict2.iteritems():
            #if len(set(val1).intersection(val2)) >=minOverlap:
                #mytup=(key1,key2)
                #myHomologs.append(mytup)
            if len(set(dict1[key1]['all']).intersection(dict2[key2]['all'])) >=minOverlap:
            #at least one of the two lncRNA share genes in the left and the right side    
#                if len(set(dict1[key1]['right']).intersection(dict2[key2]['right'])) >=minSideOverlap or len(set(dict1[key1]['right']).intersection(dict2[key2]['left'])) >=minSideOverlap:
#                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['right'])) >=minSideOverlap or len(set(dict1[key1]['left']).intersection(dict2[key2]['left'])) >=minSideOverlap:
            #forces that both lncRNA share genes in the left and the right side    
                if len(set(dict1[key1]['right']).intersection(dict2[key2]['right'])) >=minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['left'])) >=minSideOverlap: 
                        homologFound='true'
                if len(set(dict1[key1]['right']).intersection(dict2[key2]['left'])) >=minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['right'])) >=minSideOverlap:
                        homologFound='true'
                if homologFound=='true':
                        #print dict1[key1], dict2[key2]               
                        mytup=(key1,key2)
                        myHomologs.append(mytup)
                        homologFound='false'
    return myHomologs

#####################################

in1= sys.argv[1] #gene order list for sp1
in2= sys.argv[2] #gene order list for sp2
in3= sys.argv[3] #gene order list for sp3
in4= sys.argv[4] #gene order list for sp4
in5= sys.argv[5] #orthology file
out= open(sys.argv[6] , 'w') #out file
genesNearby= int(sys.argv[7])
minOverlap= int(sys.argv[8])
minSideOverlap= int(sys.argv[9])
nonMatch = sys.argv[10]
temp= open('temp', 'w')

celeList=recodeList(in1, in5, 'cele')
cbrigList=recodeList(in2, in5, 'cbrig')
cbrenList=recodeList(in3, in5, 'cbren')
cremList=recodeList(in4, in5, 'crem')

for x in celeList:
    temp.write("%s\telegans\n" %(x))       
for x in cbrigList :
    temp.write("%s\tbriggsae\n" %(x))    
for x in cbrenList:
    temp.write("%s\tbrenneri\n" %(x))
for x in cremList:
    temp.write("%s\tremanei\n" %(x))

cele_idx = [i for i, item in enumerate(celeList) if item.startswith('XLOC_')]
#list containing the positions in celeList for each lncRNA
cbrig_idx = [i for i, item in enumerate(cbrigList) if item.startswith('XLOC_')]
cbren_idx = [i for i, item in enumerate(cbrenList) if item.startswith('XLOC_')]
crem_idx = [i for i, item in enumerate(cremList) if item.startswith('XLOC_')]

for x in comparingDict('cele','cbren'):
    #print 'cele','cbren','\t'.join(x)
    out.write('cele\tcbren\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cele','cbrig'):
    #print 'cele','cbrig','\t'.join(x)
    out.write('cele\tcbrig\t%s\n'% ('\t'.join(x)))     
for x in comparingDict('cele','crem'):
    #print 'cele','crem','\t'.join(x)
    out.write('cele\tcrem\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cbren','cbrig'):
    #print 'cbren','cbrig','\t'.join(x)
    out.write('cbren\tcbrig\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cbren','crem'):     
    #print 'cbren','crem','\t'.join(x)
    out.write('cbren\tcrem\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cbrig','crem'):
    #print 'cbrig','crem','\t'.join(x)
    out.write('cbrig\tcrem\t%s\n'% ('\t'.join(x)))

for x in comparingDict('cbren','cele'):
    out.write('cbren\tcele\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cbrig','cele'):
    out.write('cbrig\tcele\t%s\n'% ('\t'.join(x)))     
for x in comparingDict('crem','cele'):
    out.write('crem\tcele\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cbrig','cbren'):
    out.write('cbrig\tcbren\t%s\n'% ('\t'.join(x)))
for x in comparingDict('crem','cbren'):     
    out.write('crem\tcbren\t%s\n'% ('\t'.join(x)))
for x in comparingDict('crem','cbrig'):
    out.write('crem\tcbrig\t%s\n'% ('\t'.join(x)))




