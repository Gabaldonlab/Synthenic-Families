# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:55:13 2016
@author: cpegueroles
"""
#to classify lncRNA from 4 different species into families
#usage: classifyFamiliesv3_Venn.py 4spv4.out 4spv4_Families.fam 4spv4_Families.txt 4spv4_venn.R >4spv4_Families.counts
#ex INPUT file: 4spv4.out
    #cele	cbren	XLOC_014828	XLOC_024666
    #cele	cbren	XLOC_014828	XLOC_024657
#ex OUTPUT file: 4spv4_Families.fam
    #fam1    XLOC_010119cbrig
    #fam1   XLOC_019544crem
    #fam1	  XLOC_024997cbren
#ex OUTPUT file: 4spv4_Families.txt
    #fam1	  bren	brig	 rem
    #fam2   	bren   rem
#ex OUTPUT file: 4spv4_Families.counts
    #EleBrenBrig 56
    #EleBrenRem 51
#ex OUTPUT file: 4spv4_venn.R -> R script to draw a venn diagram using R


import sys,os

def renameGenes(in1):#ex: 4spv4.out
    #renames genes adding the species; ex: XLOC_010119 -> XLOC_010119cbrig
    #creates a dictionary with sp1 (key) pointing to sp2 (val)
    dictTemp={}
    for row in open(in1, 'r').readlines():
        row = row.rstrip().split('\t')
        key= row[2]+row[0]
        val=[]
        if not key in dictTemp:
            dictTemp[key]=val
        dictTemp[key].append(row[2]+row[0])
        dictTemp[key].append(row[3]+row[1])
    return dictTemp
    
##################################################
in1= sys.argv[1] #out file from synteny_nematodesv4.py; ex: 4spv4.out
out= open(sys.argv[2], 'w')
out2= open(sys.argv[3], 'w')
out3= open(sys.argv[4], 'w')

#to group into families
dictTemp1= renameGenes(in1)
#for key, val in dictTemp1.iteritems():
    #print key, val 
dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1
#####################################             
#to required to include all relationships
dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1
#to ensure all relationships are included   
dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1
#####################################         
#for key, val in dictTemp1.iteritems():
#    print key, val    
           
mylist=[val for val in dictTemp1.values()]
mylistSorted=[]
for x in mylist:
    mylistSorted.append(sorted(x)) 
uniq_mylist = [list(t) for t in set(sorted(map(tuple, mylistSorted)))]
     
     
dictFam={}
#dictFam structure: fam -> {sp1:[,], sp2:[,], sp3:[,], sp4:[,]}
i=1
for x in uniq_mylist:
    key=i
    val=x
    dictFam[key]=val
    i=i+1
    
for key, val in dictFam.iteritems():
    for x in val:
        out.write("fam%s\t%s\n" % (key, x))

#to count the number of species for gene family        
EleBren=0; EleBrig=0; EleRem=0; BrenBrig=0; BrenRem=0; BrigRem=0
EleBrenBrig=0; EleBrenRem=0; EleBrigRem=0; BrenBrigRem=0
EleBrenBrigRem=0
for key, val in dictFam.iteritems():
    mySp=[]
    for x in val:
        x=x.split('c')
        sp =x[1]
        mySp.append(sp)
    #out2.write("fam%s\t%s\n" % (key,len(list(set(mySp)))))
    mySp= list(set(mySp))
    out2.write("fam%s\t%s\n" % (key,'\t'.join(list(set(mySp)))))
    if 'ele' in mySp and 'brig' in mySp and 'bren' in mySp and 'rem' in mySp:
        EleBrenBrigRem +=1
        #print key, val
    if 'ele' in mySp and 'brig' in mySp and 'bren' in mySp and not 'rem' in mySp:
        EleBrenBrig +=1
        #print key, val
    if 'ele' in mySp and 'rem' in mySp and 'bren' in mySp and not 'brig' in mySp:
        EleBrenRem +=1
        #print key, val
    if 'ele' in mySp and 'rem' in mySp and 'brig' in mySp and not 'bren' in mySp:
        EleBrigRem +=1
    if 'bren' in mySp and 'rem' in mySp and 'brig' in mySp and not 'ele' in mySp:
        BrenBrigRem +=1
        #print key, val
    if 'ele' in mySp and 'bren' in mySp and not 'rem' in mySp and not 'brig' in mySp :
        EleBren +=1
    if 'ele' in mySp and 'brig' in mySp and not 'bren' in mySp and not 'rem' in mySp:
        EleBrig +=1
    if 'ele' in mySp and 'rem' in mySp and not 'brig' in mySp and not 'bren' in mySp:
        EleRem +=1
    if 'bren' in mySp and 'brig' in mySp and not 'ele' in mySp and not 'rem' in mySp:
        BrenBrig +=1
    if 'bren' in mySp and 'rem' in mySp and not 'brig' in mySp and not 'ele' in mySp:
        BrenRem +=1
    if 'rem' in mySp and 'brig' in mySp and not 'ele' in mySp and not 'bren' in mySp:
        BrigRem +=1
        

print 'EleBrenBrig',EleBrenBrig
print 'EleBrenRem',EleBrenRem
print 'EleBrigRem', EleBrigRem
print 'BrenBrigRem',BrenBrigRem
print 'EleBren',EleBren
print 'EleBrig',EleBrig
print 'EleRem',EleRem
print 'BrenBrig',BrenBrig 
print 'BrenRem',BrenRem
print 'BrigRem',BrigRem
print 'EleBrenBrigRem',EleBrenBrigRem

out3.write('rm(list=ls())\n')
out3.write('library(VennDiagram)\n')
out3.write('setwd("%s")\n' % os.getcwd())
out3.write('ele <-1298\n')
out3.write('bren <-3175\n')
out3.write('brig <-1186\n')
out3.write('rem <-1458\n')
out3.write('overlapEleBren <- %s\n' % (EleBren + EleBrenBrig + EleBrenRem + EleBrenBrigRem))
out3.write('overlapEleBrig <- %s\n' % (EleBrig + EleBrigRem + EleBrenBrig + EleBrenBrigRem))
out3.write('overlapEleRem <- %s\n' % (EleRem + EleBrenRem + EleBrigRem + EleBrenBrigRem))
out3.write('overlapBrenBrig <- %s\n' % (BrenBrig + EleBrenBrig + BrenBrigRem + EleBrenBrigRem))
out3.write('overlapBrenRem <- %s\n' % (BrenRem + EleBrenRem + BrenBrigRem + EleBrenBrigRem))
out3.write('overlapBrigRem <- %s\n' % (BrigRem + EleBrigRem + BrenBrigRem + EleBrenBrigRem))
out3.write('overlapEleBrenBrig <- %s\n' % (EleBrenBrig + EleBrenBrigRem))
out3.write('overlapEleBrenRem <- %s\n' % (EleBrenRem + EleBrenBrigRem))
out3.write('overlapEleBrigRem <- %s\n' % (EleBrigRem + EleBrenBrigRem))
out3.write('overlapBrenBrigRem <- %s\n' % (BrenBrigRem + EleBrenBrigRem))
out3.write('overlapEleBrenBrigRem <- %s\n' % (EleBrenBrigRem))
out3.write('draw.quad.venn(area1=ele, area2=bren, area3=brig, area4=rem,\n n12=overlapEleBren, n13=overlapEleBrig, n14=overlapEleRem,\n n23=overlapBrenBrig, n24=overlapBrenRem, n34=overlapBrigRem,\n n123=overlapEleBrenBrig, n124=overlapEleBrenRem, n134=overlapEleBrigRem,\n n234=overlapBrenBrigRem, n1234=overlapEleBrenBrigRem,\n fill = c("skyblue", "pink1", "mediumorchid", "orange"))\n')
