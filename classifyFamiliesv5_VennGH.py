# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:55:13 2016

@author: cpegueroles
"""
#usage: python classifyFamiliesv5_VennGH.py cele cbrig cbren crem 4spv4.out 4spv4_Families.fam 4spv4_Families.txt 4spv4_FAMvenn.R 4spv4_GENESvenn.R >4spv4_Families.counts
#ex: 4spv4.out
#cele    cbren	XLOC_014828	XLOC_024666
#cbren   cele	XLOC_024666	XLOC_014828


import sys,os,re

def renameGenes(in1):#ex: 4spv4_6cluster3min.out
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
cele= sys.argv[1] #number of expressed (protein coding/lncRNA) genes
cbrig= sys.argv[2] #number of expressed (protein coding/lncRNA) genes
cbren= sys.argv[3] #number of expressed (protein coding/lncRNA) genes
crem= sys.argv[4] #number of expressed (protein coding/lncRNA) genes
in1= sys.argv[5] #out file from synteny_nematodesv4GH.py; ex: 4spv4_6cluster3min.out
out= open(sys.argv[6], 'w') #.fam
out2= open(sys.argv[7], 'w') #.txt
out3= open(sys.argv[8], 'w') #number of families
out4= open(sys.argv[9], 'w') #number of genes

#to group into families
dictTemp1= renameGenes(in1)

dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1

dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1
           
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

celem=0; cbrigm=0; cbrenm=0; cremm=0; #number of classified genes
for key, val in dictFam.iteritems():
    for x in val:
        out.write("fam%s\t%s\n" % (key, x))
        if re.search("cele", x):
            celem=celem+1
        if re.search("cbrig", x):
            cbrigm=cbrigm+1
        if re.search("cbren", x):
            cbrenm=cbrenm+1
        if re.search("crem", x):
            cremm=cremm+1

#to count the number of species for gene family        
EleBren=0; EleBrig=0; EleRem=0; BrenBrig=0; BrenRem=0; BrigRem=0
EleBrenBrig=0; EleBrenRem=0; EleBrigRem=0; BrenBrigRem=0
EleBrenBrigRem=0
genesEleBren=0; genesEleBrig=0; genesEleRem=0; genesBrenBrig=0; genesBrenRem=0; genesBrigRem=0
genesEleBrenBrig=0; genesEleBrenRem=0; genesEleBrigRem=0; genesBrenBrigRem=0
genesEleBrenBrigRem=0
EleF=0; BrigF=0; BrenF=0; RemF=0
EleG=0; BrigG=0; BrenG=0; RemG=0
for key, val in dictFam.iteritems():
    mySp=[]
    for x in val:
        x=x.split('c')
        sp =x[1]
        mySp.append(sp)
    out2.write("fam%s\t%s\n" % (key,'\t'.join(list(set(mySp)))))
    if 'ele' in mySp and 'brig' in mySp and 'bren' in mySp and 'rem' in mySp:
        EleBrenBrigRem +=1
        genesEleBrenBrigRem = genesEleBrenBrigRem+len(mySp)
    if 'ele' in mySp and 'brig' in mySp and 'bren' in mySp and not 'rem' in mySp:
        EleBrenBrig +=1
        genesEleBrenBrig =genesEleBrenBrig+len(mySp)
    if 'ele' in mySp and 'rem' in mySp and 'bren' in mySp and not 'brig' in mySp:
        EleBrenRem +=1
        genesEleBrenRem = genesEleBrenRem+len(mySp)
    if 'ele' in mySp and 'rem' in mySp and 'brig' in mySp and not 'bren' in mySp:
        EleBrigRem +=1
        genesEleBrigRem = genesEleBrigRem+len(mySp)
    if 'bren' in mySp and 'rem' in mySp and 'brig' in mySp and not 'ele' in mySp:
        BrenBrigRem +=1
        genesBrenBrigRem = genesBrenBrigRem+len(mySp)
    if 'ele' in mySp and 'bren' in mySp and not 'rem' in mySp and not 'brig' in mySp :
        EleBren +=1
        genesEleBren = genesEleBren+len(mySp)
    if 'ele' in mySp and 'brig' in mySp and not 'bren' in mySp and not 'rem' in mySp:
        EleBrig +=1
        genesEleBrig = genesEleBrig+len(mySp)
    if 'ele' in mySp and 'rem' in mySp and not 'brig' in mySp and not 'bren' in mySp:
        EleRem +=1
        genesEleRem = genesEleRem+len(mySp)
    if 'bren' in mySp and 'brig' in mySp and not 'ele' in mySp and not 'rem' in mySp:
        BrenBrig +=1
        genesBrenBrig= genesBrenBrig+len(mySp)
    if 'bren' in mySp and 'rem' in mySp and not 'brig' in mySp and not 'ele' in mySp:
        BrenRem +=1
        genesBrenRem = genesBrenRem+len(mySp)
    if 'rem' in mySp and 'brig' in mySp and not 'ele' in mySp and not 'bren' in mySp:
        BrigRem +=1
        genesBrigRem = genesBrigRem+len(mySp)
    if 'ele' in mySp:
        EleF +=1; EleG += len(mySp)
    if 'brig' in mySp:
        BrigF +=1; BrigG += len(mySp)
    if 'bren' in mySp:
        BrenF +=1; BrenG += len(mySp)
    if 'rem' in mySp:
        RemF +=1; RemG += len(mySp)
        
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
out3.write('ele <-%s\n' % (int(cele)-(celem -EleF)))
out3.write('bren <-%s\n' % (int(cbren)-(cbrenm-BrenF)))
out3.write('brig <-%s\n' % (int(cbrig)-(cbrigm-BrigF)))
out3.write('rem <-%s\n' % (int(crem)-(cremm-RemF)))
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


out4.write('rm(list=ls())\n')
out4.write('library(VennDiagram)\n')
out4.write('setwd("%s")\n' % os.getcwd())
out4.write('ele <-%s\n' % ((int(cele)-celem)+EleG))
out4.write('bren <-%s\n' % ((int(cbren)-cbrenm)+BrenG))
out4.write('brig <-%s\n' % ((int(cbrig)-cbrigm)+BrigG))
out4.write('rem <-%s\n' % ((int(crem)-cremm)+RemG))
out4.write('overlapEleBren <- %s\n' % (genesEleBren + genesEleBrenBrig + genesEleBrenRem + genesEleBrenBrigRem))
out4.write('overlapEleBrig <- %s\n' % (genesEleBrig + genesEleBrigRem + genesEleBrenBrig + genesEleBrenBrigRem))
out4.write('overlapEleRem <- %s\n' % (genesEleRem + genesEleBrenRem + genesEleBrigRem + genesEleBrenBrigRem))
out4.write('overlapBrenBrig <- %s\n' % (genesBrenBrig + genesEleBrenBrig + genesBrenBrigRem + genesEleBrenBrigRem))
out4.write('overlapBrenRem <- %s\n' % (genesBrenRem + genesEleBrenRem + genesBrenBrigRem + genesEleBrenBrigRem))
out4.write('overlapBrigRem <- %s\n' % (genesBrigRem + genesEleBrigRem + genesBrenBrigRem + genesEleBrenBrigRem))
out4.write('overlapEleBrenBrig <- %s\n' % (genesEleBrenBrig + genesEleBrenBrigRem))
out4.write('overlapEleBrenRem <- %s\n' % (genesEleBrenRem + genesEleBrenBrigRem))
out4.write('overlapEleBrigRem <- %s\n' % (genesEleBrigRem + genesEleBrenBrigRem))
out4.write('overlapBrenBrigRem <- %s\n' % (genesBrenBrigRem + genesEleBrenBrigRem))
out4.write('overlapEleBrenBrigRem <- %s\n' % (genesEleBrenBrigRem))
out4.write('draw.quad.venn(area1=ele, area2=bren, area3=brig, area4=rem,\n n12=overlapEleBren, n13=overlapEleBrig, n14=overlapEleRem,\n n23=overlapBrenBrig, n24=overlapBrenRem, n34=overlapBrigRem,\n n123=overlapEleBrenBrig, n124=overlapEleBrenRem, n134=overlapEleBrigRem,\n n234=overlapBrenBrigRem, n1234=overlapEleBrenBrigRem,\n fill = c("skyblue", "pink1", "mediumorchid", "orange"))\n')