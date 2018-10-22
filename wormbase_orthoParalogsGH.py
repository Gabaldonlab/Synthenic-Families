# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:34:27 2016
@author: cinta
"""
#script to create a file including all orthology relationships for the 4 nemateode species
#usage: python wormbase_orthoParalogs.py genes.list wormbase_orthologs_file out
#WARNING: prints all orthologs, including paralogs; paralogs are assigned with the same myID
#example genes.list:
    #WBGene00021681 
    #WBGene00004274 
    #XLOC_001892 
#example wormbase_orthologs_file:
    #ftp://ftp.wormbase.org/pub/wormbase/releases/WS248/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS248.orthologs.txt
#example out: 
    #elegans	remanei	briggsae	brenneri myID
    #WBGene00022277	WBGene00057179	WBGene00026679	WBGene00151451	myID1
    #WBGene00022277	WBGene00057179	WBGene00026679	WBGene00155359	myID1
    #WBGene00022276	WBGene00057180	WBGene00026678	WBGene00143538	myID2

import sys,re

in1= open(sys.argv[1] , 'r') #genes.list
in2= sys.argv[2] #wormbase_orthologs_file
out=open(sys.argv[3],'w') #out
out.write("elegans\tremanei\tbriggsae\tbrenneri\tmyID\n")

mylist=[]
for row in in1.readlines():  
    mylist.append(row.rstrip())

dictOrt={}
new=0;rem=0; brig=0; bren=0
ort= {'rem':[], 'brig':[], 'bren':[]}
for row in open(in2, 'r').readlines():
    row = row.rstrip().split("\t")
    if re.search(r'WBGene0', row[0]):
        key=row[0]
        if not key in dictOrt:
           dictOrt[key]=ort
           new=1
    if re.search(r'Caenorhabditis remanei', row[0]):
        dictOrt[key]['rem'].append(row[1])
        rem=1
    if re.search(r'Caenorhabditis briggsae', row[0]):
        dictOrt[key]['brig'].append(row[1])
        brig=1
    if re.search(r'Caenorhabditis brenneri', row[0]):
        dictOrt[key]['bren'].append(row[1])
        bren=1
    if re.search(r'=', row[0]) and new==1: 
        if rem==0:
            dictOrt[key]['rem'].append('na')
        if brig==0:
            dictOrt[key]['brig'].append('na')    
        if bren==0:
            dictOrt[key]['bren'].append('na')    
        ort= {'rem':[], 'brig':[], 'bren':[]}
        new=0;rem=0; brig=0; bren=0
        
if rem==0:
    dictOrt[key]['rem'].append('na')
if brig==0:
    dictOrt[key]['brig'].append('na')    
if bren==0:
   dictOrt[key]['bren'].append('na')

myID=1
for x in mylist:
   for key, val in dictOrt.iteritems():
        if x == key:
            length=[len(val['rem']),len(val['brig']),len(val['bren'])]
            mymax=max(length)
            i=0
            while i <mymax:
                try:
                    rem= val['rem'][i]
                except IndexError:
                    pass
                try:
                    brig= val['brig'][i]
                except IndexError:
                    pass
                try:
                    bren= val['bren'][i]
                except IndexError:
                    pass
                out.write("%s\t%s\t%s\t%s\tmyID%s\n" % (key, rem, brig, bren, myID))
                i=i+1
            myID=myID+1
            