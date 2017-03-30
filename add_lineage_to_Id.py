import os
import fileinput
import sys

#Read the csv file, parse the headers
def add_lineage(path,path1):
    f=open(path,"r").read().split('\n')
    F=f[1:]
    Exitref_Id={}
    for item in F:
        Id=item and item.split(',')[]#add index of sample_id column
        lineage=item and item.split(',')[]#add index of lineage column 
        ID=Id.replace('-','')
        QC=item and item.split(',')[]#add index of QC column  
        if QC!='FAIL' and Id!='':
                Exitref_Id[ID]=lineage
    for key,value in Exitref_Id.iteritems():
        for lines in fileinput.input(path1,inplace=True):
            if key in lines:
                lines=lines.replace(key,key+'-'+value)
            sys.stdout.write(lines)
    
            
path= #path to csv file with lineages and samples id
path1=#path to fasta file or tree file     
Final=add_lineage(path,path1)
print Final
