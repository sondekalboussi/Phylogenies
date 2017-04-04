import os
import fileinput
import sys
def add_lineage(path,path1):
    #Read the csv file, parse the headers
    f=open(path,"r").read().split('\n')
    F=f[1:]
    Exitref_Id={}
    for item in F:
        Id=item and item.split(',')[0]#IndexError problem
        lineage=item and item.split(',')[]#add column index to lineage
        ID=Id.replace('-','')
        QC=item and item.split(',')[]#add column index to QC
        if QC!='FAIL' and Id!='':
            Exitref_Id[]=lineage # if the sample name in the tree/fasta contains a - use Exitref_Id[Id] otherwise use Exitref_Id[ID]
    for key,value in Exitref_Id.iteritems():
        for lines in fileinput.input(path1,inplace=True):
            if key in lines:
                lines=lines.replace(key,key+'-'+value)
            sys.stdout.write(lines)
            
path= #path to csv file with lineages and samples id
path1=#path to fasta file or tree file     
Final=add_lineage(path,path1)
print Final


