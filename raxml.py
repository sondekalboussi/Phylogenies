
import os
import sys
import fileinput
from Bio import SeqIO
path= # add path to directory containing all the fasta files
alignment=# create a new file for the alignment
raxml_input_file=open(alignment,'a')
#check length of fasta sequences than concatenate all the fastas in one file (raxml_input_file)
for f in os.listdir(path):
        os.chdir(path)
        FastaFile = open(f, 'rU')
        for rec in SeqIO.parse(FastaFile, 'fasta'):
            name = rec.id
            seq = rec.seq
            seqLen = len(rec)
            print f,seqLen
        FastaFile.close()  
for f in os.listdir(path):
    os.chdir(path)
    x=os.path.splitext(f)[0]
    id=x.split('_')[0]
    for lines in fileinput.input(path+'/'+f,inplace=True):
        if ("gi|444893469|emb|AL123456.3|") in lines:
            lines=lines.replace("gi|444893469|emb|AL123456.3|",id)
        sys.stdout.write(lines)
    for lines in open(path+'/'+f).readlines():
            raxml_input_file.write("%s\n" %(lines.strip()))
raxml_input_file.close()    
from datetime import datetime
start_time = datetime.now()
inp_fasta='/Users/sondeskalboussi/Desktop/DATA/alignment1.fa'
out_name='T1'
os.system("raxmlHPC -s {inp} -n {out} -m GTRCAT -p12345".format(inp=inp_fasta,out=out_name))
end_time = datetime.now()
print('Duration: {}'.format((end_time - start_time)))

