import os
import random
label_Color={}
Lab=[]
column_heatmap=[]#enter the index of heatmap column separated by ,
column_label=[]#enter the index of non heatmap column separated by ,
f=open("/Users/sondeskalboussi/Desktop/rootedtrees/data1_binary.csv","r").read().split('\n')#add path to csv file between ""
F=[i.strip() for i in f]#remove return at the end of line
header=F[0].split(",")#list of header
for x in F[1:]:
    sample=x and x.split(',')[1]
    for i in range(2,len(header)):
        color = '#{:02x}{:02x}{:02x}'.format(*map(lambda x: random.randint(0, 255), range(3)))#generate random hex color
        if i not in label_Color and color not in label_Color.values():
            label_Color[header[i]]=color
        label=x and x.split(',')[i]
        if label is not "":
            new=open("script_"+str(i)+'.txt','a')
            script1="""
{} 1                        """.format(sample)
            new.write(script1)
            new.close()
for i in range(2,len(header)):
    #if os.path.exists("script_"+str(i)):
                with open("script_"+str(i)+'.txt', "r+") as f:
                        s = f.read()
                        f.seek(0,0)
                        f.write("%s\n%s\n%s\n%s%s\n%s\n%s\n%s%s\n%s" % ('DATASET_BINARY','SEPARATOR SPACE','DATASET_LABEL SNP','COLOR ',label_Color[header[i]],'FIELD_LABELS F1','FIELD_SHAPES 2','FIELD_COLORS ',label_Color[header[i]],'DATA')+ s)
