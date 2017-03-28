import os
import random
label_Color={}
path=#add path to csv file between ""
f=open(path,"r").read().split('\n')
header=f[0].split(",")[6]
F=[i.strip() for i in f[1:]]#remove the return \r at the end of list
Lab=[]
column_heatmap=[]#enter the index of heatmap column separated by ,
column_label=[]#enter the index of non heatmap column separated by ,
for x in F:
    sample=x and x.split(',')[0]
    for c in column_label:
        new=open("script_"+str(c),'a')
        label=x and x.split(',')[c]
        color = '#{:02x}{:02x}{:02x}'.format(*map(lambda x: random.randint(0, 255), range(3)))#generate random hex color
        if label not in label_Color and color not in label_Color.values():
            label_Color[label]=color
        if label is not "":
            script="""
{}	text={},fontcolor=black,linewidth=4,bkcolor={}""".format(sample,label,label_Color[label])
            new.write(script)
        new.close()
    
    for h in column_heatmap:
        new=open("script_"+str(h),'a')
        label=x and x.split(',')[h]
        if label is not "":
            Lab.append(float(label))
            group="%s\t%s" % (sample,label)
            Min=min(Lab)
            Max=max(Lab)
            L=len(Lab)
            Avg=sum(Lab)/L
            Avg1=sum(Lab[:(L/2)+1])/L/2
            Avg2=sum(Lab[(L/2)+1:])/L/2
            script="""
{}
        """.format(group)
            new.write(script)
        new.close()

for c in column_label:
    with open("script_"+str(c), "r+") as f: s = f.read(); f.seek(0); f.write("%s\n%s" % ('!grouplabel	style=2,color=pink,show=1,marginPCT=0.05,fontsize=14,fontcolor=white,fontitalic=0,textalign=middle,textorientation=horizontal,linewidth=2','!op	0.8')  + s)

for h in column_heatmap:
    with open("script_"+str(h), "r+") as f: s = f.read(); f.seek(0); f.write("%s\n%s%s%s%s%s%s%s%s%s%s\n%s\n%s"%('!colorgradient	green, yellow,red','!colorgradientMarkLabel	',Min,",",Avg1,",",Avg,",",Avg2,",",Max,
                                                                                                               '!heatmap	margin=2,colwidth=30,roundedcorne=0',

                                                                                                               '!showdataValue	show=1,fontsize=12')+ s)

