import os
import random
file=#path to the metadata
label_Color={}
Lab=[]
column_heatmap=[]#enter the index of heatmap column separated by ,
column_label=[]#enter the index of non heatmap column separated by ,
f=open(file,"r").read().split('\n')#add path to csv file between ""
F=[i.strip() for i in f]#remove return at the end of line
header=F[0].split(",")#list of header
for x in F[1:]:
    sample=x and x.split(',')[0]
    for i in range(1,len(header)):
        label=x and x.split(',')[i]
        color = '#{:02x}{:02x}{:02x}'.format(*map(lambda x: random.randint(0, 255), range(3)))#generate random hex color
        if label not in label_Color and color not in label_Color.values():
            label_Color[label]=color  
        if label is not "":
            try: 
                float(label)
                if float(label)%1!=0:
                        new=open("script_heatmap_"+header[i],'a')
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
                else: 
                    new=open("script_label_"+header[i],'a') 
                    script1="""
{}	text={},fontcolor=black,linewidth=4,bkcolor{}""".format(sample,label,label_Color[label])
                    new.write(script1)
                    new.close()
            except:        
                    new=open("script_label_"+header[i],'a') 
                    script1="""
{}	text={},fontcolor=black,linewidth=4,bkcolor{}""".format(sample,label,label_Color[label])
                    new.write(script1)
                    new.close()
for i in range(1,len(header)):
    if os.path.exists("script_label_"+header[i]):
        with open("script_label_"+header[i], "r+") as f: s = f.read(); f.seek(0); f.write("%s\n%s" % ('!grouplabel	style=2,color=pink,show=1,marginPCT=0.05,fontsize=14,fontcolor=white,fontitalic=0,textalign=middle,textorientation=horizontal,linewidth=2','!op	0.8')  + s)
    if os.path.exists("script_heatmap_"+header[i]):
        with open("script_heatmap_"+header[i], "r+") as f: s = f.read(); f.seek(0); f.write("%s\n%s%s%s%s%s%s%s%s%s%s\n%s\n%s"%('!colorgradient	green, yellow,red','!colorgradientMarkLabel	',Min,",",Avg1,",",Avg,",",Avg2,",",Max,
                                                                                                                             '!heatmap	margin=2,colwidth=30,roundedcorne=0',
                                                                                                                             
                                                                                                                             '!showdataValue	show=1,fontsize=12')+ s)

