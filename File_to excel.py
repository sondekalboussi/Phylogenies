
###########Script1: fastq path in txt file
import pandas as pd
data={"File":[""],"Machine.Lane (ID:)":[""],"Technology (PL:)":[""],"Sample (SM:)":[""],"Lib (LB:)":[""]}
fl=open("/Users/sondeskalboussi/Desktop/File.txt").readlines()#txt file contains all the fastq files paths
for lines in fl:
        #lines=lines.strip()
        file=lines
        L=lines.split("_")
        name=L[3]
        lib=L[7]
        machine="vxfgghd"
        technology="illumina"
        data["File"]=file
        data["Machine.Lane (ID:)"]=machine
        data["Technology (PL:)"]=technology
        data["Sample (SM:)"]=name
        data["Lib (LB:)"]=lib
df=pd.DataFrame([data],columns=data.keys())#create dataframe from dictionary
df.to_excel('test.xlsx', sheet_name='sheet1', index=False)#dataframe to excel sheet
df.to_csv('example.csv',index=False, encoding='utf-8',sep="\t")#dataframe to csv 


################ script2: fastq files directory path###########
import pandas as pd
data={"File":[""],"Machine.Lane (ID:)":[""],"Technology (PL:)":[""],"Sample (SM:)":[""],"Lib (LB:)":[""]}
path= #fastq files directory path
for fl in os.listdirs(path):
    os.chdir(path) 
for lines in os.join.path(path,fl):
        #lines=lines.strip()
        file=lines
        L=lines.split("_")
        name=L[3]#index
        lib=L[7]
        machine="vxfgghd"
        technology="illumina"
        data["File"]=file
        data["Machine.Lane (ID:)"]=machine
        data["Technology (PL:)"]=technology
        data["Sample (SM:)"]=name
        data["Lib (LB:)"]=lib
df=pd.DataFrame([data],columns=data.keys())#create dataframe from dictionary
df.to_excel('test.xlsx', sheet_name='sheet1', index=False)#dataframe to excel sheet
df.to_csv('example.csv',index=False, encoding='utf-8',sep="\t")#dataframe to csv 

