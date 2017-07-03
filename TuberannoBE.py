"""
Save all the fasta files of the reference genomes together in one directory.
If you want to run only one mapper please change the Boolean in "pipeline" into False for the unused mapper
when you enter the ref_genome fasta path dont put the extension
Please make sure to install java and perl 
converting the gtf file into genepred file must be run on linux OS (no mac/windows version)
"""
import os
import sys
from itertools import product
from datetime import datetime
start_time = datetime.now()
new=open(path to txt file for time recording,"w")

class mapping_to_ref_genome_functional_annotation(object):
    def __init__(self,Novoalign,BWA_mem,BWA_mem_map_output,novo_align_map_output,ref_genome,ref_gen_Dir,Trimmomatic,bwa,novoalign,GATK,fastq_Dir,dbSNP,illumina_adapters,samtools,picard,bedtools,bcftools,Delly,annovar,gtfToGenePred,MTdb)
      sel.Novoalign=Novoalign
      self.BWA_mem=BWA_mem
      self.BWA_mem_map_output=BWA_mem_map_output
      self.novo_align_map_output=novo_align_map_output
      self.ref_genome=ref_genome
      self.ref_gen_Dir=ref_gen_Dir
      self.trimmomatic=Trimmomatic
      self.bwa=bwa
      self.novoalign=novoalign
      self.GATK=GATk
      self.fastq_Dir=fastq_Dir
      self.dbSNp=dbSNP
      self.illumina_adapters=illumina_adapters
      self.samtools=samtools
      self.picard=picard
      self.bedtools=bedtools
      self.bcftools=bcftools
      self.Delly=Delly
      self.annovar=annovar
      self.gtfToGenePred=gtfToGenePred
      self.annotation_db=MTdb
    
    def process_ref_genome(self):
        #prepare the ref genome to be used for the mapping and GATK
            for self.ref_genome in os.listdir(self.ref_gen_Dir+"/"+self.genome):
                os.chdir(ref_gen_Dir+"/"+genome)
            try:
                    os.system("{}/bwa index {}.fasta""".format(self.bwa,self.ref_genome))
                    os.system("{}/novoindex {}.nix {}.fasta "format(self.novoalign,self.ref_genome,self.ref_genome))         
                    os.system("{}}/samtools faidx {}""".format(self.samtools,self.ref_genome))
                    os.system("java -jar {}/picard.jar CreateSequenceDictionary R={}.fasta O={}.dict""".format(self.picard,self.ref_genome,self.ref_genome))
            except ValueError:
                    print "Error, check the genome fasta file and/or that all the software are executable"
            
            return "Ready to start mapping!"  
    
    def trimming(self):
        global ID
        if self.BWA_mem:
            if not os.path.isdir(self.BWA_mem_map_output+"/Trimming"):
                    os.makedirs(self.BWA_mem_map_output+"/Trimming")  
            for fl in os.listdir(self.fastq_Dir):
                os.chdir(self.fastq_Dir)
                id=fl.rpartition(".fastq.gz ")[0]
                ID=id.split("_")[0]
                if "R1"in id:
                        read1=fastq_Dir+"/"+ID+"_R1.fastq.gz"
                if "R2"in id:    
                        read2=fastq_Dir+"/"+ID+"_R2.fastq.gz"
                output1=self.BWA_mem_map_output+"/Trimming/"+id+"paired_trimed.fastq")
                output2=self.BWA_mem_map_output+"/Trimming/"+id+"unpaired_trimed.fastq")
                output3=self.BWA_mem_map_output+"/Trimming/"+id+"paired_trimed.fastq")
                output=self.BWA_mem_map_output+"/Trimming/"+id+"unpaired_trimed.fastq") 
                os.system("""java -jar {}/trimmomatic.jar PE -threads 8 -phred33 {} {} {} {} {} {} ILLUMINACLIP:{}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36""".format(Trimmomatic,read1,read2,illumina_adapters,output1,output2,output3,output4))
        if self.Novoalign:   
            if not os.path.isdir(self.novoalign_mem_map_output+"/Trimming"):
                    os.makedirs(self.novoalign_mem_map_output+"/Trimming")  
            for fl in os.listdir(self.fastq_Dir):
                os.chdir(self.fastq_Dir)
                if "R1"in id:
                        read1=fastq_Dir+"/"+ID+"_R1.fastq.gz"
                if "R2"in id:    
                        read2=fastq_Dir+"/"+ID+"_R2.fastq.gz"
                output1=self.novoalign_mem_map_output+"/Trimming/"+id+"paired_trimed.fastq")
                output2=self.novoalign_mem_map_output+"/Trimming/"+id+"unpaired_trimed.fastq")
                output3=self.novoalign_mem_map_output+"/Trimming/"+id+"paired_trimed.fastq")
                output=self.novoalign_mem_map_output+"/Trimming/"+id+"unpaired_trimed.fastq") 
                os.system("""java -jar {}/trimmomatic.jar PE -threads 8 -phred33 {} {} {} {} {} {} ILLUMINACLIP:{}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36""".format(Trimmomatic,read1,read2,illumina_adapters,output1,output2,output3,output4))
        return "Trimming is done!"

       
    def run_mapping(self): 
        #run the mapping using Bwa-mem and Novoalign,the final output is sorted and indexed bam              
        if self.BWA_mem:
            if not os.path.isdir(self.BWA_mem_map_output+"/Alignment"):
                os.makedirs(self.BWA_mem_map_output+"/Alignment")           
                print "Running BWA-mem..."
            for fl in os.listdir(self.BWA_mem_map_output+"/Trimming):
                os.chdir(self.BWA_mem_map_output+"/Trimming)
                if "R1"in id:
                    read1=self.BWA_mem_map_output+"/Trimming/"+ID+"trim_R1.fastq.gz"
                if "R2"in id:    
                    read2=self.BWA_mem_map_output+"/Trimming/"+ID+"trim_R2.fastq.gz"
                readGroup="'@RG\\tID:"+IDPE[pos]+"\\tSM:"+SMPE[pos]+"\\tLB:"+LBPE[pos]+"\\tPL:Illumina'"         
                output=BWA_mem_map_output+"/Alignment/"+ID+"_sorted_bwa.bam")
                os.system("{}/bwa mem -M -t 16 -R {} {}.fasta {} {} | {}/samtools view -Sb - | {}/samtool sort - | {}/samtools index - > {}""".format(self.bwa,self.ref_genome,readGroup,read1,read2,self.samtools,self.samtools,self.samtools,output))
        if self.Novoalign:
            if not os.path.isdir(self.novo_align_map_output+"/Alignment"):
                os.makedirs(self.novo_align__map_output+"/Alignment")
                print "Running Novoalign..."
            for fl in os.listdir(self.novoalign_mem_map_output+"/Trimming):
                os.chdir(self.novoalign_mem_map_output+"/Trimming)
                id=fl.rpartition(".fastq.gz ")[0]
                ID=id.split("_")[0]
                if "R1"in id:
                    read1=self.BWA_mem_map_output+"/Trimming/"+ID+"trim_R1.fastq.gz"
                if "R2"in id:    
                    read2=self.BWA_mem_map_output+"/Trimming/"+ID+"trim_R2.fastq.gz"
                output=novo_align_map_output+"/Alignment/"+ID+"_sorted_novoalg_.bam"
                os.system("""{}/novoalign -k -d {}.nix -f {} {} -a -o SAM {} 2> {}_stats.txt | {}/samtools -bS - | {}/samtool sort - | {}/samtools index - > {}""".format(novoalign,ref_genome,readGroup,read1,read2,samtools,samtools,samtools,output))                        .
        print BWA_mem_map_output+"/Alignment"
        print novo_align__map_output+"/Alignment"
        return "mapping is done..."          #run bwa-mem and novolaign for mapping
    
    def realignment(self):
        #realign around the indels usng GATK the final output is a sorted indexed bam file
          if self.BWA_mem:
              for fl in os.listdir(self.BWA_mem_map_output):
                  os.chdir(self.BWA_mem_map_output)
                  input=self.BWA_mem_map_output+"/Alignment/"+ID+"sorted_bwa.bam"
                  output1=self.BWA_mem_map_output+"/Alignment/"+ID+"_bwa.intervals"
                  output2=self.BWA_mem_map_output+"/Alignment/"+ID+"_realign_bwa.bam"
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {}.fasta -I {} -o {}""".format(self.GATK,self.ref_genome,input,output1)
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T IndelRealigner -R {}.fasta -I {} -o {} -targetIntervals {}""".format(self.ref_genome,input,output2,input1))
                  os.system("""{}/samtools index {}""".format(self.samtools,output2))
          if self.Novoalign:
              for fl in os.listdir(self.novo_align_map_output):
                  os.chdir(self.novo_align_map_output)
                  id=os.path.splitext(fl)[0]
                  ID=id.split("_")[0]
                  input=self.novo_align_map_output+"/Alignment/"+ID+"sorted_novoalg.bam"
                  output1=self.novo_align_map_output+"/Alignment/"+ID+"_novoalg.intervals"
                  output2=self.novo_align_map_output+"/Alignment/"+ID+"_realign_novoalg.bam"
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {}.fasta -I {} -o {}""".format(self.GATK,self.ref_genome,input,output1)
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T IndelRealigner -R {}.fasta -I {} -o {} -targetIntervals {}""".format(self.GATK,self.ref_genome,input,output2,input1))
                  os.system("""{}/samtools index {}""".format(self.samtools,output2))
                
         return "realignment around indels is done!"        
    
    def base_qual_recal():
        #base quality score recalibration
          if self.BWA_mem:
              for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                  os.chdir(self.BWA_mem_map_output+"/Alignment/")
                  id=os.path.splitext(fl)[0]
                  ID=id.split("_")[0]
                  input=self.BWA_mem_map_output+"/Alignment/"+ID+"_realign_bwa.bam"
                  output1=self.BWA_mem_map_output+"/Alignment/"+ID+"bwa_recal_data.table"
                  output2=self.BWA_map_output+"/Alignment/"+ID+"_realign_recal_bwa.bam"
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T BaseRecalibrator -R {}.fasta -I {} -knownSites {} -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output1))
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T PrintReads -R {}.fasta -BQSR {} -I sample1.sorted.realg.bam -o {}""".format(self.GATK,self.ref_genome,input,output1,output2))
                  os.system("""{}/samtools index {}""".format(self.samtools,output2))
          if self.Novoalign:
              for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                  os.chdir(self.novo_align_map_output+"/Alignment/")
                  input=self.novo_align_map_output+"/Alignment/"+ID+"_realig_novoalg.bam"
                  output1=self.novo_align_map_output+"/Alignment/"+ID+"_novoalg_recal_data.table"
                  output2=self.novo_align_map_output+"/Alignment/"+ID+"_realign_recal_novoalg.bam"   
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T BaseRecalibrator -R {}.fasta -I {} -knownSites {} -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output1))
                  os.system("""java -jar {}/GenomeAnalysisTK.jar -T PrintReads -R {}.fasta -BQSR {} -I sample1.sorted.realg.bam -o {}""".format(self.GATK,self.ref_genome,input,output1,output2))
                  os.system("""{}/samtools index {}""".format(self.samtools,output2))
    
         return "Base quality score recalibration is done!""
                     
    def PCR_dup_remov(self):
        #mark and remove the duplicates
            if self.BWA_mem:
                for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                    os.chdir(self.BWA_mem_map_output+"/Alignment/")
                    input=self.BWA_mem_map_output+"/Alignment/"+ID+"_realign_bwa.bam"
                    output=self.BWA_map_output+"/Alignment/"+ID+"_realign_dedup_bwa.bam"
                    os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} OUTPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true""".format(self.picard,input,output))
                    os.system("""{}/samtools index {}""".format(self.samtools,output))
            if self.Novoalign:
                for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                    os.chdir(self.novo_align_map_output+"/Alignment/")
                    input=self.novo_align_map_output+"/Alignment/"+ID+"_realig_novoalg.bam"
                    output=self.novo_align_map_output+"/Alignment/"+ID+"_realign_dedup_novoalg.bam"
                    os.system("""java -jar {}/picard.jar MarkDuplicates INPUT={} OUTPUT={} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT M=duplicate_metrics T MP_DIR=tmp ASSUME_SORTED=true""".format(self.picard,input,output))
                    os.system("""{}/samtools index {}""".format(self.samtools,output2))
            
            return "PCR duplicates removal is done!"
    
    def merge_bam(self):
        # merge bam files of same sample from different library add read group using picard and sort and index the merged bam
            if self.BWA_mem:
                if not os.path.isdir(self.BWA_mem_map_output+"/Merged_alignment"):
                    os.makedirs(self.BWA_mem_map_output+"/Merged_alignment")   
                for fl in os.listdir(self.BWA_mem_map_output+"/Alignment/"):
                    os.chdir(self.BWA_mem_map_output+"/Alignment/")
                    multi_library=[ ]
                    to_be_merged=""
                    if ID in fl:
                        multi_library.append(self.BWA_mem_map_output+"/Alignment/"+fl)
                    #when i find the id in the file name write the file name multi
                    for i in multi_library:
                        to_be_merged+=i+" "
                    output=self.BWA_map_output+"/Merged_alignment/"+ID+"_merged_bwa.bam"
                    cmd="""{}/samtools merge {} | {}/samtool sort - | {}/samtools index - > {}""".format(self.samtools,output,to_be_merged)
                    os.system(cmd)
            if self.Novoalign:
                if not os.path.isdir(self.novo_align_map_output+"/Merged_alignment"):
                    os.makedirs(self.novo_align_map_output+"/Merged_alignment")   
                for fl in os.listdir(self.novo_align_map_output+"/Alignment/"):
                    os.chdir(self.novo_align_map_output+"/Alignment/")
                    multi_library=[ ]
                    to_be_merged=""
                    if ID in fl:
                        multi_library.appendself.(novo_align_map_output+"/Alignment/"+fl)
                    #when i find the id in the file name write the file name multi
                    for i in multi_library:
                        to_be_merged+=i+" "
                    output=self.novo_align_map_output+"/Merged_alignment/"+ID+"_merged_bwa.bam"
                    cmd="""{}/samtools merge {} | {}/samtool sort - | {}/samtools index - > {}""".format(self.samtools,output,to_be_merged)
                    os.system(cmd)        
            print  self.BWA_mem_map_output+"/Merged_alignment"
            print  self.novo_align_map_output+"/Merged_alignment"
            return "merging bam is done!" 

    def genome_coverag_mappability_stat(self.):
            if self.BWA_mem:
                genom_cov={}#sample with genome coverage >40
                read_map={}#sample with mapped reads >90
                if not os.path.isdir(self.BWA_mem_map_output+"/statistics/final_stat"):
                    os.makedirs(self.BWA_mem_map_output+"/statistics/final_stat")   
                for fl in os.listdir(self.BWA_mem_map_output+"/Merged_alignment/"):
                    os.chdir(self.BWA_mem_map_output+"/Merged_alignment/")
                    input=self.BWA_mem_map_output+"/Merged_alignment/"+ID+"_merged_bwa.bam"
                    output1=self.BWA_mem_map_output+"/statistics/"+ID+"_bwa_genomecov_bed.txt"
                    #output2=BWA_mem_map_output+"/statistics/"+ID+"_bwa_genomecov>40.txt"
                    output2=self.BWA_mem_map_output+"/statistics/"+ID+"_bwa_genomecov=0.txt"
                    os.system("""{}/bedtools genomecov -ibam -bga {} -g {}.fasta > {}""".format(self.bedtools,input,self.ref_genome,output1)
                    cov=os.system("awk '{ print $4 }' "+output1)
                    genom_cov[ID]=cov          
                    #os.system("""awk 'NF && $4>40' {} > {}""".format(output1,output2))
                    os.system("""awk 'NF && $4<2' {} > {}""".format(output1,output3))
                    output11=self.BWA_mem_map_output+"/statistics/"+ID+"_stats_bwa.txt"
                    os.system("""{}/samtools flagstat {} > {}""".format(self.samtools,input,output11))
                for fl in os.listdir(self.BWA_mem_map_output+"/statistics/"): 
                    if fl.endwith("_stats_bwa.txt"):          
                        for line in open(fl).readlines():
                            line=line.strip()
                            if "mapped (" in line:
                                    x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                    y=float(x)
                                    read_map[ID]=y                                                          
                for (i, j), (k, v) in product(genom_cov.items(), read_map.items()):   
                              if i==k and j>=90 or v>=40:
                                  os.system("""cp {} {}""".format(self.BWA_mem_map_output+"/statistics/"+ID+"_stats_bwa.txt",self.BWA_mem_map_output+"/statistics/final_stat"+ID+"_stats_bwa.txt"))
                                  
            if self.Novoalign:
                genom_cov={}
                read_map={}
                if not os.path.isdir(self.novo_align_map_output+"/statistics/final_stat"):
                    os.makedirs(self.novo_align_map_output+"/statistics/final_stat")   
                for fl in os.listdir(self.novo_align_map_output+"/Merged_alignment/"):
                    os.chdir(self.novo_align_map_output+"/Merged_alignment/")
                    id=os.path.splitext(fl)[0]
                    ID=id.split("_")[0]
                    input=self.novo_align_map_output+"/Merged_alignment/"+ID+"_merged_novo.bam"                     
                    output1=self.novo_align_map_output+"/statistics/"+ID+"_novo_genomecov_bed.txt"
                    #output2=novo_align_map_output+"/statistics/"+ID+"_novo_genomecov>40.txt"#genome coverage 40 or more
                    output2=self.novo_align_map_output+"/statistics/"+ID+"_novo_genomecov=0.txt"#genome coverage equal to zero
                    os.system("""{}/bedtools genomecov -ibam -bga {} -g {}.fasta > {}""".format(self.bedtools,input,self.ref_genome,output1)
                    cov=os.system("awk '{ print $4>=40 }' "+output1)
                    genom_cov[ID]=cov          
                    #os.system("""awk 'NF && $4>=40' {} > {}""".format(output1,output2))
                    os.system("""awk 'NF && $4<2' {} > {}""".format(output1,output2))
                    output11=self.novo_align_output+"/statistics/"+ID+"_stats_novo.txt"
                    os.system("""{}/samtools flagstat {} > {}""".format(self.samtools,input,output11))
                for fl in os.listdir(self.novo_align_mem_map_output+"/statistics"):
                   if fl.endwith("_stats_novo.txt"):           
                    for line in open(fl).readlines():
                        line=line.strip()
                        if "mapped (" in line:
                                x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                y=float(x)
                                read_map[ID]=y
                              
                for (i, j), (k, v) in product(genom_cov.items(), read_map.items()):   
                              if i==k and j>=90 or v>=40:                    
                                  os.system("""cp {} {}""".format(self.novo_align_map_output+"/statistics/"+ID+"_stats_novo.txt",self.novo_align_map_output+"/statistics/final_stat/"+ID+"_stats_novo.txt"))                          
            print self.BWA_mem_map_output+"/statistics/final_stat"
            print self.novo_align_map_output+"/statistics/final_stat"
            return ""

    def joint_variant_calling(self):
               # joint variant call SNP and indels using GATK
                if self.BWA_mem:
                    if not os.path.isdir(self.BWA_mem_map_output+"/Joint_Variants"):
                            os.makedirs(self.BWA_mem_map_output+"/Joint_Variants")
                    for fl in os.listdir(self.BWA_mem_map_output+"/Merged_alignment/"):
                            os.chdir(self.BWA_mem_map_output+"/Merged_alignment/")
                            for id in sample:
                                input=self.BWA_mem_map_output+"/Fianl_alignment/"+id+"_merged_bwa.bam"
                                output=self.BWA_map_output+"Variants/"+id+"_GATK_snps_indels.vcf" 
                                os.system("""java -jar {}/GenomeAnalysisTK.jar -R {}.fasta -T HaplotypeCaller -I {} --dbsnp {} -stand_call_conf 30 -stand_emit_conf 10.0 -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output,))
            
                if self.Novoalign:
                    if not os.path.isdir(self.novo_align_map_output+"/Joint_Variants"):
                            os.makedirs(self.novo_align_map_output+"/Joint_Variants")
                    for fl in os.listdir(self.novo_align_map_output+"/Merged_alignment"):
                        os.chdir(self.novo_align_map_output+"/Merged_alignment")
                        for id in sample:
                            input=self.novo_align_map_output+"/Merged_alignment/"+id+"_merged_novo.bam"
                            output=self.novo_align_output+"/reads_stat/"+id+"_GATK_snps_indels.vcf"
                            os.system("""java -jar {}/GenomeAnalysisTK.jar -R {}.fasta -T HaplotypeCaller -I {} --dbsnp {} -stand_call_conf 30 -stand_emit_conf 10.0 -o {}""".format(self.GATK,self.ref_genome,input,self.dbSNP,output,))
                
                print self.BWA_mem_map_output+"/Variants"
                print self.novo_align_map_output+"/Variants"
                return "" 
        
    def structural_variation_calling(self):
                #structural variant (big deletion) using DELLY
                if self.BWA_mem:
                    if not os.path.isdir(self.BWA_mem_map_output+"/struc_Variants"):
                            os.makedirs(self.BWA_mem_map_output+"/struct_Variants")
                    for fl in os.listdir(self.BWA_mem_map_output+"/Merged_alignment/"):
                            os.chdir(self.BWA_mem_map_output+"/Merged_alignment/")
                            for id in sample:
                                input=self.BWA_mem_map_output+"/Fianl_alignment/"+id+"_merged_bwa.bam"#sorted indexed
                                output=self.BWA_mem_map_output+"/struc_Variants"+id+"_bwa_SV"
                                os.system("""{}/delly call -t DEL -g {} -o {}.bcf {}""".format(self.DElly,self.ref_genome,output,input)
                                os.system("""{}/bcftools call -mvO v -o {}.vcf {}.bcf""".format(self.bcftools,output,output))#convert bcf to vcf
                                
                if self.Novoalign:
                    if not os.path.isdir(self.novo_align_map_output+"/struc_Variants"):
                            os.makedirs(self.novo_align_map_output+"/struct_Variants")
                    for fl in os.listdir(self.novo_align_map_output+"/Merged_alignment/"):
                            os.chdir(self.novo_align_map_output+"/Merged_alignment/")
                            for id in sample:
                                input=self.novo_align_map_output+"/Fianl_alignment/"+id+"_merged_bwa.bam"#sorted indexed
                                output=self.novo_align_map_output+"/struc_Variants"+id+"_novo_SV"
                                os.system("""{}/delly call -t DEL -g {} -o {}.bcf {}""".format(self.DElly,self.ref_genome,output,input)
                                os.system("""{}/bcftools call -mvO v -o {}.vcf {}.bcf""".format(self.bcftools,output,output))#convert bcf to vcf                      
                print self.BWA_mem_map_output+"/struc_Variants"
                print self.novo_align_map_output+"/struc_Variants"
                return "" 
        
    def mapping_functional_annotation_annovar(self):
             #prepare annotation_files: convert gtf into genepred file
                for fl in os.listdir(self.annotation_db):
                     if not os.path.exists(self.annotation_db+"/H37RV_refGene.txt")                   
                        os.system("""{} -genePredExt {} """.format(self.gtfToGenePred,self.H37RV.gtf,self.annotation_db+"/H37RV_refGene.txt")                             
             #run annovar annotation():
                if self.BWA_mem:
                    if not os.path.isdir(self.BWA_mem_map_output+"/Annotation"):
                            os.makedirs(self.BWA_mem_map_output+"/Annotation")
                    for fl in os.listdir(self.BWA_mem_map_output+"/Joint_Variants/"):
                            os.chdir(self.BWA_mem_map_output+"/Joint_Variants/")
                            for id in sample:
                                input=self.BWA_mem_map_output+"/Joint_Variants/"+id+".vcf"
                                output=self.BWA_mem_map_output+"/Joint_Variants_annovar/"+id+".annovar"
                                os.system("""perl {}/table_annovar.pl --vcfinput {} {} --outfile {} --buildver H37RV --nastring . --protocol refGene,dbsnp --operation g,r,f  --vcfinput""".format(self.annovar,input,TBdb,output)
            
                if self.Novoalign:
                    if not os.path.isdir(self.novo_align_map_output+"/Annotation"):
                            os.makedirs(self.novo_align_map_output+"/Annotation")
                    for fl in os.listdir(self.novo_align_map_output+"/Joint_Variants/"):
                            os.chdir(self.novo_align_map_output+"/Joint_Variants/")
                            for id in sample:
                                input=self.novo_align_map_output+"/Joint_Variants/"+id+".vcf"
                                output=self.novo_align_map_output+"/Annotation/"+id+".annovar"
                                os.system("""perl {}/table_annovar.pl --vcfinput {} {} --outfile {} --buildver H37RV --nastring . --protocol refGene,dbsnp --operation g,r,f  --vcfinput""".format(self.annovar,input,TBdb,output)          
                                
                return""
                
if __name__ == '__main__':
    
    pipeline=mapping_to_ref_genome_functional_annotation(Novoalign,BWA_mem,BWA_mem_map_output,novo_align_map_output,ref_genome,ref_gen_Dir,Trimmomatic,bwa,novoalign,GATK,fastq_Dir,dbSNP,illumina_adapters,samtools,picard,bedtools,bcftools,Delly,annovar,gtfToGenePred,MTdb)                                      
    print pipeline.process_ref_genome()
    print pipeline.run_mapping()
    print pipeline.realignment
    print pipeline.base_qual_recal()
    print pipeline.PCR_dup_remov()
    print pipeline.merge_bam()
    print pipeline.genome_coverag_mappability_stat()
    print pipeline.joint_variant_calling()
    print pipeline.structural_variation_calling()                                      
    print pipeline.mapping_functional_annotation_annovar()

end_time = datetime.now()
new.write('Pipeline execution time: {}'.format(end_time - start_time))
new.close()

#perl annotate_variation.pl -geneanno input -buildver MTB(prefixe of files in TBdb) TBdb(where all the files)
#         perl annotate_variation.pl -filter -dbtype dbsnp -buildver MTB(prefixe of files in TBdb) input TBdb(where all the files)
#          perl annotate_variation.pl -regionanno -dbtype dbsnp -buildver MTB(prefixe of files in TBdb) input TBdb(where all the files)
 # os.system("""perl {}/convert2annovar.pl {} -format vcf4 > {}""".format(self.annovar,input,output))
#os.system("""perl {}/convert2annovar.pl {} -format vcf4 > {}""".format(self.annovar,input,output))
#samtools mpileup -uf ref.fa aln.bam | bcftools call -mv > var.raw.vcf bcftools filter -s LowQual -e '%QUAL<20 || DP>00' var.raw.vcf  > var.flt.vcf                                          
