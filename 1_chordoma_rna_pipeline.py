#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2019 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run for processing of aligned RNA-seq data from UMCHOR1 and CH22 cells
#Raw reads were aligned to the HG19_ercc genome using Hisat2




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess

#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'chordoma_rna'
genome ='HG19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files


#genomeDirectory #select your genome
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/'
#genomeDirectory = '/grail/genomes/Mus_musculus/UCSC/hg19/Sequence/Chromosomes/'

genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
gtfFile = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes_ercc.gtf'



#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

ch22_rna_data_file = '%sdata_tables/190612_CH22_RNA_TABLE.txt' % (projectFolder)

umchor1_rna_data_file = '%sdata_tables/180221_UMCHOR1_RNA_TABLE.txt' % (projectFolder)


#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#=======================I. FIXING LINKS FOR BAMS=======================')
    print('#======================================================================')
    print('\n\n')


    # bam_folder = '/storage/cylin/grail/projects/chordoma_rna/190612_rna_seq/bams/'

    # def symlink_bai(bam_folder):

    #     '''
    #     resolves the symlinks of the bams to also symlink the bais
    #     '''
    #     bam_file_list = ['%s%s' % (bam_folder,fh) for fh in os.listdir(bam_folder) if fh.count('bam') > 0]
    #     print(bam_file_list)
    #     for bam_path in bam_file_list:
            
    #         #print(bam_path)
    #         #print(os.path.realpath(bam_path))
    #         bam_origin = os.path.realpath(bam_path)
            
    #         sym_origin = bam_origin.replace('.bam','.bam.bai')
    #         sym_dest = bam_path.replace('.bam','.bam.bai')
    #         #print(sym_origin)
    #         #print(sym_dest)
            
    #         sym_cmd ='ln -s %s %s' % (sym_origin,sym_dest)
    #         os.system(sym_cmd)


        

    # symlink_bai(bam_folder)

    print('\n\n')
    print('#======================================================================')
    print('#=====================II. LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for ch22 data file
    pipeline_dfci.summary(ch22_rna_data_file)

    #for umchor1
    pipeline_dfci.summary(umchor1_rna_data_file)


    print('\n\n')
    print('#======================================================================')
    print('#======================III. RUNNING CUFFNORM===========================')
    print('#======================================================================')
    print('\n\n')

    
    #for CH22 we want to do the cuff normalization in 3 batches

    # analysis_name = 'CH22_MERB_THZ1'
    # cufflinksFolder = utils.formatFolder('%scufflinks_ch22' % (projectFolder),True)
    # group_list = [['CH22_Merb_0h_1','CH22_Merb_0h_2','CH22_Merb_0h_3'],['CH22_Merb_THZ1_8h_1','CH22_Merb_THZ1_8h_2','CH22_Merb_THZ1_8h_3'],['CH22_Merb_THZ1_24h_1','CH22_Merb_THZ1_24h_2','CH22_Merb_THZ1_24h_3']]
    # bash_path = '%s%s_cuffquant.sh' % (cufflinksFolder,analysis_name)
    # pipeline_dfci.makeCuffTable(ch22_rna_data_file,analysis_name,gtfFile,cufflinksFolder,group_list,bash_path,useERCC = True)


    # analysis_name = 'CH22_MERB_DEG'
    # cufflinksFolder = utils.formatFolder('%scufflinks_ch22' % (projectFolder),True)
    # group_list = [['CH22_Merb_0h_1','CH22_Merb_0h_2','CH22_Merb_0h_3'],['CH22_Merb_Deg_8h_1','CH22_Merb_Deg_8h_2','CH22_Merb_Deg_8h_3'],['CH22_Merb_Deg_24h_1','CH22_Merb_Deg_24h_2','CH22_Merb_Deg_24h_3']]
    # bash_path = '%s%s_cuffquant.sh' % (cufflinksFolder,analysis_name)
    # pipeline_dfci.makeCuffTable(ch22_rna_data_file,analysis_name,gtfFile,cufflinksFolder,group_list,bash_path,useERCC = True)


    # analysis_name = 'CH22_WT_THZ1'
    # cufflinksFolder = utils.formatFolder('%scufflinks_ch22' % (projectFolder),True)
    # group_list = [['CH22_WT_THZ1_0h_1','CH22_WT_THZ1_0h_2','CH22_WT_THZ1_0h_3'],['CH22_WT_THZ1_8h_1','CH22_WT_THZ1_8h_2','CH22_WT_THZ1_8h_3'],['CH22_WT_THZ1_24h_1','CH22_WT_THZ1_24h_2','CH22_WT_THZ1_24h_3']]
    # bash_path = '%s%s_cuffquant.sh' % (cufflinksFolder,analysis_name)
    # pipeline_dfci.makeCuffTable(ch22_rna_data_file,analysis_name,gtfFile,cufflinksFolder,group_list,bash_path,useERCC = True)


    # analysis_name = 'CH22_WT_THZ1'
    # cufflinksFolder = utils.formatFolder('%scufflinks_ch22' % (projectFolder),True)
    # group_list = [['CH22_WT_THZ1_0h_1','CH22_WT_THZ1_0h_2','CH22_WT_THZ1_0h_3'],['CH22_WT_THZ1_8h_1','CH22_WT_THZ1_8h_2','CH22_WT_THZ1_8h_3'],['CH22_WT_THZ1_24h_1','CH22_WT_THZ1_24h_2','CH22_WT_THZ1_24h_3']]
    # bash_path = '%s%s_cuffquant.sh' % (cufflinksFolder,analysis_name)
    # pipeline_dfci.makeCuffTable(ch22_rna_data_file,analysis_name,gtfFile,cufflinksFolder,group_list,bash_path,useERCC = True)

    # analysis_name = 'UMCHOR1_THZ1'
    # cufflinksFolder = utils.formatFolder('%scufflinks_ch22' % (projectFolder),True)
    # group_list = [['UMCHOR1_THZ1_0h_1','UMCHOR1_THZ1_0h_2','UMCHOR1_THZ1_0h_3'],['UMCHOR1_THZ1_4h_1','UMCHOR1_THZ1_4h_2','UMCHOR1_THZ1_4h_3'],['UMCHOR1_THZ1_8h_1','UMCHOR1_THZ1_8h_2','UMCHOR1_THZ1_8h_3'],['UMCHOR1_THZ1_12h_1','UMCHOR1_THZ1_12h_2','UMCHOR1_THZ1_12h_3'],['UMCHOR1_THZ1_24h_1','UMCHOR1_THZ1_24h_2','UMCHOR1_THZ1_24h_3']]
    # bash_path = '%s%s_cuffquant.sh' % (cufflinksFolder,analysis_name)
    # pipeline_dfci.makeCuffTable(umchor1_rna_data_file,analysis_name,gtfFile,cufflinksFolder,group_list,bash_path,useERCC = True)



#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
