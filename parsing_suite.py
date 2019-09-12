#!/bin/env python3

#List directories:

'''Useful: https://realpython.com/working-with-files-in-python/'''

import re
import os
import shutil
from config import config as cg
from clean_pdb import clean_pdb as cp
from residue_extract import interface_residues as ir
from file_management import file_management as fm
import csv


pdbpath='/home/yobi/AbPackingAngle_V2.1/LH_Combined_Martin/'

"""
Essential: Execute each Section separately, by copying and pasting each section
into a separate script and executing

Reasons: some sections require copying contents into BASH or light formatting
in spreadsheet software such as Libre Calculation/ Excel/ Google Sheets
"""
#=======================================
#1) Generating dataset of packing angles
#=======================================

#Records list of files in folder downloaded from AbDb -Antibody Database Martin 2018
#Reads contents of folder, writes list of folder contents to a text file - pdb_list.txt
w = open('pdb_list.txt', 'a+')
w.truncate(0)
fm.write_list(cp.list_directory(pdbpath), 'pdb_list.txt')

#Write the list of downloaded pdbs to a bash executable for Ab Packing Angle
#Reads pdb_list.txt, writes to - bash_prangle_list.txt
w = open('bash_prangle_list.txt', 'a+')
w.truncate(0)
fm.write_list(cp.bash4_pack('pdb_list.txt'),'bash_prangle_list.txt')

#Copy contents of bash_prangle_list.txt into BASH, in folder with Ab Packing Angle V2.1
#Execute


#================================================
#2) Extract Free vs Complex pdbs for Pangles
#========================================
#1) Download from AB PAcking Angle Database, "Free Antibodies/Complexes List", found here:
#http://www.bioinf.org.uk/abs/abdb/Data/Martin_logs/FreeAntibody_AntibodyAntigen.list
#renamed: list_all_freevscomplex.txt

#Open list_all_freevscomplex.txt in Libre Calculator, seperate cells by space:
# 1) copy first column to list_all_free.txt
# 2) copy second column to list_all_complex.txt

#Generate clean dataset for Free:
w = open('list_all_free_4extract.txt', 'a+')
w.truncate(0)
cp.clean4pdb_extraction('list_all_free.txt','list_all_free_4extract.txt')


#Generate clean dataset for Complexed
w = open('list_all_complex_4extract.txt', 'a+')
w.truncate(0)
cp.clean4pdb_extraction('list_all_complex.txt','list_all_complexed_4extract.txt')



#Free: Generate BASH executable code for free pdbs Ab Packing Angle:
w = open('free_4bash.txt', 'a+')
w.truncate(0)
fm.write_list(cp.bash4_pack('list_all_free_4extract.txt','free_pangles.txt'),'free_4bash.txt')

#Complexed: Generates BASH executable code for complexed pdbs Ab Packing Angle:
w = open('complex_4bash.txt', 'a+')
w.truncate(0)
fm.write_list(cp.bash4_pack('list_all_complexed_4extract.txt', 'complex_pangles.txt'),'complex_4bash.txt')


#NOW: Copy both lists into BASH and run in folder with Ab Packing Angle
#Creates two files
#1) Free = free_pangles.txt
#2) Complex = complex_pangles.txt




#===============================================
#3) Formatting All, Free and Complex for RStudio
#===============================================
#all_pdbs_pangles
#complex_pangles.txt
#free_pangles.txt


#Open text file, clean(garbled lines, approx 50)
#and write to a csvfile (works for Libre calc)
w = open('full_prangle_list.csv', 'a+')
w.truncate(0)
cp.prangleclean4_csv('/home/yobi/AbPackingAngle_V2.1/all_pdbs_pangles.txt','full_prangle_list.csv')


#converts txt file pdb to have commas for input to R and stat analysis
#copy csv pangle column to text file and add commas
w = open('free_commas4R.txt', 'a+')
w.truncate(0)
fm.write_list(cp.add_commas('free_nocommas.txt'),'free_commas4R.txt')

w = open('complex_commas4R.txt', 'a+')
w.truncate(0)
fm.write_list(cp.add_commas('complex_nocommas.txt'),'complex_commas4R.txt')


#converts all plain pdbs into '','' and identifies 5+ pdbs per line
w=open('fiveplus_pdbs.txt', 'a+')
w.truncate(0)
fm.write_list(cp.apost_commas_4pdbs('list_all_freevscomplex.txt'), 'fiveplus_pdbs.txt')



#===========================================================
#4) Generating scripts to derive mean, range, sd for RStudio
#===========================================================

#creates dataframe of all pdbs
df=fm.csv_2df('all_pdbs.csv') 
gen_4r.all_pangle_groupjoin()

#All
gen_4r.mean_gen4_allgroups()
gen_4r.range_gen4_allgroups()
gen_4r.sd_gen4_allgroups()

#Free
gen_4r.mean_gen4_freegroups()
gen_4r.range_gen4_freegroups()
gen_4r.sd_gen4_freegroups()

#Complex
gen_4r.mean_gen4_complexgroups()
gen_4r.range_gen4_complexgroups()
gen_4r.sd_gen4_complexgroups()

#Copy generated scripts from shell into R(Studio)

#====================================
#5) Dataframes for Pandas and Seaborn
#====================================
#Creates:
#free_g1=df_all_pangles.loc[['1MHH_1','1NLB_1','1YMH_1']]
#complex_g1=df_all_pangles.loc[['1A6V_1','1A6V_2','1A6W_1']]

#Stage 1: Replace comma duplicates with singles
w=open('fiveplus_pdbs_s1.txt', 'a+')
w.truncate(0)
fm.write_list(cp.replace_commas_withcomma('fiveplus_pdbs.txt'), 'fiveplus_pdbs_s1.txt')


#Stage 2: Remove comma that appears at the end of some lines
w=open('fiveplus_pdbs_s2.txt', 'a+')
w.truncate(0)
fm.write_list(cp.remove_commas_end('fiveplus_pdbs.txt'), 'fiveplus_pdbs_s2.txt')

#stage 3: add apostraphes around commas
w=open('fiveplus_pdbs_s3.txt', 'a+')
w.truncate(0)
fm.write_list(cp.apost_aroundcommas('fiveplus_pdbs_s1.txt'), 'fiveplus_pdbs_s3.txt')

#Stage 4: homogenise " into '
w=open('fiveplus_pdbs_s4.txt', 'a+')
w.truncate(0)
fm.write_list(cp.same_apostraphes('fiveplus_pdbs_s3.txt'), 'fiveplus_pdbs_s4.txt')

#Stage 5: Format for pandas dataframe
w=open('fiveplus_pdbs_s5.txt', 'a+')
w.truncate(0)
fm.write_list(cp.split_semicolon_4pdbs('fiveplus_pdbs_s4.txt'), 'fiveplus_pdbs_s5.txt')

#Stage 6: homogenise " into '
w=open('ready4df_pangle.txt', 'a+')
w.truncate(0)
fm.write_list(cp.alternate_input_4df('fiveplus_pdbs_s5.txt'), 'ready4df_pangle.txt')

#Creates file that will be read into datasets using
#df_pangles.py, df_hydro_plots.py, df_flex_plots


#=========================================================
#1.6. Flexibility grade: extracting PDBs by list in textfile
#=========================================================
#Reads list of pdb files within text file and extracts listed pdbs from folder

#provides bash script for extracting pdbs linked to Grade 1 flexibility
w = open('flex_code/bash_flex_g1_pdbs.txt', 'a+')
w.truncate(0)
fm.write_list(cp.write2_bash('flex_code/flex_grade1_list.txt', 'Grade_1'), 'flex_code/bash_flex_g1_pdbs.txt')

#provides bash script for extracting pdbs linked to Grade 2 flexibility
w = open('bash_flex_g2_pdbs.txt', 'a+')
w.truncate(0)
fm.write_list(cp.write2_bash('flex_code/flex_grade2_list.txt', 'Grade_2'), 'flex_code/bash_flex_g2_pdbs.txt')

#provides bash script for extracting pdbs linked to Grade 3 flexibility
w = open('bash_flex_g3_pdbs.txt', 'a+')
w.truncate(0)
fm.write_list(cp.write2_bash('flex_code/flex_grade3_list.txt', 'Grade_3'), 'flex_code/bash_flex_g3_pdbs.txt')

#provides bash script for extracting pdbs linked to Grade 3 flexibility
w = open('bash_flex_gR_pdbs.txt', 'a+')
w.truncate(0)
fm.write_list(cp.write2_bash('flex_code/flex_gradeR_list.txt', 'Grade_R'), 'flex_code/bash_flex_gR_pdbs.txt')

#Text has been copied into df_flex_hydroplots





#=========
#Extras
#=========
"""
cp.add_commas('complex_commas.txt','complexAb_commas.txt')

cp.chain_extractor('raw_complexed.txt')
w=open('Redundant_over2.txt','w')
w.truncate(0)
w.write(workshop.redund_over2(cg.redundant_pdbs))

cp.bash4_pack(cg.complex_list, cg.complex_pangle)

cp.chain_extractor('raw_complexed.txt')

cp.bash4_pack(cg.allpdbs)

clean_pdb.write2_csv(cg.complx_file)

clean_pdb.write2_bash()

#Creating groups of pdbs and extracting pangles from dataframe
w=open('extracted_group_pangles.txt', 'a+')

free_g=''
complex_g=''
group_no=0


#for index, words in enumerate(dfc.free_all_groups):

    # group_no+=1
    # print('fg_',str(group_no),'<-c(',words,')')
    
for index, words in enumerate(dfc.complex_all_groups):
    print('cg_',str(index+1),'<-c(',words,')\n')
"""
