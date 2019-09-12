#!/usr/bin/env python
 
import csv
import sys
import pprint
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sy
import matplotlib.pyplot as plt

import re
import os
import shutil
from clean_pdb import clean_pdb as cp
from residue_extract import interface_residues as ir
from file_management import file_management as fm
from config import config as cg
from Bio.PDB import *


#=======================
#1) Extracting Residues
#=======================
# Generates list of functions for extracting Heavy and Light chains from each flexibility grade

#Prints lists found below in Grades 1/2/3/R extraction 
print(interface_residues.gen_heavy_chain_parser(cg.path_2flexcode+'flex_grade1_list.txt','Grade_1/','Grade_1L'))
print(interface_residues.gen_light_chain_parser(cg.path_2flexcode+'flex_grade1_list.txt','Grade_1/','Grade_1H'))

print(interface_residues.gen_light_chain_parser(cg.path_2flexcode+'flex_grade2_list.txt','Grade_2/','Grade_2L'))
print(interface_residues.gen_heavy_chain_parser(cg.path_2flexcode+'flex_grade2_list.txt','Grade_2/','Grade_2H'))

print(interface_residues.gen_light_chain_parser(cg.path_2flexcode+'flex_grade3_list.txt','Grade_3/','Grade_3L'))
print(interface_residues.gen_heavy_chain_parser(cg.path_2flexcode+'flex_grade3_list.txt','Grade_3/','Grade_3H'))

print(interface_residues.gen_heavy_chain_parser(cg.path_2flexcode+'flex_gradeR_list.txt','Grade_R/','Grade_RL'))
print(interface_residues.gen_light_chain_parser(cg.path_2flexcode+'flex_gradeR_list.txt','Grade_R/','Grade_RH'))


print(interface_residues.replace_aa_hydro_regex('/home/yobi/AbPackingAngle_V2.1/code/flex_code/grade2_residues.txt'))

#==================
#2) Grade 1 extraction
#==================
w = open(cg.path_2flexcode+'grade1_residues.txt', 'a+')
w.truncate(0)

Grade_1L=''
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'4M6M_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'4M6N_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'4M5Y_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'4M5Z_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'4NUG_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'5FUU_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'6DCQ_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'6MAR_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'5WTG_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'5WK2_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'5WK3_1.pdb')
Grade_1L += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_1/'+'5WTH_1.pdb')

Grade_1H=''
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'4M6M_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'4M6N_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'4M5Y_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'4M5Z_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'4NUG_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'5FUU_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'6DCQ_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'6MAR_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'5WTG_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'5WK2_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'5WK3_1.pdb')
Grade_1H += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_1/'+'5WTH_1.pdb')

print(Grade_1L)
print(Grade_1H)

fm.write_list(Grade_1L ,cg.path_2flexcode+'grade1_residues.txt')
fm.write_list(Grade_1H ,cg.path_2flexcode+'grade1_residues.txt')

#==================
#3) Grade 2 extraction
#==================
w = open(cg.path_2flexcode+'grade2_residues.txt', 'a+')
w.truncate(0)

Grade_2L=''
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4JDV_2.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4JAM_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4JAN_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4JDV_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4RFE_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4R97_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4R9Y_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4RFN_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'5FCU_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'6EA5_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'6EA7_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'6DZL_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'6DZM_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'6EA5_2.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'6EA7_2.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'5H30_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'5H32_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4UT9_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'5H37_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4M7Z_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4M93_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4MA1_1.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4MA1_3.pdb')
Grade_2L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_2/'+'4M7J_1.pdb')

Grade_2H=''
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4JDV_2.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4JAM_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4JAN_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4JDV_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4RFE_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4R97_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4R9Y_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4RFN_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'5FCU_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'6EA5_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'6EA7_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'6DZL_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'6DZM_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'6EA5_2.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'6EA7_2.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'5H30_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'5H32_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4UT9_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'5H37_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4M7Z_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4M93_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4MA1_1.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4MA1_3.pdb')
Grade_2H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_2/'+'4M7J_1.pdb')

print(Grade_2L)
print(Grade_2H)

fm.write_list(Grade_2L ,cg.path_2flexcode+'grade2_residues.txt')
fm.write_list(Grade_2H ,cg.path_2flexcode+'grade2_residues.txt')

#==================
#4) Grade 3 extraction
#==================
w = open(cg.path_2flexcode+'grade3_residues.txt', 'a+')
w.truncate(0)

Grade_3L=''
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'6AL4_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'6A76_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'6A77_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'6A78_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'6AL5_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1MQK_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1AR1_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1QLE_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'3EHB_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'3HB3_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'4R26_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'4R2G_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'5T3S_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'6IEQ_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1CFQ_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1BOG_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1CFN_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1CFS_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1CFT_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1HH6_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1HH9_1.pdb')
Grade_3L += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_3/'+'1HI6_1.pdb')

Grade_3H=''
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'6AL4_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'6A76_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'6A77_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'6A78_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'6AL5_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1MQK_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1AR1_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1QLE_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'3EHB_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'3HB3_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'4R26_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'4R2G_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'5T3S_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'6IEQ_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1CFQ_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1BOG_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1CFN_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1CFS_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1CFT_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1HH6_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1HH9_1.pdb')
Grade_3H += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_3/'+'1HI6_1.pdb')


print(Grade_3L)
print(Grade_3H)

fm.write_list(Grade_3L ,cg.path_2flexcode+'grade3_residues.txt')
fm.write_list(Grade_3H ,cg.path_2flexcode+'grade3_residues.txt')

#================
#5) Grade_R - Rigid
#================

w = open(cg.path_2flexcode+'gradeR_residues.txt', 'a+')
w.truncate(0)

Grade_RL=''
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'1RFD_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'1Q72_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'1QYG_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'1RIU_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'1RIV_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'1RZ7_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'3JWD_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'3JWO_1.pdb')
Grade_RL += ir.find_aa_heavychain(cg.path_2flexpdbs+'Grade_R/'+'4DVR_1.pdb')

Grade_RH=''
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'1RFD_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'1Q72_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'1QYG_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'1RIU_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'1RIV_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'1RZ7_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'3JWD_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'3JWO_1.pdb')
Grade_RH += ir.find_aa_lightchain(cg.path_2flexpdbs+'Grade_R/'+'4DVR_1.pdb')

print(Grade_RL)
print(Grade_RH)

fm.write_list(Grade_RL ,cg.path_2flexcode+'gradeR_residues.txt')
fm.write_list(Grade_RH ,cg.path_2flexcode+'gradeR_residues.txt')



sns.set(style="darkgrid", palette="pastel", color_codes=True)

df1= pd.read_csv(hydropath+'/grade_1_hydrophob.csv')
df2= pd.read_csv(hydropath+'/grade_2_hydrophob.csv')
df3= pd.read_csv(hydropath+'/grade_3_hydrophob.csv')
dfR= pd.read_csv(hydropath+'/grade_R_hydrophob.csv')

df_mega = pd.read_csv(hydropath+'/mega_df.csv')

df_hydro = pd.read_csv(hydropath+'/hydrodfdf.csv')

df_melt = pd.melt(df_mega, id_vars=['Grade'], value_vars=['L','H'])
#df1= pd.melt(df1)


print(df_mega.head())
