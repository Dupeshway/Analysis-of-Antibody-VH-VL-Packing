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

hydropath=cg.path_2hydrocode


sns.set(style="darkgrid", palette="pastel", color_codes=True)

#==================
#1) Melting Dataframes
#==================

df1= pd.read_csv(hydropath+'/grade_1_hydrophob.csv')
df2= pd.read_csv(hydropath+'/grade_2_hydrophob.csv')
df3= pd.read_csv(hydropath+'/grade_3_hydrophob.csv')
dfR= pd.read_csv(hydropath+'/grade_R_hydrophob.csv')

df_mega = pd.read_csv(hydropath+'/mega_df.csv')

df_hydro = pd.read_csv(hydropath+'/hydrodfdf.csv')

df_melt = pd.melt(df_mega, id_vars=['Grade'], value_vars=['L','H'])
#df1= pd.melt(df1)


print(df_mega.head())

#Hydroscale:
h_value = df_aa_scale.h_value.rename('H-value')

#Sets
g1_set = df1.ab_set_id.rename('Grade_1 Ab Set ID')
g2_set = df2.ab_set_id.rename('Grade_2 Ab Set ID')
g3_set = df3.ab_set_id.rename('Grade_3 Ab Set ID')
gR_set = dfR.ab_set_id.rename('Grade_R Ab Set ID')


#Names
g1_pdb = df1.Grade_1.rename('Grade_1 PDB')
g2_pdb = df2.Grade_2.rename('Grade_2 PDB')
g3_pdb = df3.Grade_3.rename('Grade_3 PDB')
gR_pdb = dfR.Grade_R.rename('Grade_R PDB')

#Hydrophobicity
g1_light = df1.L.rename('Grade_1 Light_chain')
g1_heavy = df1.H.rename('Grade_1 Heavy Chain')
g1_total = df1.Total.rename('Grade_1 Total')


g2_light = df2.L.rename('Grade_2 Light Chain')
g2_heavy = df2.H.rename('Grade_2 Heavy Chain')
g2_total = df2.Total.rename('Grade_2 Total')

g3_light = df3.L.rename('Grade_3 Light Chain')
g3_heavy = df3.H.rename('Grade_3 Heavy Chain')
g3_total = df3.Total.rename('Grade_3 Total')


gR_light = dfR.L.rename('Grade_R Light Chain')
gR_heavy = dfR.H.rename('Grade_R Heavy Chain')
gR_total = dfR.Total.rename('Grade_R Total')


#Weighted Hydrophobicity

g1_weight_light = df1.L_weighted.rename('Grade_1 Weight Light Chain')
g1_weight_heavy = df1.H_weighted.rename('Grade_1 Weight Heavy Chain')

g2_weight_light = df2.L_weighted.rename('Grade_2 Weight Light Chain')
g2_weight_heavy = df2.H_weighted.rename('Grade_2 Weight Heavy Chain')

g3_weight_light = df3.L_weighted.rename('Grade_3 Weight Light Chain')
g3_weight_heavy = df3.H_weighted.rename('Grade_3 Weight Heavy Chain')

gR_weight_light = dfR.L_weighted.rename('Grade_R Weight Light Chain')
gR_weight_heavy = dfR.H_weighted.rename('Grade_R Weight Heavy Chain')

#===============================
#2) Total vs Flexibility (Grades)
#================================

sns.swarmplot(x='Grade', y='H-Value', hue='Chain', data=df_hydro)
#sns.violinplot(x='Grade', y='H-Value', hue='Chain', palette='muted', data=df_hydro)
#sns.violinplot(x='Grade',y='Total', palette='muted', data=df_mega)
#sns.violinplot(x='Grade',y='Total_weighted', data=df_mega)
sns.violinplot(x='Total',y='Grade', data=df_mega, legend=False)
sns.violinplot(x='Total',y='Grade', data=df_mega, legend=False)

plt.title('Weighted Total (Bright) and Total (Shaded)'
          + 'H-Values by Grade Flexibility')
plt.legend(loc='upper right')

plt.xlabel('Flexibility Grades')
plt.ylabel('Hydrophobicity value')
plt.show()



#==================================
#3) L and H correlation - no density
#===================================


#sns.jointplot(gR_light, gR_heavy, label='Grade_R Light_chain')
sns.jointplot(g2_light, g2_heavy, label='Grade_2 Light_chain')

sns.lmplot(g1_light, g1_set, data=df1['Total'])
#sns.lmplot(g2_light, g2_heavy)

sns.lmplot(x='L_weighted', y='H_weighted', hue='Grade', data=df_mega, legend=False)


plt.title('Correlation of Weighted Light and Heavy chain H-Values')
plt.legend(loc='upper right')

plt.xlabel('Weighted Light chain H-Value')
plt.ylabel('Weighted Heavy chain H-Value')
plt.show()


#===============================
#4) Weighted L and H correlation
#===============================


#sns.jointplot(gR_light, gR_heavy, label='Grade_R Light_chain')
sns.jointplot(g2_light, g2_heavy, label='Grade_2 Light_chain')

sns.lmplot(g1_light, g1_set, data=df1['Total'])
#sns.lmplot(g2_light, g2_heavy)

sns.lmplot(x='L_weighted', y='H_weighted', hue='Grade', data=df_mega, legend=False)


plt.title('Correlation of Weighted Light and Heavy chain H-Values')
plt.legend(loc='upper right')

plt.xlabel('Weighted Light chain H-Value')
plt.ylabel('Weighted Heavy chain H-Value')
plt.show()






