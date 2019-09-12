#!/usr/bin/env python
 
import csv
import sys
import pprint
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


#=====================
#1) Melting Dataframes
#=====================

sns.set(style="whitegrid", palette="muted")
#Initial DF, used for Mean densities and correltaion
df= pd.read_csv('dataframe_all_groups.csv')
#m_df = pd.melt(df, id_vars=['all_ab_groups'])
#df_index=df.set_index('all_ab_groups',drop=True) #purely by the index


#Mean DFs
dfMean      = pd.read_csv('melted_mean_datasets.csv')
dfM_FvsC    = pd.read_csv('melted_mean_freecomp_dataset.csv')
dfM_red_467 = pd.read_csv('melted_mean_red46_7.csv')


#Range Dfs
dfRange     = pd.read_csv('melted_range_datasets.csv')
dfR_FvsC    = pd.read_csv('melted_range_freecomp_dataset.csv')
dfR_red_467 = pd.read_csv('melted_range_red46_7.csv')

melted_dfR_all = pd.melt(dfRange,
                         id_vars="Range",
                         var_name="All_pdb")

melted_dfR_free = pd.melt(dfRange,
                         id_vars="Range",
                         var_name="Free_pdb")

melted_dfR_complexed = pd.melt(dfRange,
                         id_vars="Range",
                         var_name="Complexed_pdb")


print(melted_dfR_all.head())

#===================================
#2) MEAN free vs comp across all sets
#===================================

#sns.swarmplot(x='Ab_set_id', y='Mean', hue='Data_set', data=dfM_FvsC)
sns.distplot(x='Data_set', y='Mean', data=dfM_FvsC)
#sns.swarmplot(x='Mean', y='Redundancy', hue='Data_set', data=dfM_FvsC)


plt.title('Mean Variance across Datasets with redundancy count')
#plt.legend(loc='upper left')

plt.xlabel('Ab Data Sets')
plt.ylabel('Mean P-Angle')

plt.show()


#============================
#3) Mean of largest redundancies
#============================

sns.violinplot(x='Data_set', y='Mean', data=dfM_red_467)

plt.title('Mean Variance across All_pdb Redundancy counts')
#plt.legend(loc='upper right')

plt.xlabel('Ab Sets - 66 total')
plt.ylabel('Mean P-Angle')

#=============================
#4) MEAN variance across all sets
#=============================
sns.swarmplot(x='Ab_set_id', y='Mean', hue='Data_set', data=dfMean)


plt.title('Mean Variance across all sets')
plt.legend(loc='upper right')

plt.xlabel('Ab Sets - 66 total')
plt.ylabel('Mean P-Angle')




#==========================
#5) Range: Free / Comp / All
#==========================

sns.distplot(melted_dfR_free.Range, label='All')
sns.distplot(melted_dfR_free.Range, label='Free')
sns.distplot(melted_dfR_free.Range, label='Complexed')
#sns.swarmplot(x='Data_set', y='Range', hue='Redundancy Category', data=dfR_FvsC)

plt.title('Range Variance across Datasets and Redundancy counts')

plt.legend(loc='upper left')


plt.xlabel('All Ab sets by Redundancy count')
plt.ylabel('P-Angle Range')





#==========================
#6) Range Free / Comp / All
#=========================

sns.violinplot(x='Data_set', y='Range', hue='Redundancy Category', data=dfR_FvsC, legend_out=True)
#sns.swarmplot(x='Data_set', y='Range', hue='Redundancy Category', data=dfR_FvsC)

plt.title('Range Variance across Datasets and Redundancy counts')

plt.legend(loc='upper left')


plt.xlabel('All Ab sets by Redundancy count')
plt.ylabel('P-Angle Range')



#======================================
#7) Range redundancy plot (all 4-6 and 7+)
#======================================

sns.swarmplot(x='Data_set', y='Range', data=dfR_red_467)

plt.title('Range Variance across All Redundancy counts')
plt.legend(loc='upper right')

plt.xlabel('All Ab sets by Redundancy count')
plt.ylabel('P-Angle Range')








#===============================
#8) Range density plot (all 3 in 1)
#==============================

sns.kdeplot(range_all, bins=30, shade=True)
sns.kdeplot(range_free_nan, shade=True)
sns.kdeplot(range_complex_nan, shade=True)



plt.title('Density of Range for All Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('Density')
plt.show()


#=======================================
#9) Range against redund count (all 3 in 1)
#=======================================

sns.jointplot(all_redund_count, range_all, kind="reg", space=0, color='blue')
sns.jointplot(free_redund_count, range_free, kind="reg", space=0, color='g')
sns.jointplot(complex_redund_count, range_complex, kind="reg", space=0, color='r')

plt.show()

#============================
#10) Standard deviations
#===========================

sns.barplot(x='Data_set', y='SD', hue='Redundancy Category', data=dfR_FvsC)
#sns.swarmplot(x='Data_set', y='Range', hue='Redundancy Category', data=dfR_FvsC)

plt.title('SD across Datasets and Redundancy counts')

plt.legend(loc='upper left')


plt.xlabel('All Ab sets by Redundancy count')
plt.ylabel('SD P-Angle')


plt.show()

#==================================
#11) SD of Density plot(all 3 in 1)
#==================================


sns.kdeplot(sd_all)
sns.distplot(sd_all, label='All 4+ Redundancies')
sns.distplot(sd_all_10plus, label='All 10+ Redundancies', color='violet')


sns.jointplot(range_free_nan, sd_free_nan, kind='reg', dropna=True, color='g')
sns.jointplot(range_complex_nan, sd_complex_nan, kind='reg', dropna=True, color='r')

sns.distplot(range_free_nan, kind='reg', color='purple', dropna=True)


plt.title('Density of SD for Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('SD in degrees')
plt.ylabel('Density')
plt.show()


#===========================
#12) SD and Range Scatter plot
#===========================

#Plots
#sns.jointplot(range_all, sd_all, kind="reg", space=0, color='seagreen')
sns.kdeplot(sd_all, shade=True)
sns.kdeplot(range_all, shade=True)

sns.jointplot(range_all, sd_all, kind='reg', space=0)

plt.title('Correlation: Range and SD for all Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('SD in degrees')
plt.show()


#===============================
#13) SD Range Free vs Complex
#===============================
sns.kdeplot(sd_free, shade=True)
sns.kdeplot(sd_complex, shade=True)

plt.title('Correlation: Range and SD for all Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('SD in degrees')
plt.show()





#time_temp2=sns.swarmplot(x='Time step (secs)',y='mosfet_temperature (degc)', data=df, color='green')

#plt.show(all_pdb_pangles)

"""
#==========
#14) ARCHIVE:
#==========
#Redundancy count of All, Free, Complex
all_redund_count      =df.all_redund_count.rename('All Redundancies Count')
free_redund_count     =df.free_redund_count.rename('Free Redundancies Count')
complex_redund_count  =df.complex_redund_count.rename('Complex Redundancies Count')

#Mean of All, Free, Complex
mean_all          = df.mean_all_groups.rename('Mean of All Ab sets')
mean_free         = df.mean_free_groups.rename('Mean of Free sets')
mean_complex      = df.mean_complex_groups.rename('Mean of Complex sets')
median_all        = df.all_median.rename('Median of Each Abs')

#Range of All, Free, Complex
range_all         = df.range_all_groups.rename('Range All sets')

range_free        = df.range_free_groups.rename('Range Free sets')



range_complex     = df.range_complex_groups.rename('Range Complex Sets')

#SD of All, Free, Complex
sd_all            = df.sd_all_groups.rename('SD All Sets')


sd_free           = df.sd_free_groups.rename('SD Free Sets')


sd_complex        = df.sd_complex_groups.rename('SD Complex Sets')


#create range vs sd df
#data_sd_range   = [[range_all, sd_all]]
#df_sd_range     = pd.DataFrame(data_sd_range)



#Mean density plot (all 3 in 1)
#=============================

sns.distplot(mean_all, bins=15)
#sns.kdeplot(mean_free)
#sns.kdeplot(mean_complex)

plt.title('Density of Means for All Antibody groups')
plt.legend(loc='upper left')

plt.xlabel('Mean in degrees')
plt.ylabel('Density')
plt.show()

#Mean Scatter plot
#=================

sns.jointplot(mean_free, mean_complex, kind="reg", space=0, color='purple')

#plt.title('Scatter of Means for Antibody groups')
#plt.legend(loc='upper left')


#plt.xlabel('Mean in degrees')
#plt.ylabel('Density')
plt.show()



#Mean across redundancy numbers:
#==============================

sns.jointplot(all_redund_count, mean_all, kind="reg", space=0, color='blue')
sns.jointplot(free_redund_count, mean_free, kind="reg", space=0, color='g')
sns.jointplot(complex_redund_count, mean_complex, kind="reg", space=0, color='r')


plt.title('Scatter of Means for Antibody groups')
plt.legend(loc='lower left')


#plt.xlabel('Mean in degrees')
#plt.ylabel('Density')
plt.show()



#Median plot
#==========

sns.kdeplot(no_groups, median_all)
plt.show()




#Range density plot (all +4-9+ and 10+)
#==============================

print(melted_df['range_all_groups'])

sns.distplot(range_all, bins=7, label='All 4+ Redundancies')
sns.distplot(range_all_10plus, bins=7, label='All 10+ Redundancies', color='violet')

plt.title('Density of Range for All Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('Density')
plt.show()




#Range density plot (all 3 in 1)
#==============================

#sns.kdeplot(range_all, bins=30, shade=True)
#sns.kdeplot(range_free_nan, shade=True)
#sns.kdeplot(range_complex_nan, shade=True)



plt.title('Density of Range for All Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('Density')
plt.show()


print(df.head())
#Range against redund count (all 3 in 1)
#==============================

sns.jointplot(all_redund_count, range_all, kind="reg", space=0, color='blue')
sns.jointplot(free_redund_count, range_free, kind="reg", space=0, color='g')
sns.jointplot(complex_redund_count, range_complex, kind="reg", space=0, color='r')

plt.show()


#SD of Density plot(all 3 in 1)
#===============================


#sns.kdeplot(sd_all)
#sns.distplot(sd_all, label='All 4+ Redundancies')
#sns.distplot(sd_all_10plus, label='All 10+ Redundancies', color='violet')


sns.jointplot(range_free_nan, sd_free_nan, kind='reg', dropna=True, color='g')
#sns.jointplot(range_complex_nan, sd_complex_nan, kind='reg', dropna=True, color='r')

#sns.distplot(range_free_nan, kind='reg', color='purple', dropna=True)


#plt.title('Density of SD for Antibody groups')
#plt.legend(loc='upper right')

#plt.xlabel('SD in degrees')
#plt.ylabel('Density')
plt.show()



#SD and Range Scatter plot
#===========================

#Plots
#sns.jointplot(range_all, sd_all, kind="reg", space=0, color='seagreen')
sns.kdeplot(sd_all, shade=True)
sns.kdeplot(range_all, shade=True)

sns.jointplot(range_all, sd_all, kind='reg', space=0)

plt.title('Correlation: Range and SD for all Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('SD in degrees')
plt.show()



#SD Range Free vs Complex
#===============================
sns.kdeplot(sd_free, shade=True)
sns.kdeplot(sd_complex, shade=True)

plt.title('Correlation: Range and SD for all Antibody groups')
plt.legend(loc='upper right')

plt.xlabel('Range in degrees')
plt.ylabel('SD in degrees')
plt.show()





#time_temp2=sns.swarmplot(x='Time step (secs)',y='mosfet_temperature (degc)', data=df, color='green')

#plt.show(all_pdb_pangles)
"""

