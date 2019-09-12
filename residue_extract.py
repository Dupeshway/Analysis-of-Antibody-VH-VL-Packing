#!/bin/env python3

'''Useful: https://realpython.com/working-with-files-in-python/'''

import re
import os
import shutil
import pandas as pd
from file_management import file_management as fm
from config import config as cg
import csv

path='/home/yobi/AbPackingAngle_V2.1/code'
path2flex='/home/yobi/AbPackingAngle_V2.1/code/flex_pdbs/'

class interface_residues:


    def replace_aa_hydro_regex(text_file):


        dic={"ILE": 00.730,
            "PHE":  00.610,
            "VAL":  00.540,
            "LEU":  00.530,
            "TRP":  00.370,
            "MET":  00.260,
            "ALA":  00.250,
            "GLY":  00.160,
            "CYS":  00.040,
            "TYR":  00.020,
            "PRO":  -0.070,
            "THR":  -0.180,
            "SER":  -0.260,
            "HIS":  -0.400,
            "GLU":  -0.620,
            "ASN":  -0.640,
            "GLN":  -0.690,
            "ASP":  -0.720,
            "LYS":  -1.100,
            "ARG":  -1.800}



        clean_data=''
        f = open(text_file)
        raw_data = f.read().splitlines()

        for index, item in enumerate(raw_data):
            print(item[:3])

            if item in dic.keys():
                clean_data += raw_data.replace(item, dic[item])

        return clean_data



    def replace_aa_hydro_v1(text_file):

        clean_data=''
        jeff=''

        f = open(text_file)
        raw_data = f.read().splitlines()

        for index, words in enumerate(raw_data):

            jeff = words.replace('TRP', '0.370', 1)+'\n'
            jeff = words.replace('MET', '0.260', 1)+'\n'
            jeff = words.replace('ALA', '0.250', 1)+'\n'
            jeff = words.replace('GLY', '0.160', 1)+'\n'
            jeff = words.replace('CYS', '0.040', 1)+'\n'
            jeff = words.replace('TYR', '0.020', 1)+'\n'

            clean_data+=jeff

        return clean_data



    def find_aa_lightchain(pdb_file):

        residues=pdb_file[-10:]+': L \n'

        f= open(pdb_file, 'r')
        raw_data = f.read().splitlines()
        count=0
        
        for index, words in enumerate(raw_data):

            if re.findall('L  43',words) and count==0:
                residues+=words[17:20]+' L43\n'
                count=1

            if re.findall('L  44',words) and count==1:
                residues+=words[17:20]+' L44\n'
                count=2


            if re.findall('L  45',words) and count==2:
                residues+=words[17:20]+' L45\n'
                count=3


            if re.findall('L  46',words) and count==3:
                residues+=words[17:20]+' L46\n'
                count=4



            if re.findall('L  96',words) and count==4:
                residues+=words[17:20]+' L96\n'
                count=5


            if re.findall('L  98',words) and count==5:
                residues+=words[17:20]+' L98\n'
                count=6

        return residues 


    def find_aa_heavychain(pdb_file):


        residues=pdb_file[-10:]+': H \n'

        f= open(pdb_file, 'r')
        raw_data = f.read().splitlines()
        count=0

        for index, words in enumerate(raw_data):


            if re.findall('H  44',words) and count==0:
                residues+=words[17:20]+' H44\n'
                count=1


            if re.findall('H  45',words) and count==1:
                residues+=words[17:20]+' H45\n'
                count=2


            if re.findall('H  46',words) and count==2:
                residues+=words[17:20]+' H46\n'
                count=3


            if re.findall('H  47',words) and count==3:
                residues+=words[17:20]+' H47\n'
                count=4


            if re.findall('H 101',words) and count==4:
                residues+=words[17:20]+' H101\n'
                count=5


            if re.findall('H 105',words) and count==5:
                residues+=words[17:20]+' H105\n'
                count=6


        return residues




    def gen_light_chain_parser(file_list, flex_folder, flex_grade):
        """
        Parses file_list 
        """

        f= open(file_list, 'r')
        raw_data = f.read().splitlines()

        print(flex_grade+"=''")
        
        residue_parser=''
        
        for i in raw_data:
            residue_parser+=flex_grade+' += ir.find_aa_lightchain(cg.path_2flexpdbs+'+"'"+flex_folder+"'+'"+i+"')\n"

        return residue_parser



    def gen_heavy_chain_parser(file_list, flex_folder, flex_grade):

        f= open(file_list, 'r')
        raw_data = f.read().splitlines()

        print(flex_grade+"=''")
        
        residue_parser=''
        
        for i in raw_data:
            residue_parser+=flex_grade+' += ir.find_aa_heavychain(cg.path_2flexpdbs+'+"'"+flex_folder+"'+'"+i+"')\n"

        return residue_parser





            
