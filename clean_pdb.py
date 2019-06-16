#!/bin/env python3

'''Useful: https://realpython.com/working-with-files-in-python/'''

import re
import os
import shutil
from config import config as cg
import csv

class clean_pdb:

    def write2_csv(data):
        '''WIP
        to extract and place in cells in a csv, complexed and non-complexed pdbs
        input: txt file with list of complex vs non-complex
        output: csv file with this information in seperate cells
        '''

        f = open(data, 'r')
        raw_data = f.read().splitlines()


        w = open(cg.csv_file,'a+')
        w.truncate(0)

        with open(cg.csv_file, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(raw_data)
            
        print('Data written to:', cg.csv_file)

clean_pdb.write2_csv(cg.complx_file)

"""

        for i in data4csv:
            with open(cg.csv_file, 'w') as csvFile:
                writer = csv.writer(csvFile)
                writer.writerows(i)

        csvFile.close()

    def write2_bash():
        '''COMPLETE
        Convert the file list of pdbs into bash script that would copy the files
        into another folder
        input: file with list of pdb files
        output: bash script for copying pdbs
        '''
        f=open(cg.w_file, 'r')
        data=f.read().splitlines()
        bashful=''
        
        w = open(cg.b_script,'a+')
        w.truncate(0)
        for i in data:
            w.write('cp -v '+i+' /home/yobi/Documents/Bioinformatics/Project/LH_combo/Redundant_pdbs\n')


        print('Data written to:', cg.b_script)


    def run_bash():
        '''Script to open and run a bash file, for cohesion when running the entire
        program, MAY NOT BE NECESSARY
        '''

    def capture_pdb():
        '''COMPLETE
        capture files from folder than only have redundancies
        poss.input1: file with list of wanted pdbs,
        poss.input2: folder with mixed pdb files
        output: folder with redundant pdb files
        '''
        file= cg.redund_file
        f=open(file,'r')
        raw_data=f.read().splitlines()
        
        clean_data=[]
        rough_data=[]
        seperator=', '
        pdb='.pdb'

        for line in raw_data:
            rough_data+=line.split(seperator)

        w = open(cg.w_file,'a+')
        w.truncate(0)
        for i in rough_data:
            w.write(i[0:6]+'.pdb\n')
        print('Data written to:', cg.w_file)

    

    def redund_3plus(file):
        '''COMPLETED
        open a txt file with lists of redundant and non-redundant pdbs
        input: text file with mixed pdbs with redundant pdbs on the same line
        output: text file with only redundant pdbs
        '''
        f=open(file,'r')
        raw_data=f.read().splitlines()
        redund_pdbs=''
        p=re.compile(r'1,\s')
        
        for line in raw_data:
            if p.findall(line):
                redund_pdbs+=line+'\n'
        
        w_file='sort_redund_pdbs.txt'

        w = open(w_file,'a+')
        w.truncate(0)
        w.write(redund_pdbs)
        print('written to: ', w_file)


    def list_directory():
        '''WiP
        input: should be folder directory, unlikely to use anything other than the same folder
        capture files from folder than only have redundancies
        output: folder with redundant pdb files
        '''

        folder= '/home/yobi/Documents/Bioinformatics/Project/LH combined CHothia/LH_Combined_Chothia/'
        with os.scandir(folder) as entries:
            gary=''
        with os.scandir(folder) as entries:
            for entry in entries:
                gary+=entry.name+'\n'
            print(gary)



"""
