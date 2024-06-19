import argparse
import warnings
import subprocess
import pkg_resources
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
from collections import Counter
import pickle
import re
import glob
import time
import uuid
from time import sleep
from tqdm import tqdm
from sklearn.ensemble import RandomForestRegressor
import zipfile


# Function to check the seqeunce
def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
        seqid.append('>'+name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append(">Seq_"+str(i))
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    return df1,df2

# Function to check the length of seqeunces
def lenchk(file1):
    cc = []
    df1 = file1
    df1.columns = ['seq']
    for i in range(len(df1)):
        if len(df1['seq'][i])>30:
            cc.append(df1['seq'][i][0:30])
        else:
            cc.append(df1['seq'][i])
    df2 = pd.DataFrame(cc)
    df2.to_csv('out_len', index = None , header = None)
    df2.columns = ['Seq']
    return df2
# Function to generate the features out of seqeunces
# file = pd.read_csv('example.csv')


# Function to read and implement the model
def ML_run(file1, out):
    a=[]
    df = pd.read_csv(file1)
    ff = pd.read_csv('Data/selected_features_mrmr1000_new.csv')
    ff2 = ff["SelectedFeatures"].tolist()
    aa = pd.concat([df[ff2]], axis =1)
#     aa.to_csv('out_selected', index = None)
    clf = pickle.load(open('Data/eippred.sav','rb'))
    data_test = aa#pd.read_csv(file_name)
    X_test = data_test
    y_p_score1=clf.predict(X_test)
    y_p_s1=y_p_score1.tolist()
    a.extend(y_p_s1)
    df = pd.DataFrame(a)
    df1 = df.iloc[:,-1].round(3)
    df2 = pd.DataFrame(df1)
    df2.columns = ['MIC']
    dd = pd.concat([aa,df2], axis =1)
    dd.to_csv(out, index = None)
    return df2

def emb_process(file):
    df = pd.read_csv(file)
    df.insert(0, 'seq_ID', ['seq_' + str(i) for i in range(1, len(df) + 1)])
    ss = df[['seq_ID']]#.to_csv(sys.argv[2], header = None, index = None)
    df2 = df.drop(['seq_ID'], axis =1)
    colNumber = df2.shape[1]
    headerRow=[]
    for i in range(colNumber):
        headerRow.append('prot'+str(i))
    df2.columns=headerRow
    df3 = pd.concat([ss,df2], axis =1)
    return df3

def generate_mutant(original_seq, residues, position):
    std = "ACDEFGHIKLMNPQRSTVWY"
    if all(residue.upper() in std for residue in residues):
        if len(residues) == 1:
            mutated_seq = original_seq[:position-1] + residues.upper() + original_seq[position:]
        elif len(residues) == 2:
            mutated_seq = original_seq[:position-1] + residues[0].upper() + residues[1].upper() + original_seq[position+1:]
        else:
            print("Invalid residues. Please enter one or two of the 20 essential amino acids.")
            return None
    else:
        print("Invalid residues. Please enter one or two of the 20 essential amino acids.")
        return None
    return mutated_seq

def generate_mutants_from_dataframe(df, residues, position):
    mutants = []
    for index, row in df.iterrows():
        original_seq = row['Seq']
        mutant_seq = generate_mutant(original_seq, residues, position)
        if mutant_seq:
            mutants.append((original_seq, mutant_seq,position))
    return mutants

print('############################################################################################')
print('# This program EIPPred is developed for predicting, desigining and scanning MIC of peptides #')
print('# mellitus causing  peptides, developed by Prof G. P. S. Raghava group.               #')
print('# Please cite: EIPPred; available at https://webs.iiitd.edu.in/raghava/eippred/  #')
print('############################################################################################')

parser = argparse.ArgumentParser(description='Please provide following arguments')

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3], help="Job Type: 1:Predict, 2: Design, by default 1")
parser.add_argument("-p",'--Position', type=int, help='Position of mutation (1-indexed)')
parser.add_argument("-r",'--Residues', type=str, help='Mutated residues (one or two of the 20 essential amino acids in upper case)')


args = parser.parse_args()



# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 

# Output file 

if args.output == None:
    result_filename= "outfile.csv"
else:
    result_filename = args.output

# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)

position = args.Position
residues = args.Residues
#======================= Prediction Module start from here =====================
if Job == 1:
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    os.system('python3 '+'Data/composition_calculate.py -i'+'out_len -o'+ 'out2')
    mlres = ML_run('out2', 'out4')
    df3 = pd.concat([df_2,df1,mlres],axis=1)
    df3.to_csv(result_filename,index = None)

    os.remove('out_len')
    os.remove('out2')
    os.remove('out4')
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")

#===================== Design Model Start from Here ======================
elif Job == 2:
    print('\n======= Thanks for using Design module of EIPPred. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)

    while not all(residue.upper() in "ACDEFGHIKLMNPQRSTVWY" for residue in residues) or len(residues) > 2:
        print("Invalid input. Please enter one or two of the 20 essential amino acids in upper case.")
        residues = input("Enter the mutated residues: ")

    mutants = generate_mutants_from_dataframe(df1, residues, position)
    result_df = pd.DataFrame(mutants, columns=['Original Sequence','Mutant Sequence','Position'])
    result_df['Mutant Sequence'].to_csv('out_len_mut', index = None, header =None)
    os.system('python3 '+'Data/composition_calculate.py -i'+'out_len -o'+ 'out2')
    os.system('python3 '+'Data/composition_calculate.py -i'+'out_len_mut -o'+ 'out3')
    mlres = ML_run('out2', 'out22')
    mlres_m = ML_run('out3','out33')
    df3 = pd.concat([df_2,result_df['Original Sequence'],mlres,result_df[['Mutant Sequence', 'Position']], mlres_m ],axis=1)
    df3.columns = [['ID', 'Original Sequence', 'Ori_MIC','Mutant Sequence','Position','Mut_MIC']]
    df3.to_csv(result_filename, index = None)
    os.remove('out_len')
    os.remove('out_len_mut')
    os.remove('out2')
    os.remove('out22')
    os.remove('out3')
    os.remove('out33')
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")

print('\n======= Thanks for using EIPPred. Your results are stored in file :',result_filename,' =====\n\n')
