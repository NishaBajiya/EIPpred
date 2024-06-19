# **EIPpred**
A computational approach to predict the MIC values of inhibitory peptides against E.coli using the amino acid sequence information.
## Introduction
EIPpred is developed to predict and design the inhibitory peptides. In the standalone version, the Random Forest regressor-based model has been implemented. 
EIPpred is also available as a web server at https://webs.iiitd.edu.in/raghava/eippred. Please read/cite the content about the EIPpred for complete information, including the algorithm behind the approach.

## Standalone
The Standalone version of EIPpred is written in python3, and the following libraries are necessary for the successful run:
- scikit-learn
- Pandas
- Numpy

## Minimum USAGE
To know about the available option for the standalone, type the following command:
```
python eippred.py -h
```
To run the example, type the following command:
```
python3 eippred.py -i example_input.fa
```
This will predict the MIC values of the submitted sequences, which will help identify the inhibitory activity of the peptides against E.coli. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma-separated variables).

## Full Usage
```
usage: eippred.py [-h] 
                  [-i INPUT]
                  [-o OUTPUT]
                  [-j {1,2,3}]
                  [-d {1,2}]
```
```
Please provide the following arguments for the successful run.

Optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -j {1,2,3}, --job {1,2,3}
                        Job Type: 1:Predict, 2: Design, by default 1
  -p POSITION, --Position POSITION
                        Position: Define the position of the mutation (1-indexed)
  -r RESIDUES, --Residues RESIDUES
                        Residues: define the residue to be mutated (one or two from the 20 essential amino acids in upper case)
```

**Input File:** It allows users to provide input in the FASTA format.

**Output File:** The Program will save the results in the CSV format; if the user does not provide the output file name, it will be stored in "outfile.csv".

**Job:** User is allowed to choose between two different modules, such as 1 for prediction and 2 for Designing; by default, it's 1.

**Position**: User can choose any position from the long sequence for mutation. This option is available for only the Design module.

**Residues:** This option allows users to incorporate the mutation of single amino-acid and dipeptide amino-acid residue from the original peptide sequences at the specific position defined by the user.

EIPpred Package Files
=======================
It contains the following files; a brief description of these files is given below

INSTALLATION                    : Installations instructions

LICENSE                         : License information

README.md                       : This file provides information about this package

eippred.py                      : Main Python program

example_input.fa                : Example file contains peptide sequences in FASTA format

example_predict_output.csv      : Example output file for predict module

example_design_output.csv       : Example output file for design module
