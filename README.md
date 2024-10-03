# **EIPpred**
A computational approach to predict the MIC values of inhibitory peptides against E.coli using the amino acid sequence information.

## Introduction
EIPpred is developed to predict and design the inhibitory peptides and scan the inhibitory regions in a protein sequence. In the standalone version, the Random Forest regressor-based model has been implemented. 
EIPpred is also available as a web server at https://webs.iiitd.edu.in/raghava/eippred. Please read/cite the content about the EIPpred for complete information, including the algorithm behind the approach.

## PIP installation
The PIP version is also available for easy installation and usage of this tool. The following command is required to install the package

```
pip install eippred
```

To know about the available option for the pip package, type the following command:

```
eippred -h
```
## Standalone
The Standalone version of EIPpred is written in python3, and the following libraries are necessary for the successful run:
- scikit-learn ( Version = 1.4.2 or less )
- Pandas (Version = 1.5.3 or less )
- Numpy (Version = 1.26.4 or less )

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
                  [-w {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}]
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
                        Job Type: 1:Predict, 2: Design, 3: Protein scan, by default 1
   -w {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30},  --winleng {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
                       Window Length: 8 to 30 (scan mode only), by default 9
```

**Input File:** It allows users to provide input in the FASTA format.

**Output File:** The Program will save the results in the CSV format; if the user does not provide the output file name, it will be stored in "outfile.csv".

**Job:** User is allowed to choose between two different modules, such as 1 for prediction, 2 for Designing and 3 for Protein scanning; by default, it's 1.

**

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
