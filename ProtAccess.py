#!/usr/bin/env python

import sys, re
import os
from argparse import ArgumentParser
import Bio
from Bio import Entrez
from Bio import SeqIO
import Bio.SeqUtils
from io import StringIO
from collections import OrderedDict
import requests
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

parser = ArgumentParser(description = 'Extract information from UniProt and NCBI websites based on UniProt ID or list of UniProt IDs obtained from website.')
parser.add_argument("-u", "--uniprotID", type = str, required = False, help = "Single UniProt ID. Used to obtain information on one protein only.")
parser.add_argument("-i", "--inputfile", type = str, required = False, help = "File with UniProt IDs downloaded from uniprot.org as list file. Used to obtain information on a group of proteins.")
parser.add_argument("-n", "--ncbiSeq", action='store_true', required = False, help = "Generate fasta file(s) with RNA trascript sequences.")
parser.add_argument("-p", "--proteinSeq", action='store_true', required = False, help = "Generate fasta file(s) with protein sequences. Must be used with -s and -l.")
parser.add_argument("-s", "--isoform", type = str, required = False, help = "Used with -p and -l. Set as 'canonical' to obtain canonical protein sequence only and 'all' to obtain protein sequences of all isoforms.")
parser.add_argument("-l", "--aa_length", type = str, required = False, help = "Used with -p and -s. Set as 'y' to obtain lenght of amino acid sequence of each isoform and 'n' if such information is not needed.")
parser.add_argument("-g", "--GO", action='store_true', required = False, help = "Stores Gene Ontology terms from the UniProt page of each protein. Needs flag -t to specify GOterm. Needs flag -c OR -o to generate output")
parser.add_argument("-t", "--GOterm", type = str, required = False, help = "Used with -g. Set as F for GO Molecular Function, P for GO Biological Process, or C for Cellular Component")
parser.add_argument("-c", "--countplot", action='store_true', required = False, help = "Used with -g and -t to create a countplot of GO terms associated with proteins in list. Only used for list of UniProt IDs.")
parser.add_argument("-x", "--GOnum", type = int, required = False, help = 'Used with -g and -t. Accepts integers. When given a number x, it gives the top x GO terms for each protein from its UniProt page. If less than x GO terms are found for a protein, it gives all GO terms')


if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

#This is needed to access Entrez via Biopython
Entrez.email = 'danasalim.fakhreddine@uvic.cat'
handle_main = Entrez.einfo()
record_main = Entrez.read(handle_main)

#Functions 
def uniprotSum(code):
    ''' Function that gives a list of information on protein given a uniprot ID. 
    
    First element in list is the description of the protein, second element is the gene name, 
    
    and the third until the last are the NCBI transcript codes. 
    
    code  -> str
    '''
    url = "http://www.uniprot.org/uniprot/" + code + ".txt"
    req = requests.get(url) 
    req
    lines = req.text.splitlines() 
    gene_sum = []
    for line in lines:
        if re.search('DE   RecName: Full=', line):
            if re.search('{', line): 
                gene_sum.append(line[19:line.index('{')])
            else:
                gene_sum.append(line[19:line.index(';')])
        if re.search('GN   Name=', line):
            if re.search('{', line):
                gene_sum.append(line[10:line.index('{')])
            else:
                gene_sum.append(line[10:line.index(';')])
        if re.search('NM_', line):
            if re.search(code, line):
                gene_sum.append(line[line.index('NM_'): line.index(code) -3])
            else:
                gene_sum.append(line[line.index('NM_'):-1])
    return(gene_sum)

def ncbiSeq(code):
    '''Function that gives a fasta file with the trascript RNA sequence(s) given the uniprot ID. 
    
    It also provides the protein name, Gene ID, and total number of trascripts 
    
    code -> str
    '''
    sum_code = uniprotSum(code)
    NM = sum_code[2:]
    NM_len = len(sum_code) - 2
    title = str(code) + '/'+str(code) + '_transcripts.fasta'
    f= open(title,"w+") 
    for i, trascripts in enumerate(NM):
        handle = Entrez.efetch(db="nuccore", id=trascripts, rettype="fasta", retmode="text")
        nucSeq = handle.read()
        final ='\n'+nucSeq
        f.write(final)
    print("\nUniprot ID: ", code, '\nProtein name: ',sum_code[0], "\nGene: ", sum_code[1], "\nTranscripts: ", NM_len, '\n---------------------------------\n', "File ", title, " created!")
    f.close()
    

def uniprotSeq(code, seq , AA_len):
    ''' 
    Function that poduces fasta file with protein sequences of the uniprot accession number (code) provided.
    
    If seq = "canonical" on the canonical sequence is analyzed. If seq = "all", all sequences are provided. 
    
    If AA_len is 'y', length info of each sequnece is provided. If AA_len is 'n', no information is provided.
    
    
    code -> str
    seq -> str
    AA_len -> str'''
    if seq == 'canonical':
        url = "http://www.uniprot.org/uniprot/" + code + ".fasta?include=no"
    if seq == 'all':
        url = "http://www.uniprot.org/uniprot/" + code + ".fasta?include=yes"

    req = requests.get(url) 
    lines = req.text.split('\n')
    title = str(code) + '/'+code + '_protein_seq' + '_' + seq + '.fasta'
    f = open(title, "w+")
    string = ''
    import re
    for line in lines:
        if '>' in line:
            string= string + '\n' + line + '\n'
        if '>' not in line:
            string = string + line + '\n'
    f.write(string)
    f.close()
    
    if AA_len is 'y':
        from Bio import SeqIO
        import Bio.SeqUtils
        from io import StringIO
        fasta_io = StringIO(string)
        fasta_sequences = SeqIO.parse(fasta_io,'fasta')
        fasta_sequences
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            AA = 'No. of AA: '+ str(len(sequence))
            print(name, '\n', AA)
        print('\n-----------------------\n')
        print('File ', title, ' created!')
        print('\n- - - - - - - - - - - - - - - - - - - - - - - -\n')
    if AA_len is 'n':
        print('File ',title, ' created!')
        print('\n- - - - - - - - - - - - - - - - - - - - - - - -\n')    

def removeDuplicates(listofElements):
    '''Function used to remove duplicates within a list of strings
    
    listofElements -> list 
    '''
    uniqueList = []
    for elem in listofElements:
        if elem not in uniqueList:
            uniqueList.append(elem)

    return uniqueList


def uniprotGO(code, GOterm):
    '''code is the uniprot ID and GOterm is either 'F' 'C' or 'P'
    
    F => GO Molecular Function 
    P => GO Biological Process
    C => Cellular Component
    
    code -> str
    GOterm -> str 
    
    '''
    url = "http://www.uniprot.org/uniprot/" + code + ".txt"
    req = requests.get(url)
    lines = req.text.splitlines()
    GO = []
    uniprot = []
    for index,line in enumerate(lines):
        if re.search(' GO;', line):
            GO.append(line)
    GO_term = GOterm+ ':'
    for i in GO:
        j = i[21:]
        if re.search(GO_term, j):
            if ';' in j:
                uniprot.append(j[2:j.index(';')])
            else: 
                uniprot.append(j[2:], '\n')
    return(uniprot)

#To give full name of GO terms
if args.GOterm == "F":
    term = 'GO Molecular Function'
if args.GOterm == "P":
    term = 'GO Biological Process'
if args.GOterm == "C":
    term = 'Cellular Component'


#To obtain information on a single protein
if args.uniprotID:
    try:
        os.mkdir(args.uniprotID)
    except:
        pass

    if args.ncbiSeq:

        print("\n~~~NCBI transcript information~~~\n")

        ncbiSeq(code = args.uniprotID)
    if args.proteinSeq:
        print("\n~~~Uniprot protein information~~~\n")
        uniprotSeq(code =args.uniprotID, seq = args.isoform ,AA_len = args.aa_length)
    if args.GO:
        print("\n~~~Gene Ontology information~~~\n")
        print('\nThe first', args.GOnum, term, "terms.\nIf less than", args.GOnum, "terms appear, these represent all", term, "terms for this protein.\n")
        GOs = uniprotGO(args.uniprotID, GOterm = args.GOterm)
        if len(GOs) < args.GOnum:
            print(GOs)
        else:
            print(GOs[:args.GOnum])


#To obtain information on a list of proteins with UniProt ID stored each on a line in a file
if args.inputfile:
    dictGO = OrderedDict()
    with open(args.inputfile) as list:
        for line in list:
            line =line.strip('\n')
            if args.ncbiSeq:
                try:
                    os.mkdir(line)
                except:
                    pass
                print("--------",line,"--------\n")
                print("~~NCBI~~")
                ncbiSeq(line)
                print('\n')
            if args.proteinSeq:
                try:
                    os.mkdir(line)
                except:
                    pass
                print("~~UniProt~~\n")
                uniprotSeq(line, seq = args.isoform ,AA_len = args.aa_length)
            if args.GO:
                dictGO[line.strip('\n')] = uniprotGO(line.strip('\n'), GOterm = args.GOterm)
            




    if args.countplot:
        df = pd.DataFrame(dictGO.values())
        df.index = dictGO.keys()
        all_GO = []
        for values in dictGO.values():
            for value in values:
                all_GO.append(value)
        unique_GO =removeDuplicates(all_GO)
        occurence = []

        for GO in unique_GO:
            x = all_GO.count(GO)
            occurence.append(int(x))
        df1 = pd.DataFrame(data= occurence, index= unique_GO)
        plt.subplots(figsize=(df1.iloc[0].max()*4, len(df1)/5))
        plot = sns.countplot(y=all_GO, data=df1, palette = 'pastel')
        fig_title = args.GOterm + "_GO_terms.png"
        plot_title = term + " countplot"
        plt.title(plot_title)

        plt.savefig(fig_title, bbox_inches='tight')
        print(fig_title, "created!")
    
    
    if args.GOnum:
        print('\nThe first', args.GOnum, term, "terms.\nIf less than", args.GOnum, "terms appear, these represent all", term, "terms for this protein.\n")
        for key in dictGO:
            print(key)
            if len(dictGO[key]) < args.GOnum:
                print(dictGO[key],'\n')
                
            else:
                print(dictGO[key][:args.GOnum], '\n')
    



            
