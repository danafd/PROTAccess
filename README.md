# ProtAccess

Online biological databases have provided scientists with a plethora of information that can be extracted, manipulated and studied. However, the huge amount of data can sometimes be daunting to deal with. ProtAccess is a Python tool that helps researchers access two widely used online databases, UniProt and NCBI, to extract key information such as protein and nucleotide sequences as well as Gene Ontology terms by inputting UniProt IDs only (Single ID or list of IDs stored in file). 

## Overview 

The `--uniprotID` flag can be used to input a single UniProt ID, and the `--inputfile` can be used to input a file containing a list of UniProt IDs downloaded directly as a **list** from UniProt. 

This tool creates a folder for each UniProt ID (protein) in your directory. Upon using the `--ncbiSeq` flag, a fasta file containing the trascript sequences extracted from NCBI is stored. Upon using the `--proteinSeq` flag followed by the `--isoform` (canonical or all) and `--aa_length` (y or n) flags, a fasta file containing the amino acid sequences of the specified isoforms is created and the lenght of each isoform is displayed (for '`--aa_length` y'). This information is extracted from UniProt. 


The Gene Ontology terms associated with each protein can also be extracted from UniProt using the `--GO` flag. the `--GOterm` flag can be set as F to extract the GO Molecular Function terms, P to extract the Biological Process terms, and C to extract the Cellular Component terms. For input files containing several Uniprot IDs, the `--countplot` flag can be used to create and store a countplot that displays the GO terms that appear for the given group of proteins as a PNG file. The `--GOnum` file can be used to obtain a list of the GO terms associated with each protein. 

Please use the `--help` flag for more information. 

This tools operates on Python 3.7.0 and runs from command line. 

## Detailed information on each flag 

```
usage: ProtAccess.py [-h] [-u UNIPROTID] [-i INPUTFILE] [-n] [-p] [-s ISOFORM]
                     [-l AA_LENGTH] [-g] [-t GOTERM] [-c] [-x GONUM]

Extract information from UniProt and NCBI websites based on UniProt ID or list
of UniProt IDs obtained from website.

optional arguments:
  -h, --help            show this help message and exit
  -u UNIPROTID, --uniprotID UNIPROTID
                        Single UniProt ID. Used to obtain information on one
                        protein only.
  -i INPUTFILE, --inputfile INPUTFILE
                        File with UniProt IDs downloaded from uniprot.org as
                        list file. Used to obtain information on a group of
                        proteins.
  -n, --ncbiSeq         Generate fasta file(s) with RNA trascript sequences.
  -p, --proteinSeq      Generate fasta file(s) with protein sequences. Must be
                        used with -s and -l.
  -s ISOFORM, --isoform ISOFORM
                        Used with -p and -l. Set as 'canonical' to obtain
                        canonical protein sequence only and 'all' to obtain
                        protein sequences of all isoforms.
  -l AA_LENGTH, --aa_length AA_LENGTH
                        Used with -p and -s. Set as 'y' to obtain length of
                        amino acid sequence of each isoform and 'n' if such
                        information is not needed.
  -g, --GO              Stores Gene Ontology terms from the UniProt page of
                        each protein. Needs flag -t to specify GOterm. Needs
                        flag -c OR -o to generate output
  -t GOTERM, --GOterm GOTERM
                        Used with -g. Set as F for GO Molecular Function, P
                        for GO Biological Process, or C for Cellular Component
  -c, --countplot       Used with -g and -t to create a countplot of GO terms
                        associated with proteins in list. Only used for list
                        of UniProt IDs.
  -x GONUM, --GOnum GONUM
                        Used with -g and -t. Accepts integers. When given a
                        number x, it gives the top x GO terms for each protein
                        from its UniProt page. If less than x GO terms are
                        found for a protein, it gives all GO terms
```

## Example file

The example file *ID_list.txt* found on this github page can be cloned and used to test this program. Here, the short flags will be used instead of the long flags for simplicity. Please refer to the above section for a description of each flag. 

The below commad would create a folder for each UniProt ID found in the file and store the NCBI RNA transcripts in it. It would also display the gene name and ID and the number of trascripts for this UniProt ID as extracted from NCBI. 

```python ProtAccess.py -i ID_list.txt –ncbiSeq``` 

The below script would add a fasta files to each folder with the canonical sequence of this UniProt ID. It would also display the length of the canonical isoform. This information is extracted from UniProt. 

```
python ProtAccess.py -i ID_list.txt -p -s canonical -l y
```
On the other hand the below command would creata a fasta file to each folder with all isoforms associated with each protein. 

```
python ProtAccess.py -i ID_list.txt -p -s all -l y
```

The following command can be used to obtain a countplot of all Gene Ontology Molecular function terms associated with the proteins in the given list. a png file will be stored in your directory. 

```
python ProtAccess.py -i ID_list.txt -g -t F -c
```

The `--x` command can be used to display the an x number of GO terms for each protein. 

```
python ProtAccess.py -i ID_list.txt -g -t F -x 5
```

All the above commands can be used on a single UniProt ID (Use `-u` to input ID instead of `-i` to input file) except the `-c` flag that can only be used on a set of UniProt IDs. 

```
python ProtAccess.py -u A5YKK6 -n -p -s all -l y -g -t F -x 4
```


