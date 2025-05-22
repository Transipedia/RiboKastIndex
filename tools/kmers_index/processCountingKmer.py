#!/usr/bin/env python
import argparse
import fastqTools
import functionsCountingKmer_min_3 as fc
import re

parser=argparse.ArgumentParser()

parser.add_argument("-i", "--input",help="-i Path and name of the fastq ", 
                    type=str, required=True)
parser.add_argument("-m", "--mode",help="-m Mode to count the kmer : normal (count all the kmer of the length k), phase (count the kmer of length k corresponding to the phase)", 
                    type=str, required=True)
parser.add_argument("-k", "--kmerSize",help="-k Indicate the size of the kmer", 
                    type=int, required=True)
parser.add_argument("-o", "--output",help="-o Location of the resulting file", 
                    type=str, required=True)
parser.add_argument("-n", "--name",help="-n name of the output file", 
                    type=str, required=True)
parser.add_argument("-t", "--tab",help="-t  File containing the information about the phase of the kmer separated by a  whitespace", 
                    type=str, required=False)

args=parser.parse_args()

if args.mode=="phase": #Get the list of all srr with atleast a length of read phased
    listPhasedSRR = fc.processTab(args.tab)

dictKmer={}

fastq=fastqTools.readFastq(args.input)
sequences=fc.processSequence(fastq)
nameFastq=re.search(r'/([^/]+)$', args.input) #Get the full name of the file
nameBeforeFirstDot=nameFastq.group(1).split('.')
    
if args.mode=="phase":
    print("Nom fichier : " + nameBeforeFirstDot[0])
    sublistPhasedSRR = fc.phaseFastq(args.name, listPhasedSRR)
    print("Sous liste :")
    print(sublistPhasedSRR)
    dictKmer=fc.phasedKmerCounting(sequences, args.kmerSize, sublistPhasedSRR )
          
elif args.mode=="normal":
        dictKmer=fc.kmerCounting(sequences,args.kmerSize)
        
if (dictKmer):
    fc.createTSV(args.output + "/Kmer/" + args.name + ".tsv",dictKmer)
