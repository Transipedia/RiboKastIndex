import gzip

def readFastq(nomFastq):
    try:
        with gzip.open(nomFastq, 'rt') as file:
        #Read and decode the contents of fastq.gz
        #Need to modify this implementation if reading big fastq
            lines=file.read().splitlines() 
            return lines
    except IOError:
        print("Wrong file name")

def getSequences(fastq):#Take the input from read fastq and return a list containing only the sequence
    sequences=[]
    for i in range (1,len(fastq)-1,4):
        sequences.append(fastq[i])
    return sequences

def compareFastq(seq1,seq2):#Take a 2 list of sequences and compare hem to see if there is diffrence between the 2
    if len(seq1)!=len(seq2):
        print("length are not the same")
        return False
    seq1Sorted=sorted(seq1)
    seq2Sorted=sorted(seq2)
    for i in range(len(seq1)):
        if seq1Sorted[i]!=seq2Sorted[i]:
            print("Sequences are different on the line :" + str((i+1)*4-2))
            return False
    print("All sequences are the same")
    return True
