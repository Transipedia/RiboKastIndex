import csv

def processSequence(fastq):#Process only second line of each fastq and add it to a list
    sequences=[]  
    for i in range(1, len(fastq)-1, 4):
        sequences.append(fastq[i])
    return sequences

def processTab(tab):#Return a list containing the SSR and its length and phase when the phasing score>70%
    with open(tab, 'r', newline='') as file:
        reader = csv.reader(file, delimiter=' ')
        listSSR=[]
        for row in reader:
            ssr = row[0]
            length = int(row[1])
            phase = int(row[2])
            info=[ssr,length,phase]
            listSSR.append(info)
        return(listSSR)

def phaseFastq(fastqName,listSSR):#Generate a sublist containing the tupple with the phase and kmer size corresponding to the current fastq
    sublistSSR=[]
    for ssr in listSSR:
        print(fastqName)
        if fastqName.find(ssr[0]) != -1:
            print("find substring")
            sublistSSR.append(ssr)
    return sublistSSR
def phasedKmerCounting(sequences, kmerSize, sublistSSR):
    dictKmer = {}

    # Extraire la valeur minimum de phase dans sublistSSR
    min_phase = min(sub[2] for sub in sublistSSR)

    # Calculer les phases dynamiques basées sur la valeur minimum trouvée
    phase1 = min_phase
    phase2 = phase1 + 3
    phase3 = phase2 + 3

    for seq in sequences:
        for sub in sublistSSR:
            seq_id, length, phase = sub
            if int(length) == len(seq) and len(seq) >= kmerSize:
                # Déterminer la phase choisie en fonction des plages spécifiques
                if phase1 <= phase < phase2:
                    chosenPhase = phase1
                elif phase2 <= phase < phase3:
                    chosenPhase = phase2
                elif phase >= phase3:
                    chosenPhase = phase3
                else:
                    raise Exception(f"Phase {phase} pour {seq_id} ne rentre pas dans les plages spécifiées")

                phaseCut = int(phase) - chosenPhase
                if phaseCut < 0:
                    raise Exception("La phase souhaitée est inférieure à la phase de la séquence kmer actuelle")

                endCut = len(seq) - kmerSize - phaseCut
                if endCut == 0:  # Si pas de coupure, prendre le kmer entier
                    endCut = -len(seq)

                kmer = seq[phaseCut:-endCut]
                if kmer in dictKmer:  # Ajouter le kmer ou incrémenter son nombre d'occurrences
                    dictKmer[kmer] += 1
                else:
                    dictKmer[kmer] = 1

    sorted_dictKmer = dict(sorted(dictKmer.items()))
    return sorted_dictKmer

def kmerCounting(sequences,kmerSize):
    dictKmer={}
    for seq in sequences:
        kmer=[]
        if len(seq)>=kmerSize:#Check if seq is equal or greater than the wanted kmer size
            for i in range(len(seq)-kmerSize+1):#Formula for all the possible kmer
                kmer=seq[i:i+kmerSize]
                if kmer in dictKmer:#Add the kmer in the dict or increment his number of occurence
                    dictKmer[kmer]+=1
                else:
                    dictKmer[kmer]=1
    sorted_dictKmer = dict(sorted(dictKmer.items()))
    return sorted_dictKmer

def createTSV(filename,dico):#Write the Kmer file in a tsv format
    with open(filename,"w") as file:
        for key in dico:
            file.write(key + "\t" + str(dico[key])+ "\n")
