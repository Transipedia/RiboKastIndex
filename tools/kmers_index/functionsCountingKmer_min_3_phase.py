import csv
from collections import defaultdict

def processSequence(fastq_lines):
    """
    Extract only the sequence lines (2nd line in each 4-line FASTQ block).
    Strips trailing newlines.
    """
    sequences = []
    for i in range(1, len(fastq_lines) - 1, 4):
        sequences.append(fastq_lines[i].strip())
    return sequences

def processTab(tab):
    """
    Read a space-separated table:
      <SSR_id> <length> <phase>
    Ignores empty lines and those starting with '#'.
    """
    listSSR = []
    with open(tab, 'r', newline='') as file:
        for raw in file:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            # Robust split on any whitespace (one or more spaces, tabs, etc.)
            parts = line.split()
            if len(parts) < 3:
                continue
            ssr = parts[0]
            length = int(parts[1])
            phase = int(parts[2])
            listSSR.append([ssr, length, phase])
    return listSSR

def phaseFastq(fastqName, listSSR):
    """
    Build a sublist of SSR rows whose SSR_id is contained in 'fastqName'.
    """
    sublistSSR = []
    for ssr in listSSR:
        if fastqName.find(ssr[0]) != -1:
            sublistSSR.append(ssr)
    return sublistSSR

def phasedKmerCounting(sequences, kmerSize, sublistSSR, forced_phase=None):
    """
    Count phased k-mers.
    - If forced_phase is None: use original bucket logic based on min phase:
        chosenPhase in {min_phase, min_phase+3, min_phase+6} depending on SSR phase.
    - If forced_phase is set: impose that chosenPhase for all SSRs (no buckets).
    Extract exactly ONE k-mer per matching read: seq[phaseCut : phaseCut + kmerSize]
    when 0 <= phaseCut and phaseCut + kmerSize <= len(seq).
    """
    if not sublistSSR:
        return {}

    dictKmer = defaultdict(int)

    if forced_phase is None:
        # Original dynamic buckets (min, min+3, min+6)
        min_phase = min(sub[2] for sub in sublistSSR)
        phase1 = min_phase
        phase2 = phase1 + 3
        phase3 = phase2 + 3

    for seq in sequences:
        L = len(seq)
        if L < kmerSize:
            continue

        for ssr_id, length, phase in sublistSSR:
            if int(length) != L:
                continue

            if forced_phase is None:
                # Bucket assignment
                if phase1 <= phase < phase2:
                    chosenPhase = phase1
                elif phase2 <= phase < phase3:
                    chosenPhase = phase2
                elif phase >= phase3:
                    chosenPhase = phase3
                else:
                    # This shouldn't happen if buckets cover all >= min_phase
                    continue
            else:
                # Imposed phase: exactly what the user asked for
                chosenPhase = forced_phase

            phaseCut = int(phase) - chosenPhase
            # We only accept non-negative cut and the slice must fit entirely
            if phaseCut < 0:
                continue
            end_pos = phaseCut + kmerSize
            if end_pos > L:
                continue

            kmer = seq[phaseCut:end_pos]
            dictKmer[kmer] += 1

    return dict(sorted(dictKmer.items()))

def kmerCounting(sequences, kmerSize):
    dictKmer = defaultdict(int)
    for seq in sequences:
        L = len(seq)
        if L >= kmerSize:
            for i in range(L - kmerSize + 1):
                kmer = seq[i:i+kmerSize]
                dictKmer[kmer] += 1
    return dict(sorted(dictKmer.items()))

def createTSV(filename, dico):
    with open(filename, "w") as file:
        for key, val in dico.items():
            file.write(f"{key}\t{val}\n")
