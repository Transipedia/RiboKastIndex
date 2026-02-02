#!/usr/bin/env python3
import argparse
import os
import re

import fastqTools
import functionsCountingKmer_min_3_phase as fc  # contains phasedKmerCounting(..., forced_phase=None)

def main():
    parser = argparse.ArgumentParser(
        description="K-mer counting: normal or phased. Can optionally impose a fixed phase."
    )

    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Path to FASTQ file")
    parser.add_argument("-m", "--mode", type=str, required=True, choices=["normal", "phase"],
                        help="Counting mode: 'normal' (all k-mers) or 'phase' (phased k-mers)")
    parser.add_argument("-k", "--kmerSize", type=int, required=True,
                        help="K-mer size (k)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output folder (TSV will be written under <output>/Kmer/)")
    parser.add_argument("-n", "--name", type=str, required=True,
                        help="Base name for the output TSV (no extension)")
    parser.add_argument("-t", "--tab", type=str, required=False,
                        help="Table with: <SSR_id> <length> <phase> (space-separated)")

    # NEW: allow forcing a fixed phase (disables min/min+3/min+6 logic)
    parser.add_argument("--forced-phase", type=int, default=None,
                        help="Force a specific phase (integer). If unset, use the original min/min+3/min+6 buckets.")
    # Optional: let user provide a 'logical' FASTQ name for SSR filtering (fallback: file basename)
    parser.add_argument("--fastq-name", type=str, default=None,
                        help="Logical FASTQ name used to match SSR entries. Default: input basename without extension.")

    args = parser.parse_args()

    if args.mode == "phase" and not args.tab:
        parser.error('In "phase" mode, --tab is required.')

    # Read FASTQ (expects fastqTools.readFastq to return list of lines)
    fastq_lines = fastqTools.readFastq(args.input)
    sequences = fc.processSequence(fastq_lines)

    # Derive FASTQ logical name (if not provided)
    if args.fastq_name:
        fastq_name_key = args.fastq_name
    else:
        m = re.search(r'/([^/]+)$', args.input)
        base = m.group(1) if m else os.path.basename(args.input)
        fastq_name_key = base.split('.')[0]

    dictKmer = {}

    if args.mode == "phase":
        listPhasedSRR = fc.processTab(args.tab)
        sublistPhasedSRR = fc.phaseFastq(fastq_name_key, listPhasedSRR)

        dictKmer = fc.phasedKmerCounting(
            sequences=sequences,
            kmerSize=args.kmerSize,
            sublistSSR=sublistPhasedSRR,
            forced_phase=args.forced_phase  # NEW
        )

    elif args.mode == "normal":
        dictKmer = fc.kmerCounting(sequences, args.kmerSize)

    if dictKmer:
        out_dir = os.path.join(args.output, "Kmer")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{args.name}.tsv")
        fc.createTSV(out_path, dictKmer)
        print(f"[OK] Written: {out_path}")
    else:
        print("[INFO] No k-mers counted. Check inputs/SSR filtering.")

if __name__ == "__main__":
    main()
