import csv
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="-i Path and file name of the input", 
                    type=str, required=True)
parser.add_argument("-o", "--output", help="-o Path and file name of the output", 
                    type=str, required=True)
args = parser.parse_args()

input_filename = args.input
output_filename = args.output

# Verify if the sequence is valid
def is_valid_sequence(sequence):
    valid_letters = set('ATGC')
    # Verifyif the sequence begins with A,T,C,G
    return sequence and sequence[0] in valid_letters and all(letter in valid_letters for letter in sequence)

# Remove unvalid sequences
def remove_non_atgc_lines(input_file, output_file):
    with open(input_file, 'r', newline='') as input_file, open(output_file, 'w', newline='') as output_file:
        reader = csv.reader(input_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')
        for row in reader:
            sequence = row[0]
            if is_valid_sequence(sequence):
                writer.writerow(row)

remove_non_atgc_lines(input_filename, output_filename)

