import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    current_tag = None
    current_sum = []
    current_count = 0

    for line in f_in:
        parts = line.strip().split('\t')
        tag = parts[0]
        values = list(map(int, parts[1:]))

        # New tag encountered
        if tag != current_tag:
            # If we were already accumulating a previous tag, write its result
            if current_tag is not None:
                f_out.write(f"{current_tag}\t" + "\t".join(map(str, current_sum)) + "\n")

            # Start accumulation for new tag
            current_tag = tag
            current_sum = values
            current_count = 1
        else:
            # Same tag — accumulate the values
            current_sum = [a + b for a, b in zip(current_sum, values)]
            current_count += 1

    # Write the last accumulated tag
    if current_tag is not None:
        f_out.write(f"{current_tag}\t" + "\t".join(map(str, current_sum)) + "\n")
