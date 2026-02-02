import pandas as pd
import sys

def filter_and_merge(input_file, merge_file, output_file, final_output_file, pct_threshold):
    # Read the TSV file containing p-sites per region
    df = pd.read_csv(input_file, sep='\t')

    # Filter rows where region == "CDS" and percentage >= threshold
    filtered_df = df[(df['region'] == 'CDS') & (df['percentage'] >= pct_threshold)]

    # Read the second TSV file containing offsets
    df2 = pd.read_csv(merge_file, sep='\t')

    # Merge on ["length", "sample"]
    merged_df = pd.merge(filtered_df, df2, on=['length', 'sample'])

    # Save the intermediate merged table
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"[INFO] Intermediate file saved: {output_file}")

    # Build the final output for KmerCount
    final_output = merged_df[['sample', 'length', 'corrected_offset_from_5']].copy()

    # Optional: clean sample name (keep part after the first dot)
    final_output['sample'] = final_output['sample'].apply(
        lambda x: x.split('.')[1] if '.' in x and len(x.split('.')) > 1 else x
    )

    # Save final output without header, space-separated
    final_output.to_csv(final_output_file, sep=' ', index=False, header=False)
    print(f"[INFO] Final KmerCount file generated: {final_output_file}")

if __name__ == "__main__":
    # Expected args: script.py input merge output final_output pct_threshold
    if len(sys.argv) != 6:
        print("Usage: python3 filter_merge_for_KmerCount.py <input_file> <merge_file> <output_file> <final_output_file> <pct_threshold>")
        sys.exit(1)

    # Parse threshold as a number (allows 60 or 60.0)
    try:
        pct_threshold = float(sys.argv[5])
    except ValueError:
        print("[ERROR] pct_threshold must be a number (e.g., 60 or 60.0).")
        sys.exit(1)

    filter_and_merge(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], pct_threshold)
