import pandas as pd
import sys

def filter_and_merge(input_file, merge_file, output_file):
    # Lire le fichier CSV (input_file)
    df = pd.read_csv(input_file, sep='\t')  # Remplacez sep='\t' par sep=',' si le fichier est séparé par des virgules

    # Filtrer les lignes où "region" est "CDS" et "percentage" est supérieur à 50%
    filtered_df = df[(df['region'] == 'CDS') & (df['percentage'] > 50)]

    # Lire le second fichier CSV (merge_file)
    df2 = pd.read_csv(merge_file, sep='\t')  # Remplacez sep='\t' par sep=',' si nécessaire

    # Effectuer la fusion sur les colonnes "length" et "sample"
    merged_df = pd.merge(filtered_df, df2, on=['length', 'sample'])

    # Enregistrer les résultats dans un fichier CSV
    merged_df.to_csv(output_file, sep='\t', index=False)  # Remplacez sep='\t' par sep=',' si nécessaire

    print(f"Filtrage et fusion terminés. Les résultats sont enregistrés dans {output_file}")

if __name__ == "__main__":
    # Vérification des arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <merge_file> <output_file>")
    else:
        input_file = sys.argv[1]
        merge_file = sys.argv[2]
        output_file = sys.argv[3]
        filter_and_merge(input_file, merge_file, output_file)
