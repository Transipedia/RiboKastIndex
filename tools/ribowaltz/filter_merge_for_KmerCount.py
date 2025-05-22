import pandas as pd
import sys

def filter_and_merge(input_file, merge_file, output_file, final_output_file):
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

    # Générer le fichier final avec les colonnes "sample", "length", et "corrected_offset_from_5"
    final_output = merged_df[['sample', 'length', 'corrected_offset_from_5']].copy()

    # Extraire la partie du sample entre les deux points
    final_output['sample'] = final_output['sample'].apply(lambda x: x.split('.')[1])

    # Enregistrer le résultat dans le fichier final avec des séparateurs d'espace
    final_output.to_csv(final_output_file, sep=' ', index=False, header=False)

    print(f"Fichier final généré avec les colonnes souhaitées : {final_output_file}")

if __name__ == "__main__":
    # Vérification des arguments
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_file> <merge_file> <output_file> <final_output_file>")
    else:
        input_file = sys.argv[1]
        merge_file = sys.argv[2]
        output_file = sys.argv[3]
        final_output_file = sys.argv[4]
        filter_and_merge(input_file, merge_file, output_file, final_output_file)
