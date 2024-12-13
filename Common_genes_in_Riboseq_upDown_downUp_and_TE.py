import pandas as pd

# Function to read two tab-delimited text files and find common row names
def find_common_rownames(file1, file2):
    # Read the canonical genes file into a DataFrame
    tsv_file = "m39_canonicalTranscriptsGeneName.txt"
    canonical_genes = pd.read_csv(tsv_file, sep='\t', header=None, names=['transcript_id', 'gene_name', 'gene_id'])

    # Read the first file into a DataFrame and sort by 'logFC'
    df1 = pd.read_csv(file1, sep='\t', index_col=0)
    df1_sorted = df1.sort_values(by='logFC', ascending=False)

    # Filter the DataFrame based on absolute 'logFC' and 'adj.P.Val'
    df_filtered = df1_sorted[(abs(df1_sorted['logFC']) > 0.585) & (df1_sorted['adj.P.Val'] < 0.05)]

    # Read the second file into a DataFrame and select the last 8 columns
    df2 = pd.read_csv(file2, sep='\t', index_col=0)
    df2_last_8_columns = df2.iloc[:, -8:]

    # Filter rows where 'reversible.translation' column has 'upDown' or 'downUp' pattern only
    df2_filtered = df2_last_8_columns[df2_last_8_columns['reversible.translation'].str.contains('upDown|downUp', na=False)]
    df2_filtered.columns = df2_filtered.columns.str.replace('.translation|.DEtranslation', '', regex=True)

    # Rename columns by removing the '.translation' or '.DEtranslation' suffix if it exists
    df2_last_8_columns.columns = df2_last_8_columns.columns.str.replace('.translation|.DEtranslation', '', regex=True)

    # Filter rows where 'reversible.translation' column has 'upDown' or 'downUp' pattern only
    df2_filtered = df2_last_8_columns[df2_last_8_columns['reversible'].str.contains('upDown|downUp', na=False)]

    # Print the head of both DataFrames
    print("Head of file1:")
    print(df_filtered.head())
    print("\nHead of file2:")
    print(df2_filtered.head())

    # Find the common row names between df1 and df2
    common_rownames = df_filtered.index.intersection(df2_filtered.index)

    # Find the common row names between common_rownames and canonical_genes
    common_rownames = common_rownames.intersection(canonical_genes['gene_name'].str.strip())

    # Get the canonical genes that are in common_rownames
    canonical_de = canonical_genes[canonical_genes['gene_name'].str.strip().isin(common_rownames)]

    # Get the filtered df2 rows that are in common_rownames
    df2_canonical = df2_filtered[df2_filtered.index.str.strip().isin(common_rownames)]

    return canonical_de, df2_canonical

# Example usage
file1 = 'Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb.txt'
file2 = 'Riboseq_total_RNAseq_DEGs_protein_coding.txt'

# Find common row names and get the DataFrames
canonical_de, df2_canonical = find_common_rownames(file1, file2)

# Print the common row names
print("Common row namesa Riboseq_Unt_Vs_TGFb_TE_RNAseq_DEGs_upDown_downUp:")
print(canonical_de['gene_name'].tolist())

# Save the DataFrames to an Excel file with separate sheets
with pd.ExcelWriter('Riboseq_Unt_Vs_TGFb_TE_RNAseq_DEGs_upDown_downUp.xlsx') as writer:
    canonical_de.to_excel(writer, sheet_name='Canonical_DE_Genes')
    df2_canonical.to_excel(writer, sheet_name='DF2_Canonical')

print("Common genes in Riboseq upDown/downUp and DE TE genes and output saved in excel file :: Riboseq_Unt_Vs_TGFb_TE_RNAseq_DEGs_upDown_downUp.xlsx")


