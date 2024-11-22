from Bio import SeqIO
import pandas as pd

# File paths
file_path = 'Polysome_Monosome_Davide_Like_topTable_Unt_Vs_TGFb.txt'
tsv_file = "m39_canonicalTranscriptsGeneName.txt"
fasta_file = "../m39_transcript_region_fasta/protein_coding_CDS.fasta"

# Read the TSV file
df1 = pd.read_csv(file_path, sep='\t')

# Read the given TSV file
df_genes = pd.read_csv(tsv_file, sep='\t', header=None, names=['transcript_id', 'gene_name', 'gene_id'])

# Strip any leading or trailing spaces from the column names
df_genes.columns = df_genes.columns.str.strip()

# Sort the DataFrame based on logFC column
df1_sorted = df1.sort_values(by='logFC', ascending=False)

# Filter the DataFrame for Up and Down regulated genes
df_up = df1_sorted[(df1_sorted['logFC'] > 0.585) & (df1_sorted['adj.P.Val'] < 0.05)]
df_down = df1_sorted[(df1_sorted['logFC'] < -0.585) & (df1_sorted['adj.P.Val'] < 0.05)]

# Find common gene names/transcript IDs in df_up and df_down
common_genes_up = df_up.index.intersection(df_genes['gene_name'].str.strip())
common_genes_down = df_down.index.intersection(df_genes['gene_name'].str.strip())

# Print the common gene names/transcript IDs
print("Common genes in Up-regulated genes:")
print(common_genes_up)
print("\nCommon genes in Down-regulated genes:")
print(common_genes_down)

# Print the transcript IDs of common genes
print("Transcript IDs of common genes in Up-regulated genes:")
print(df_genes[df_genes['gene_name'].str.strip().isin(common_genes_up)]['transcript_id'])
print("\nTranscript IDs of common genes in Down-regulated genes:")
print(df_genes[df_genes['gene_name'].str.strip().isin(common_genes_down)]['transcript_id'])

# Get the transcript IDs of common genes
common_transcripts_up = df_genes[df_genes['gene_name'].str.strip().isin(common_genes_up)]['transcript_id']
common_transcripts_down = df_genes[df_genes['gene_name'].str.strip().isin(common_genes_down)]['transcript_id']

# Read the FASTA file and extract sequences for common transcripts
common_sequences_up = []
common_sequences_down = []

for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id in common_transcripts_up.values:
        common_sequences_up.append(record)
    elif record.id in common_transcripts_down.values:
        common_sequences_down.append(record)

# Save the common sequences to separate FASTA files
SeqIO.write(common_sequences_up, "Polysome_Monosome_Davide_Unt_Vs_TGFb_transcript_up.fasta", "fasta")
SeqIO.write(common_sequences_down, "Polysome_Monosome_Davide_Unt_Vs_TGFb_transcript_down.fasta", "fasta")

print("Common transcripts FASTA files have been created: 'Polysome_Monosome_Davide_Unt_Vs_TGFb_transcript_up.fasta' and 'Polysome_Monosome_Davide_Unt_Vs_TGFb_transcript_down.fasta'")

