from Bio import SeqIO
import pandas as pd

# Define the codon table globally in the specified order, including stop codons
codon_table = {
    'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0,
    'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0,
    'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0,
    'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0,
    'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0,
    'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0,
    'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0,
    'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0,
    'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0,
    'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0,
    'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0,
    'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0,
    'TAC': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0,
    'TCG': 0, 'TCT': 0, 'TGC': 0, 'TGG': 0,
    'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0,
    'TTT': 0,
    # Stop codons
    'TAA': 0, 'TAG': 0, 'TGA': 0
}

def count_codons(sequence):
    codon_counts = {codon: 0 for codon in codon_table.keys()}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_counts:
            codon_counts[codon] += 1
    return codon_counts

def get_codon_frequencies_for_each_sequence(fasta_file):
    sequence_frequencies = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        codon_counts = count_codons(sequence)
        total_codons = sum(codon_counts.values())
        codon_frequencies = {codon: round(count / total_codons, 4) for codon, count in codon_counts.items()}
        sequence_frequencies[record.id] = codon_frequencies

    return sequence_frequencies

def convert_to_codon_matrix(sequence_frequencies):
    # Create a DataFrame from the sequence frequencies dictionary
    df = pd.DataFrame.from_dict(sequence_frequencies, orient='index')

    # Fill NaN values with zeros (in case some codons are missing in some sequences)
    df.fillna(0, inplace=True)

    # Reorder columns according to the specified codon order
    df = df[codon_table.keys()]

    return df

def calculate_sum_of_frequencies(codon_matrix):
    # Calculate the sum of frequencies for each codon
    sum_frequencies = codon_matrix.sum(axis=0)

    return sum_frequencies

# Example usage
protein_coding = '../m39_transcript_region_fasta/protein_coding_CDS.fasta'
protein_coding_up = 'Polysome_Input_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta'
protein_coding_down = 'Polysome_Input_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta'

protein_coding_sequence_frequencies = get_codon_frequencies_for_each_sequence(protein_coding)
protein_coding_codon_matrix = convert_to_codon_matrix(protein_coding_sequence_frequencies)

protein_coding_up_sequence_frequencies = get_codon_frequencies_for_each_sequence(protein_coding_up)
protein_coding_up_codon_matrix = convert_to_codon_matrix(protein_coding_up_sequence_frequencies)

protein_coding_down_sequence_frequencies = get_codon_frequencies_for_each_sequence(protein_coding_down)
protein_coding_down_codon_matrix = convert_to_codon_matrix(protein_coding_down_sequence_frequencies)

# Calculate the sum of frequencies for each codon
sum_frequencies_protein_coding = calculate_sum_of_frequencies(protein_coding_codon_matrix)
sum_frequencies_protein_coding_up = calculate_sum_of_frequencies(protein_coding_up_codon_matrix)
sum_frequencies_protein_coding_down = calculate_sum_of_frequencies(protein_coding_down_codon_matrix)

# Write all outputs to an Excel file with different sheets
with pd.ExcelWriter('codon_analysis_results.xlsx') as writer:
    protein_coding_codon_matrix.to_excel(writer, sheet_name='Protein Coding')
    protein_coding_up_codon_matrix.to_excel(writer, sheet_name='Upregulated')
    protein_coding_down_codon_matrix.to_excel(writer, sheet_name='Downregulated')
    sum_frequencies_protein_coding.to_excel(writer, sheet_name='Sum Frequencies Protein Coding')
    sum_frequencies_protein_coding_up.to_excel(writer, sheet_name='Sum Frequencies Upregulated')
    sum_frequencies_protein_coding_down.to_excel(writer, sheet_name='Sum Frequencies Downregulated')

print("All outputs have been written to 'codon_analysis_results.xlsx'")

