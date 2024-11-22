from Bio import SeqIO
from scipy.stats import chi2_contingency
import pandas as pd
import numpy as np

# Define the codon table globally
codon_table = {
    'ATA': 0, 'ATC': 0, 'ATT': 0, 'ATG': 0,
    'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0,
    'AAC': 0, 'AAT': 0, 'AAA': 0, 'AAG': 0,
    'AGC': 0, 'AGT': 0, 'AGA': 0, 'AGG': 0,
    'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0,
    'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0,
    'CAC': 0, 'CAT': 0, 'CAA': 0, 'CAG': 0,
    'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0,
    'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0,
    'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0,
    'GAC': 0, 'GAT': 0, 'GAA': 0, 'GAG': 0,
    'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0,
    'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0,
    'TTC': 0, 'TTT': 0, 'TTA': 0, 'TTG': 0,
    'TAC': 0, 'TAT': 0, 'TAA': 0, 'TAG': 0,
    'TGC': 0, 'TGT': 0, 'TGA': 0, 'TGG': 0,
}

def count_codons(sequence):
    codon_counts = {codon: 0 for codon in codon_table.keys()}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_counts:
            codon_counts[codon] += 1
    return codon_counts

def get_codon_counts_from_fasta(fasta_file):
    total_counts = {codon: 0 for codon in codon_table.keys()}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        codon_counts = count_codons(sequence)
        for codon, count in codon_counts.items():
            total_counts[codon] += count
    return total_counts

def process_files(fasta_file1_path, fasta_file2_path, output_file):
    # Get codon counts from the two FASTA files
    codon_counts1 = get_codon_counts_from_fasta(fasta_file1_path)
    codon_counts2 = get_codon_counts_from_fasta(fasta_file2_path)

    # Create a contingency table
    contingency_table = np.array([
        [codon_counts1[codon] for codon in codon_table.keys()],
        [codon_counts2[codon] for codon in codon_table.keys()]
    ])

    # Check for zero counts to avoid division by zero
    if np.any(contingency_table == 0):
        print("Warning: Zero counts found in contingency table. Adjusting to avoid division by zero.")
        contingency_table = contingency_table + np.finfo(float).eps

    # Perform chi-square test
    chi2, p_global, dof_global, expected_global = chi2_contingency(contingency_table)

    # Print results
    print(f"Global Chi-square statistic: {round(chi2, 3)}")
    print(f"Global P-value: {round(p_global, 5)}")
    print(f"Degrees of freedom: {dof_global}")

    # Save results to a file
    with open(output_file, "w") as out_file:
        out_file.write("Codon\tObserved_File1\tObserved_File2\tExpected_File1\tExpected_File2\n")
        for i, codon in enumerate(codon_table.keys()):
            out_file.write(f"{codon}\t{round(contingency_table[0][i], 3)}\t{round(contingency_table[1][i], 3)}\t{round(expected_global[0][i], 3)}\t{round(expected_global[1][i], 3)}\n")

    print(f"Chi-square test results have been saved to {output_file}")

    # Read the input file
    df = pd.read_csv(output_file, sep='\t')

    # Calculate chi-square and p-values for each row separately
    results = []
    for index, row in df.iterrows():
        observed = [[row['Observed_File1']], [row['Observed_File2']]]
        expected = [[row['Expected_File1']], [row['Expected_File2']]]
        chi2_row, p_row, dof_row, exp_row = chi2_contingency([observed, expected])
        p_value_str = f"{p_row}" + ("*" if p_row < 0.05 else "")
        results.append((row['Codon'], round(row['Observed_File1'], 3), round(row['Observed_File2'], 3), round(row['Expected_File1'], 3), round(row['Expected_File2'], 3), p_value_str))

    # Create a DataFrame for the results
    results_df = pd.DataFrame(results, columns=['Codon', 'Observed_File1', 'Observed_File2', 'Expected_File1', 'Expected_File2', '*P-value'])

    # Save the results to a TSV file
    output_file_with_p_values = output_file.replace(".tsv", "_with_p_values.tsv")
    results_df.to_csv(output_file_with_p_values, sep='\t', index=False)

    print(f"Chi-square test results with p-values have been saved to {output_file_with_p_values}")

# Input and output files
fasta_files1 = [
    'Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta',
    'Riboseq_Total_RNAseq_Davide_Like_topTable_CX5461_Vs_TGFb_up.fasta',
    'Polysome_Monosome_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta',
    'Polysome_Monosome_Davide_Like_topTable_CX5461_Vs_TGFb_up.fasta',
    'Polysome_Input_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta',
    'Polysome_Input_Davide_Like_topTable_CX5461_Vs_TGFb_up.fasta'
]
fasta_files2 = [
    'Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta',
    'Riboseq_Total_RNAseq_Davide_Like_topTable_CX5461_Vs_TGFb_down.fasta',
    'Polysome_Monosome_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta',
    'Polysome_Monosome_Davide_Like_topTable_CX5461_Vs_TGFb_down.fasta',
    'Polysome_Input_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta',
    'Polysome_Input_Davide_Like_topTable_CX5461_Vs_TGFb_down.fasta'
]
output_files = [
    'chi_square_results_Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb.tsv',
    'chi_square_results_Riboseq_Total_RNAseq_Davide_Like_topTable_CX5461_Vs_TGFb.tsv',
    'chi_square_results_Polysome_Monosome_Davide_Like_topTable_Unt_Vs_TGFb.tsv',
    'chi_square_results_Polysome_Monosome_Davide_Like_topTable_CX5461_Vs_TGFb.tsv',
    'chi_square_results_Polysome_Input_Davide_Like_topTable_Unt_Vs_TGFb.tsv',
    'chi_square_results_Polysome_Input_Davide_Like_topTable_CX5461_Vs_TGFb.tsv'
]

# Process each pair of files
for fasta_file1_path, fasta_file2_path, output_file in zip(fasta_files1, fasta_files2, output_files):
    process_files(fasta_file1_path, fasta_file2_path, output_file)

