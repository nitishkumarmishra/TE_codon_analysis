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

# Read the two FASTA files
fasta_file1 = "temp_CDS.fa"
fasta_file2 = "temp_CDS1.fa"

codon_counts1 = get_codon_counts_from_fasta(fasta_file1)
codon_counts2 = get_codon_counts_from_fasta(fasta_file2)

# Create a contingency table
contingency_table = np.array([
    [codon_counts1[codon] for codon in codon_table.keys()],
    [codon_counts2[codon] for codon in codon_table.keys()]
])

# Perform chi-square test
chi2, p, dof, expected = chi2_contingency(contingency_table)

# Print results
print(f"Global Chi-square statistic: {chi2}")
print(f"Global P-value: {p}")
print(f"Degrees of freedom: {dof}")
#print(f"Expected frequencies: \n{expected}")

# Save results to a file
# Output file name
output_file = "chi_square_results.tsv"


with open(output_file, "w") as out_file:
    out_file.write("Codon\tObserved_File1\tObserved_File2\tExpected_File1\tExpected_File2\n")
    for i, codon in enumerate(codon_table.keys()):
        out_file.write(f"{codon}\t{contingency_table[0][i]}\t{contingency_table[1][i]}\t{expected[0][i]}\t{expected[1][i]}\n")

print(f"Chi-square test results have been saved to {output_file}")

# Read the input file
df = pd.read_csv(output_file, sep='\t')

# Calculate chi-square and p-values for each row separately
results = []
for index, row in df.iterrows():
    observed = [[row['Observed_File1']], [row['Observed_File2']]]
    expected = [[row['Expected_File1']], [row['Expected_File2']]]
    chi2_row, p_row, dof_row, exp_row = chi2_contingency([observed, expected])
    results.append((row['Codon'], row['Observed_File1'], row['Observed_File2'], row['Expected_File1'], row['Expected_File2'], p_row))

# Create a DataFrame for the results
results_df = pd.DataFrame(results, columns=['Codon', 'Observed_File1', 'Observed_File2', 'Expected_File1', 'Expected_File2', 'P-value'])

# Save the results to a TSV file
output_file_with_p_values = "chi_square_p_values_with_results.tsv"
results_df.to_csv(output_file_with_p_values, sep='\t', index=False)

print(f"Chi-square test results with p-values have been saved to {output_file_with_p_values}")

# Print lines with p-value < 0.05
significant_results_df = results_df[results_df['P-value'] < 0.05]
print(f"These are the Codons with p-value < 0.05")
print(significant_results_df)
