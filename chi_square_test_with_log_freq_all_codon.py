from Bio import SeqIO
from scipy.stats import chi2_contingency
import pandas as pd
import numpy as np

"""
This script performs a Chi-square test on codon counts from two FASTA files, calculates codon frequencies, computes the log ratio of these frequencies, and saves the results to a TSV file. The key steps are as follows:

1. Import Libraries:
   - SeqIO from Biopython for parsing FASTA files.
   - chi2_contingency from SciPy for performing the Chi-square test.
   - pandas for data manipulation and storage.
   - numpy for numerical operations.

2. Define Codon Table:
   - A dictionary `codon_table` is defined to store codon counts.

3. Functions:
   - `count_codons(sequence)`: Counts the occurrences of each codon in a given sequence.
   - `get_codon_counts_from_fasta(fasta_file)`: Parses a FASTA file and returns the total codon counts.
   - `calculate_codon_frequency(codon_counts)`: Calculates the frequency of each codon based on the total counts.

4. Read FASTA Files:
   - Two FASTA files (`temp_CDS.fa` and `temp_CDS1.fa`) are read to obtain codon counts using the `get_codon_counts_from_fasta` function.

5. Calculate Codon Frequencies:
   - Codon frequencies are calculated for each file using the `calculate_codon_frequency` function.

6. Calculate Log Ratio of Codon Frequencies:
   - The log ratio of codon frequencies between the two files is computed.

7. Create Contingency Table:
   - A contingency table is created using the codon counts from both files.

8. Perform Chi-square Test:
   - A global Chi-square test is performed on the contingency table.
   - The Chi-square statistic, p-value, degrees of freedom, and expected frequencies are printed, with expected frequencies rounded to three decimal places.

9. Save Results to File:
   - The results (observed and expected counts) are saved to a TSV file (`chi_square_results.tsv`).

10. Calculate Chi-square and P-values for Each Codon:
    - For each codon, a separate Chi-square test is performed.
    - The log frequency ratio is calculated.
    - The results, including rounded p-values, frequencies, and log frequency ratios, are appended to a list.

11. Save Detailed Results to File:
    - The detailed results, including p-values and log frequency ratios, are saved to another TSV file (`chi_square_p_values_with_log_ratio.tsv`).

12. Print Significant Results:
    - Codons with p-values less than 0.05 are printed.
"""


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

def calculate_codon_frequency(codon_counts):
    total_codons = sum(codon_counts.values())
    codon_frequency = {codon: count / total_codons for codon, count in codon_counts.items()}
    return codon_frequency

# Read the two FASTA files
fasta_file1 = "temp_CDS.fa"
fasta_file2 = "temp_CDS1.fa"

codon_counts1 = get_codon_counts_from_fasta(fasta_file1)
codon_counts2 = get_codon_counts_from_fasta(fasta_file2)

# Calculate codon frequencies
codon_frequency1 = calculate_codon_frequency(codon_counts1)
codon_frequency2 = calculate_codon_frequency(codon_counts2)

# Calculate log ratio of codon frequencies
log_ratio = {codon: np.log2(codon_frequency1[codon] / codon_frequency2[codon]) if codon_frequency2[codon] != 0 else np.inf for codon in codon_table.keys()}

# Create a contingency table
contingency_table = np.array([
    [codon_counts1[codon] for codon in codon_table.keys()],
    [codon_counts2[codon] for codon in codon_table.keys()]
])

# Perform chi-square test
chi2, p_global, dof_global, expected_global = chi2_contingency(contingency_table)

# Print results
print(f"Chi-square statistic: {chi2:.3f}")
print(f"P-value: {p_global:.3f}")
print(f"Degrees of freedom: {dof_global}")
print(f"\n")
print(f"Expected frequencies: \n{np.round(expected_global, 3)}")

# Save results to a file
output_file = "chi_square_results.tsv"
with open(output_file, "w") as out_file:
    out_file.write("Codon\tObserved_File1\tObserved_File2\tExpected_File1\tExpected_File2\n")
    for i, codon in enumerate(codon_table.keys()):
        out_file.write(f"{codon}\t{contingency_table[0][i]}\t{contingency_table[1][i]}\t{expected_global[0][i]}\t{expected_global[1][i]}\n")

print(f"Chi-square test results have been saved to :: {output_file}")

# Read the input file
df = pd.read_csv(output_file, sep='\t')

# Calculate chi-square and p-values for each row separately and add log frequency ratio
results = []
for index, row in df.iterrows():
    observed = [[row['Observed_File1']], [row['Observed_File2']]]
    expected = [[row['Expected_File1']], [row['Expected_File2']]]
    chi2_row, p_row, dof_row, exp_row = chi2_contingency([observed, expected])
    
    # Calculate log frequency ratio
    freq1 = row['Observed_File1'] / sum(df['Observed_File1'])
    freq2 = row['Observed_File2'] / sum(df['Observed_File2'])
    log_freq_ratio = np.log2(freq1 / freq2) if freq2 != 0 else np.inf
    
    #results.append((row['Codon'], row['Observed_File1'], row['Observed_File2'], row['Expected_File1'], row['Expected_File2'], p_row, round(freq1, 3), round(freq2, 3), round(log_freq_ratio, 3)))a
    results.append((
        row['Codon'],
        row['Observed_File1'],
        row['Observed_File2'],
        round(row['Expected_File1'], 3),
        round(row['Expected_File2'], 3),
        round(p_row, 3),
        round(freq1, 3),
        round(freq2, 3),
        round(log_freq_ratio, 3)
    ))

# Create a DataFrame for the results
results_df = pd.DataFrame(results, columns=['Codon', 'Observed1', 'Observed2', 'Expected1', 'Expected2', 'P-value', 'Freq1', 'Freq2', 'Log_Freq_Ratio'])

# Save the results to a TSV file
output_file_with_p_values_and_log_ratio = "chi_square_p_values_with_log_ratio.tsv"
results_df.to_csv(output_file_with_p_values_and_log_ratio, sep='\t', index=False)

print(f"Chi-square test results with p-values and log frequency ratios have been saved to :: {output_file_with_p_values_and_log_ratio}")

# Print lines with p-value < 0.05
significant_results_df = results_df[results_df['P-value'] < 0.05]
print(f"\n")
print(significant_results_df)

