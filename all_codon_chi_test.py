from Bio import SeqIO
from scipy.stats import chi2_contingency

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

#fasta_file1 = "random_1000_first_5000_sequences.fasta"
#fasta_file2 = "random_1000_last_5000_sequences.fasta"



codon_counts1 = get_codon_counts_from_fasta(fasta_file1)
codon_counts2 = get_codon_counts_from_fasta(fasta_file2)



# Print codon counts for debugging
print("Codon counts for file 1:")
print(codon_counts1)
print("Codon counts for file 2:")
print(codon_counts2)


# Perform chi-square test for each codon
results = []
for codon in codon_table.keys():
    observed = [[codon_counts1[codon]], [codon_counts2[codon]]]
    chi2, p, dof, expected = chi2_contingency(observed)
    results.append((codon, chi2, p))
print(f"Observed: {observed}")
print(f"Expected : {expected}")
print(f"p-value : {p}")
# Print results
for codon, chi2, p in results:
    print(f"Codon: {codon}, Chi-square: {chi2}, P-value: {p}")
    if p < 0.05:
        print(f"The codon {codon} usage is significantly different between the two FASTA files.")
    else:
        print(f"The codon {codon} usage is not significantly different between the two FASTA files.")


# Perform chi-square test for each codon and save results
output_file = "chi_square_results.tsv"
with open(output_file, "w") as out_file:
    out_file.write("Codon\tChi-square\tP-value\n")
    for codon in codon_table.keys():
        observed = [[codon_counts1[codon]], [codon_counts2[codon]]]
        chi2, p, dof, expected = chi2_contingency(observed)
        out_file.write(f"{codon}\t{chi2}\t{p}\n")

print(f"Chi-square test results have been saved to {output_file}")


