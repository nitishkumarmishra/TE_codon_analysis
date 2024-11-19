from Bio import SeqIO

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
    # Initialize a local codon count dictionary
    codon_counts = {codon: 0 for codon in codon_table.keys()}

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_counts:
            codon_counts[codon] += 1

    return codon_counts

# Read the FASTA file and save the output to a TSV file
fasta_file = "temp_CDS.fa"
output_file = "codon_counts.tsv"

with open(output_file, "w") as out_file:
    header = ["Sequence_ID"] + list(codon_table.keys())
    out_file.write("\t".join(header) + "\n")

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        codon_counts = count_codons(sequence)
        counts = [str(codon_counts[codon]) for codon in codon_table.keys()]
        out_file.write(f"{record.id}\t" + "\t".join(counts) + "\n")

print(f"Codon counts have been saved to {output_file}")

