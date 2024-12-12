from Bio import SeqIO
from collections import defaultdict

def find_headers_with_stop_codons(fasta_file):
    stop_codons = {"TAA", "TAG", "TGA"}
    headers_with_stop_codons = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        for i in range(0, len(sequence) - 5, 3):
            codon = sequence[i:i+3]
            if codon in stop_codons:
                headers_with_stop_codons.append(record.id)
                break  # Stop checking further once a stop codon is found

    return headers_with_stop_codons

def write_combined_headers_to_tsv(fasta_files, output_file):
    with open(output_file, 'w') as f:
        for fasta_file, label in fasta_files:
            headers = find_headers_with_stop_codons(fasta_file)
            for header in headers:
                f.write(f"{label}\t{header}\n")

# Paths to your FASTA files with labels
fasta_files = [
    ("Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta", "up"),
    ("Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta", "down"),
    ("../m39_transcript_region_fasta/protein_coding_CDS.fasta", "total")
]

# Write combined headers with in-frame stop codons to a single TSV file
write_combined_headers_to_tsv(fasta_files, "combined_headers_with_stop_codons.tsv")

print("Combined headers with in-frame stop codons have been written to combined_headers_with_stop_codons.tsv file.")

