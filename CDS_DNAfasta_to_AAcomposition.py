from Bio import SeqIO
from collections import Counter


# Function to convert CDS FASTA sequences to protein sequences
def convert_cds_fasta_to_protein(input_fasta, output_fasta):
    # Open the input and output FASTA files
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Translate the CDS sequence to protein sequence
            protein_sequence = record.seq.translate()
            # Create a new SeqRecord for the protein sequence
            protein_record = record[:]
            protein_record.seq = protein_sequence
            # Write the protein sequence to the output FASTA file
            SeqIO.write(protein_record, output_handle, "fasta")


# Function to count all 20 amino acids in the CDS FASTA file and write counts to a separate file
def sum_amino_acids_in_fasta(input_fasta, output_file):
    amino_acid_counter = Counter()
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # All 20 amino acids
    
    # Open the input FASTA file
    with open(input_fasta, "r") as input_handle:
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Translate the CDS sequence to protein sequence
            protein_sequence = str(record.seq.translate())
            # Update the counter with the amino acids in the protein sequence
            amino_acid_counter.update(protein_sequence)
    
    # Write the amino acid counts to the output file
    with open(output_file, "w") as output_handle:
        for amino_acid in amino_acids:
            count = amino_acid_counter.get(amino_acid, 0)
            output_handle.write(f"Amino Acid: {amino_acid}, Count: {count}\n")


def avg_amino_acids_in_fasta(input_fasta, output_file):
    amino_acid_counter = Counter()
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # All 20 amino acids
    total_proteins = 0

    # Open the input FASTA file
    with open(input_fasta, "r") as input_handle:
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            total_proteins += 1
            # Translate the CDS sequence to protein sequence
            protein_sequence = str(record.seq.translate())
            # Update the counter with the amino acids in the protein sequence
            amino_acid_counter.update(protein_sequence)

    # Write the average amino acid counts to the output file
    with open(output_file, "w") as output_handle:
        for amino_acid in amino_acids:
            avg_count = amino_acid_counter.get(amino_acid, 0) / total_proteins
            output_handle.write(f"{amino_acid}, Average Count: {avg_count:.2f}\n")


# Function to count all 20 amino acids in each protein sequence in the CDS FASTA file and write counts to a separate file
def count_amino_acids_in_fasta(input_fasta, output_file):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # All 20 amino acids
    
    # Open the input FASTA file
    with open(input_fasta, "r") as input_handle, open(output_file, "w") as output_handle:
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Translate the CDS sequence to protein sequence
            protein_sequence = str(record.seq.translate())
            # Count the amino acids in the protein sequence
            amino_acid_counter = Counter(protein_sequence)
            
            # Write the protein ID and amino acid counts to the output file
            output_handle.write(f"Protein ID: {record.id}\n")
            for amino_acid in amino_acids:
                count = amino_acid_counter.get(amino_acid, 0)
                output_handle.write(f"{amino_acid}, Count: {count}\n")
            output_handle.write("\n")


# Amino acid counts in matrix format
def count_amino_acids_matrix(input_fasta, output_file):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # All 20 amino acids

    # Open the input FASTA file
    with open(input_fasta, "r") as input_handle:
        # Initialize a list to store the counts for each protein
        protein_counts = []
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Translate the CDS sequence to protein sequence
            protein_sequence = str(record.seq.translate())
            # Count the amino acids in the protein sequence
            amino_acid_counter = Counter(protein_sequence)
            # Create a list of counts for the current protein
            counts = [record.id] + [str(amino_acid_counter.get(amino_acid, 0)) for amino_acid in amino_acids]
            # Append the counts to the protein_counts list
            protein_counts.append(counts)

    # Write the matrix to the output file
    with open(output_file, "w") as output_handle:
        # Write the header row
        output_handle.write("Transcript_ID\t" + " ".join(amino_acids) + "\n")
        # Write the counts for each protein
        for counts in protein_counts:
            output_handle.write(" ".join(counts) + "\n")


# Function to calculate the percentage of each amino acid in the CDS FASTA file and write counts to a matrix file
def percentage_amino_acids_matrix(input_fasta, output_file):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # All 20 amino acids
    total_amino_acids = Counter()
    
    # Open the input FASTA file
    with open(input_fasta, "r") as input_handle:
        # Iterate over each record in the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Translate the CDS sequence to protein sequence
            protein_sequence = str(record.seq.translate())
            # Update the counter with the amino acids in the protein sequence
            total_amino_acids.update(protein_sequence)
    
    # Calculate the total number of amino acids
    total_count = sum(total_amino_acids.values())
    
    # Write the percentage amino acid counts to the output file
    with open(output_file, "w") as output_handle:
        output_handle.write("Amino_Acid\tPercentage\n")
        for amino_acid in amino_acids:
            percentage = (total_amino_acids.get(amino_acid, 0) / total_count) * 100
            output_handle.write(f"{amino_acid}\t{percentage:.2f}\n")


# Example usage
up_fasta = "Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta"
down_fasta = "Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta"

output_up_fasta = "protein_sequences_Riboseq_Unt_Vs_TGFb_up.fasta"
output_down_fasta = "protein_sequences_Riboseq_Unt_Vs_TGFb_down.fasta"

output_up_count = "protein_sequences_Riboseq_Unt_Vs_TGFb_up.txt"
output_down_count = "protein_sequences_Riboseq_Unt_Vs_TGFb_down.txt"

output_up_count_matrix = "protein_sequences_Riboseq_Unt_Vs_TGFb_up_matrix.txt"
output_down_count_matrix = "protein_sequences_Riboseq_Unt_Vs_TGFb_down_matrix.txt"

output_up_sum = "protein_sequences_Riboseq_Unt_Vs_TGFb_up_sum.txt"
output_down_sum = "protein_sequences_Riboseq_Unt_Vs_TGFb_down_sum.txt"


output_up_avg = "protein_sequences_Riboseq_Unt_Vs_TGFb_up_avg.txt"
output_down_avg = "protein_sequences_Riboseq_Unt_Vs_TGFb_down_avg.txt"

output_up_percent = "protein_sequences_Riboseq_Unt_Vs_TGFb_up_percent.txt"
output_down_percent = "protein_sequences_Riboseq_Unt_Vs_TGFb_down_percent.txt"

# Convert CDS FASTA sequences to protein sequences
convert_cds_fasta_to_protein(up_fasta, output_up_fasta)
convert_cds_fasta_to_protein(down_fasta, output_down_fasta)


# Count each protein in the CDS FASTA file and write counts to a separate file
count_amino_acids_in_fasta(up_fasta, output_up_count)
count_amino_acids_in_fasta(down_fasta, output_down_count)

# Count each protein in the CDS FASTA file and write counts to a separate matrix file
count_amino_acids_matrix(up_fasta, output_up_count_matrix)
count_amino_acids_matrix(down_fasta, output_down_count_matrix)

# Sum each protein in the CDS FASTA file and write sum to a separate file
sum_amino_acids_in_fasta(up_fasta, output_up_sum)
sum_amino_acids_in_fasta(down_fasta, output_down_sum)


# Average each protein in the CDS FASTA file and write average to a separate file
avg_amino_acids_in_fasta(up_fasta, output_up_avg)
avg_amino_acids_in_fasta(down_fasta, output_down_avg)

# Percent each protein in the CDS FASTA file and write average to a separate file
percentage_amino_acids_matrix(up_fasta, output_up_percent)
percentage_amino_acids_matrix(down_fasta, output_down_percent)

print(f"Protein sequences have been written to {output_up_fasta}")
print(f"Protein sequences have been written to {output_down_fasta}")

print(f"Protein counts have been written to {output_up_count}")
print(f"Protein counts have been written to {output_down_count}")

print(f"Protein counts have been written to {output_up_count_matrix}")
print(f"Protein counts have been written to {output_down_count_matrix}")

print(f"Protein counts have been written to {output_up_sum}")
print(f"Protein counts have been written to {output_down_sum}")

print(f"Protein counts have been written to {output_up_avg}")
print(f"Protein counts have been written to {output_down_avg}")


print(f"Protein counts have been written to {output_up_percent}")
print(f"Protein counts have been written to {output_down_percent}")

