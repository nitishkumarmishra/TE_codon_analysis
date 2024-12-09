from Bio import SeqIO
from collections import defaultdict
import scipy.stats as stats
import math
from tabulate import tabulate
import matplotlib.pyplot as plt
from adjustText import adjust_text


def count_codons(fasta_file):
    codon_counts = defaultdict(int)

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                codon_counts[codon] += 1

    return codon_counts

def count_codons_per_gene(fasta_file):
    gene_codon_counts = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        codon_counts = defaultdict(int)
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                codon_counts[codon] += 1
        gene_codon_counts[record.id] = codon_counts

    return gene_codon_counts

def calculate_frequencies(codon_counts):
    total_codons = sum(codon_counts.values())
    frequencies = {codon: count / total_codons for codon, count in codon_counts.items()}
    return frequencies

def sum_frequencies(frequencies_list):
    summed_frequencies = defaultdict(float)
    for frequencies in frequencies_list:
        for codon, freq in frequencies.items():
            summed_frequencies[codon] += freq
    return summed_frequencies

def hypergeometric_test(codon, target_counts, background_counts, total_background):
    a = target_counts.get(codon, 0)
    b = sum(target_counts.values()) - a
    c = background_counts.get(codon, 0)
    d = total_background - c

    # Perform the hypergeometric test
    p_val = stats.hypergeom.sf(a-1, total_background, c, sum(target_counts.values()))

    return p_val, a, b, c, d

def calculate_log_odds_ratio(a, b, c, d):
    # Calculate odds ratio (OR)
    odds_ratio = (a / b) / (c / d) if b != 0 and d != 0 else float('inf')

    # Calculate log odds ratio (LOR)
    log_odds_ratio = math.log(odds_ratio) if odds_ratio != float('inf') else float('inf')

    return log_odds_ratio

# Paths to your FASTA files
upregulated_fasta = "Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_up.fasta"
downregulated_fasta = "Riboseq_Total_RNAseq_Davide_Like_topTable_Unt_Vs_TGFb_down.fasta"
all_genes_fasta = "../m39_transcript_region_fasta/protein_coding_CDS.fasta"

# Count codons in upregulated, downregulated, and all protein-coding genes
upregulated_codon_counts = count_codons(upregulated_fasta)
downregulated_codon_counts = count_codons(downregulated_fasta)
all_genes_codon_counts = count_codons(all_genes_fasta)

# Calculate frequencies for each gene
upregulated_codon_counts_per_gene = count_codons_per_gene(upregulated_fasta)
downregulated_codon_counts_per_gene = count_codons_per_gene(downregulated_fasta)
all_genes_codon_counts_per_gene = count_codons_per_gene(all_genes_fasta)

# Calculate frequencies for each gene
upregulated_frequencies_list = [calculate_frequencies(counts) for counts in upregulated_codon_counts_per_gene.values()]
downregulated_frequencies_list = [calculate_frequencies(counts) for counts in downregulated_codon_counts_per_gene.values()]
all_genes_frequencies_list = [calculate_frequencies(counts) for counts in all_genes_codon_counts_per_gene.values()]

# Sum frequencies across all genes
upregulated_summed_frequencies = sum_frequencies(upregulated_frequencies_list)
downregulated_summed_frequencies = sum_frequencies(downregulated_frequencies_list)
all_genes_summed_frequencies = sum_frequencies(all_genes_frequencies_list)

# Total number of codons in all protein-coding genes
total_background = sum(all_genes_codon_counts.values())

# Perform hypergeometric test and calculate log odds ratio for each codon
results = []
for codon in set(upregulated_codon_counts.keys()).union(downregulated_codon_counts.keys()):
    up_p_val, up_a, up_b, up_c, up_d = hypergeometric_test(codon, upregulated_codon_counts, all_genes_codon_counts, total_background)
    down_p_val, down_a, down_b, down_c, down_d = hypergeometric_test(codon, downregulated_codon_counts, all_genes_codon_counts, total_background)

    log_odds_ratio = calculate_log_odds_ratio(up_a, up_b, down_a, down_b)

    # Add * if p-value < 0.001
    p_val_up_str = f"{up_p_val:.3e}" + ("*" if up_p_val < 0.001 else "")
    p_val_down_str = f"{down_p_val:.3e}" + ("*" if down_p_val < 0.001 else "")

    # Stop codons
    if codon in ["TAA", "TAG", "TGA"]:
        codon = codon + "*"

    up_freq_sum = upregulated_summed_frequencies.get(codon, 0)
    down_freq_sum = downregulated_summed_frequencies.get(codon, 0)
    all_genes_freq_sum = all_genes_summed_frequencies.get(codon, 0)
    #all_genes_freq_sum = up_freq_sum + down_freq_sum
    total_freq_sum = up_freq_sum + down_freq_sum

    results.append([codon, up_a, up_b,
                    down_a, down_b,
                    p_val_up_str, p_val_down_str,
                    log_odds_ratio,
                    up_freq_sum,
                    down_freq_sum,
                    all_genes_freq_sum,
                    total_freq_sum])

# Sort the results alphabetically by codon
results.sort(key=lambda x: x[0])


# Convert results to a dictionary format
results_dict = {
    "codon": [row[0] for row in results],
    "log_odds_ratios": [row[7] for row in results],
    "total_freq_sum": [row[11] for row in results]
}



# Extract data for plotting using the dictionary keys
codons = results_dict["codon"]
log_odds_ratios = results_dict["log_odds_ratios"]
total_freq_sums = results_dict["total_freq_sum"]

# Create a dot plot
plt.figure(figsize=(10, 6))
#plt.scatter(log_odds_ratios, total_freq_sums, alpha=0.3)

# Annotate each point with the codon name
#for i, row in enumerate(results):
#    plt.annotate(row[0], (log_odds_ratios[i], total_freq_sums[i]), fontsize=8, alpha=0.7)

# Plot each point with the specified color conditions
texts = []
for i, (log_odds_ratio, total_freq_sum) in enumerate(zip(log_odds_ratios, total_freq_sums)):
    if total_freq_sum > 10 and log_odds_ratio < -0.2:
        color = 'red'
    elif total_freq_sum > 10 and log_odds_ratio > 0.2:
        color = 'blue'
    else:
        color = 'gray'
    plt.scatter(log_odds_ratio, total_freq_sum, color=color, alpha=0.6)
    # Convert results to a dictionary format
results_dict = {
    "codon": [row[0] for row in results],
    "log_odds_ratios": [row[7] for row in results],
    "total_freq_sum": [row[11] for row in results]
}



# Extract data for plotting using the dictionary keys
codons = results_dict["codon"]
log_odds_ratios = results_dict["log_odds_ratios"]
total_freq_sums = results_dict["total_freq_sum"]

# Create a dot plot
plt.figure(figsize=(12, 10))
#plt.scatter(log_odds_ratios, total_freq_sums, alpha=0.3)

# Annotate each point with the codon name
#for i, row in enumerate(results):
#    plt.annotate(row[0], (log_odds_ratios[i], total_freq_sums[i]), fontsize=8, alpha=0.7)

# Plot each point with the specified color conditions
texts = []
for i, (log_odds_ratio, total_freq_sum) in enumerate(zip(log_odds_ratios, total_freq_sums)):
    if total_freq_sum > 10 and log_odds_ratio < -0.2:
        color = 'red'
    elif total_freq_sum > 10 and log_odds_ratio > 0.2:
        color = 'blue'
    else:
        color = 'gray'
    plt.scatter(log_odds_ratio, total_freq_sum, color=color, alpha=0.4)
    texts.append(plt.text(log_odds_ratio, total_freq_sum, results[i][0], fontsize=8, alpha=0.75))

# Adjust text to avoid overlapping and add lines
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.75))


plt.xlabel('Log Odds Ratio (Log OR)')
plt.ylabel('Total Frequency Sum (Freq_Sum)')
plt.title('Riboseq Unt Vs TGFb up Frequency Sum vs Log OR')

# Add dotted lines at specified positions
plt.axhline(y=10, color='gray', linestyle='dotted')
plt.axvline(x=0.1, color='gray', linestyle='dotted')
plt.axvline(x=-0.1, color='gray', linestyle='dotted')

# Remove grid
plt.grid(False)

# Set x-axis scale from -1 to 1
plt.xlim(-1, 1)


# Save the plot as a PNG file
plt.savefig("Riboseq_Unt_Vs_TGFb_up_dot_plot_with_lines.png")

print("The dot plot with lines has been saved as Ribisqe_Unt_Vs_TGFb_up_dot_plot_with_lines.png")


# Print the results in a table format
print(tabulate(results,
               headers=["Codon", "Up_A", "Up_B",
                        "Down_A", "Down_B",
                        "Up_Pval", "Down_Pval",
                        "Log OR",
                        "Up_Freq_Sum",
                        "Dn_Freq_Sum",
                        "PC_Freq_Sum",
                        "DE_Freq_Sum"]))

# Write the results to a text file
with open("Riboseq_Unt_Vs_TGFb_up_codon_log_odds_results.txt", "w") as f:
    f.write(tabulate(results,
        headers=["Codon", "Up_A", "Up_B",
            "Down_A", "Down_B",
            "Up_Pval", "Down_Pval",
            "Log OR", "Up_Freq_Sum",
            "Dn_Freq_Sum", "PC_Freq_Sum", "DE_Freq_Sum"]))

print("Results have been written to codon_log_odds_results.txt")
