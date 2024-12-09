from Bio import SeqIO
from collections import defaultdict
import scipy.stats as stats
import math
from tabulate import tabulate
import matplotlib.pyplot as plt


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

def average_frequencies(frequencies_list):
    summed_frequencies = sum_frequencies(frequencies_list)
    num_genes = len(frequencies_list)
    avg_frequencies = {codon: freq / num_genes for codon, freq in summed_frequencies.items()}
    return avg_frequencies

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

def calculate_standard_error(a, b, c, d):
    return math.sqrt(1/a + 1/b + 1/c + 1/d)

def calculate_z_score(log_odds_ratio, se):
    return log_odds_ratio / se

def calculate_p_value(z):
    return 2 * (1 - stats.norm.cdf(abs(z)))

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

# Calculate average frequencies across all genes
upregulated_avg_frequencies = average_frequencies(upregulated_frequencies_list)
downregulated_avg_frequencies = average_frequencies(downregulated_frequencies_list)
all_genes_avg_frequencies = average_frequencies(all_genes_frequencies_list)

# Total number of codons in all protein-coding genes
total_background = sum(all_genes_codon_counts.values())

# Perform hypergeometric test and calculate log odds ratio for each codon
results = []
for codon in set(upregulated_codon_counts.keys()).union(downregulated_codon_counts.keys()):
    up_p_val, up_a, up_b, up_c, up_d = hypergeometric_test(codon, upregulated_codon_counts, all_genes_codon_counts, total_background)
    down_p_val, down_a, down_b, down_c, down_d = hypergeometric_test(codon, downregulated_codon_counts, all_genes_codon_counts, total_background)

    log_odds_ratio = calculate_log_odds_ratio(up_a, up_b, down_a, down_b)

    # Calculate standard error and p-value for log odds ratio
    se_log_odds_ratio = calculate_standard_error(up_a, up_b, down_a, down_b)
    z_score_log_odds_ratio = calculate_z_score(log_odds_ratio, se_log_odds_ratio)
    p_value_log_odds_ratio = calculate_p_value(z_score_log_odds_ratio)

    # Add * if p-value < 0.001
    p_val_up_str = f"{up_p_val:.3e}" + ("*" if up_p_val < 0.001 else "")
    p_val_down_str = f"{down_p_val:.3e}" + ("*" if down_p_val < 0.001 else "")

    # Stop codons
    if codon in ["TAA", "TAG", "TGA"]:
        codon = codon + "*"

    up_freq_sum = upregulated_summed_frequencies.get(codon, 0)
    down_freq_sum = downregulated_summed_frequencies.get(codon, 0)
    all_genes_freq_sum = up_freq_sum + down_freq_sum
    total_freq_sum = up_freq_sum + down_freq_sum

    up_avg_freq = upregulated_avg_frequencies.get(codon, 0)
    down_avg_freq = downregulated_avg_frequencies.get(codon, 0)
    all_genes_avg_freq = all_genes_avg_frequencies.get(codon, 0)

    results.append([codon, up_a, up_b,
                    down_a, down_b,
                    p_val_up_str, p_val_down_str,
                    log_odds_ratio,
                    up_freq_sum,
                    down_freq_sum,
                    all_genes_freq_sum,
                    total_freq_sum,
                    p_value_log_odds_ratio,
                    up_avg_freq,
                    down_avg_freq,
                    all_genes_avg_freq])

# Sort the results alphabetically by codon
results.sort(key=lambda x: x[0])



#####################################################
#####################################################
# Plotting
#####################################################
#####################################################

################    For TE down   ###################
down_avg_freqs = [row[14] for row in results]
up_avg_freqs = [row[13] for row in results]
all_genes_avg_freqs = [row[15] for row in results]

# Calculate log ratios
log_down_up_ratios = [math.log(down / up) if up != 0 else float('inf') for down, up in zip(down_avg_freqs, up_avg_freqs)]
log_down_all_ratios = [math.log(down / all_genes) if all_genes != 0 else float('inf') for down, all_genes in zip(down_avg_freqs, all_genes_avg_freqs)]

# Create a dot plot
plt.figure(figsize=(10, 6))
plt.scatter(log_down_up_ratios, log_down_all_ratios, alpha=0.6)

# Annotate each point with the codon name
for i, row in enumerate(results):
    plt.annotate(row[0], (log_down_up_ratios[i], log_down_all_ratios[i]), fontsize=8, alpha=0.7)


# Add dotted lines at specified positions
plt.axhline(y=0, color='gray', linestyle='dotted')
plt.axvline(x=0, color='gray', linestyle='dotted')


plt.xlabel('Log(Down_Avg_Freq / Up_Avg_Freq)')
plt.ylabel('Log(Down_Avg_Freq / All_Genes_Avg_Freq)')
plt.title('Dot Plot of Log Ratios Riboseq_Total_RNAseq Unt_Vs_TGFb_up')

# Save the plot as a PNG file
plt.savefig("Riboseq_Total_RNAseq_Unt_Vs_TGFb_up_dot_plot_down_log_ratios.png")

print("The dot plot of log ratios has been saved as Riboseq_Total_RNAseq_dot_plot_down_log_ratios.png")


#################    For TE up   #####################
down_avg_freqs = [row[14] for row in results]
up_avg_freqs = [row[13] for row in results]
all_genes_avg_freqs = [row[15] for row in results]

# Calculate log ratios
log_up_down_ratios = [math.log(up / down) if up != 0 else float('inf') for up, down in zip(up_avg_freqs, down_avg_freqs,)]
log_up_all_ratios = [math.log(up / all_genes) if all_genes != 0 else float('inf') for up, all_genes in zip(up_avg_freqs, all_genes_avg_freqs)]

# Create a dot plot
plt.figure(figsize=(10, 6))
plt.scatter(log_up_down_ratios, log_up_all_ratios, alpha=0.6)

# Annotate each point with the codon name
for i, row in enumerate(results):
    plt.annotate(row[0], (log_up_down_ratios[i], log_up_all_ratios[i]), fontsize=8, alpha=0.7)


# Add dotted lines at specified positions
plt.axhline(y=0, color='gray', linestyle='dotted')
plt.axvline(x=0, color='gray', linestyle='dotted')


plt.xlabel('Log(Up_Avg_Freq / Down_Avg_Freq)')
plt.ylabel('Log(Up_Avg_Freq / All_Genes_Avg_Freq)')
plt.title('Dot Plot of Log Ratios Riboseq_Total_RNAseq Unt_Vs_TGFb_up')

# Save the plot as a PNG file
plt.savefig("Riboseq_Total_RNAseq_Unt_Vs_TGFb_up_dot_plot_up_log_ratios.png")

print("The dot plot of log ratios has been saved as Riboseq_Total_RNAseq_dot_plot_up_log_ratios.png")

##########################################################
##########################################################
# Print the results in a table format
##########################################################
##########################################################

print(tabulate(results,
               headers=["Codon", "Up_A", "Up_B",
                        "Down_A", "Down_B",
                        "Up_P-value", "Down_P-value",
                        "Log OR",
                        "Up_Freq_Sum",
                        "Down_Freq_Sum",
                        "All_Genes_Freq_Sum",
                        "Total_Freq_Sum",
                        "Log OR P-value",
                        "Up_Avg_Freq",
                        "Down_Avg_Freq",
                        "All_Genes_Avg_Freq"]))

# Write the results to a text file
with open("codon_log_odds_results.txt", "w") as f:
    f.write(tabulate(results,
                      headers=["Codon", "Up_A", "Up_B",
                               "Down_A", "Down_B",
                               "Up_P-value", "Down_P-value",
                               "Log OR",
                               "Up_Freq_Sum",
                               "Down_Freq_Sum",
                               "All_Genes_Freq_Sum",
                               "Total_Freq_Sum",
                               "Log OR P-value",
                               "Up_Avg_Freq",
                               "Down_Avg_Freq",
                               "All_Genes_Avg_Freq"]))

print("Results have been written to Riboseq_Total_RNAseq_codon_Unt_Vs_TGFb_log_odds_results.txt")
