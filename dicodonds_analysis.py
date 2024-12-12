from Bio import SeqIO
from collections import defaultdict
import scipy.stats as stats
import math
from tabulate import tabulate
import matplotlib.pyplot as plt


# Define stop codons
stop_codons = {'TAA', 'TAG', 'TGA'}

def count_dicodons(fasta_file):
    dicodon_counts = defaultdict(int)
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        for i in range(0, len(sequence) - 5, 3):
            dicodon = sequence[i:i+6]
            if len(dicodon) == 6 and dicodon[:3] not in stop_codons:
                dicodon_counts[dicodon] += 1
    return dicodon_counts

def count_dicodons_per_gene(fasta_file):
    gene_dicodon_counts = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        dicodon_counts = defaultdict(int)
        for i in range(0, len(sequence) - 5, 3):
            dicodon = sequence[i:i+6]
            if len(dicodon) == 6 and dicodon[:3] not in stop_codons:
                dicodon_counts[dicodon] += 1
        gene_dicodon_counts[record.id] = dicodon_counts
    return gene_dicodon_counts

def calculate_frequencies(dicodon_counts):
    total_dicodons = sum(dicodon_counts.values())
    frequencies = {dicodon: count / total_dicodons for dicodon, count in dicodon_counts.items()}
    return frequencies

def sum_frequencies(frequencies_list):
    summed_frequencies = defaultdict(float)
    for frequencies in frequencies_list:
        for dicodon, freq in frequencies.items():
            summed_frequencies[dicodon] += freq
    return summed_frequencies

def average_frequencies(frequencies_list):
    summed_frequencies = sum_frequencies(frequencies_list)
    num_genes = len(frequencies_list)
    avg_frequencies = {dicodon: freq / num_genes for dicodon, freq in summed_frequencies.items()}
    return avg_frequencies

def hypergeometric_test(dicodon, target_counts, background_counts, total_background):
    a = target_counts.get(dicodon, 0)
    b = sum(target_counts.values()) - a
    c = background_counts.get(dicodon, 0)
    d = total_background - c
    p_val = stats.hypergeom.sf(a-1, total_background, c, sum(target_counts.values()))
    return p_val, a, b, c, d

def calculate_log_odds_ratio(a, b, c, d):
    if b == 0 or d == 0 or (a / b) == 0 or (c / d) == 0:
        return float('inf')
    odds_ratio = (a / b) / (c / d)
    log_odds_ratio = math.log(odds_ratio)
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

# Count dicodons in upregulated, downregulated, and all protein-coding genes
upregulated_dicodon_counts = count_dicodons(upregulated_fasta)
downregulated_dicodon_counts = count_dicodons(downregulated_fasta)
all_genes_dicodon_counts = count_dicodons(all_genes_fasta)

# Calculate frequencies for each gene
upregulated_dicodon_counts_per_gene = count_dicodons_per_gene(upregulated_fasta)
downregulated_dicodon_counts_per_gene = count_dicodons_per_gene(downregulated_fasta)
all_genes_dicodon_counts_per_gene = count_dicodons_per_gene(all_genes_fasta)

# Calculate frequencies for each gene
upregulated_frequencies_list = [calculate_frequencies(counts) for counts in upregulated_dicodon_counts_per_gene.values()]
downregulated_frequencies_list = [calculate_frequencies(counts) for counts in downregulated_dicodon_counts_per_gene.values()]
all_genes_frequencies_list = [calculate_frequencies(counts) for counts in all_genes_dicodon_counts_per_gene.values()]

# Sum frequencies across all genes
upregulated_summed_frequencies = sum_frequencies(upregulated_frequencies_list)
downregulated_summed_frequencies = sum_frequencies(downregulated_frequencies_list)
all_genes_summed_frequencies = sum_frequencies(all_genes_frequencies_list)

# Calculate average frequencies across all genes
upregulated_avg_frequencies = average_frequencies(upregulated_frequencies_list)
downregulated_avg_frequencies = average_frequencies(downregulated_frequencies_list)
all_genes_avg_frequencies = average_frequencies(all_genes_frequencies_list)

# Total number of dicodons in all protein-coding genes
total_background = sum(all_genes_dicodon_counts.values())

# Perform hypergeometric test and calculate log odds ratio for each dicodon
results = []
for dicodon in set(upregulated_dicodon_counts.keys()).union(downregulated_dicodon_counts.keys()):
    up_p_val, up_a, up_b, up_c, up_d = hypergeometric_test(dicodon, upregulated_dicodon_counts, all_genes_dicodon_counts, total_background)
    down_p_val, down_a, down_b, down_c, down_d = hypergeometric_test(dicodon, downregulated_dicodon_counts, all_genes_dicodon_counts, total_background)
    
    # Ensure no division by zero occurs and avoid math domain error
    if up_b == 0 or down_b == 0 or (up_a / up_b) == 0 or (down_a / down_b) == 0:
        log_odds_ratio = float('inf')
        se_log_odds_ratio = float('inf')
        z_score_log_odds_ratio = float('inf')
        p_value_log_odds_ratio = float('inf')
    else:
        log_odds_ratio = calculate_log_odds_ratio(up_a, up_b, down_a, down_b)
        se_log_odds_ratio = calculate_standard_error(up_a, up_b, down_a, down_b)
        z_score_log_odds_ratio = calculate_z_score(log_odds_ratio, se_log_odds_ratio)
        p_value_log_odds_ratio = calculate_p_value(z_score_log_odds_ratio)

    p_val_up_str = f"{up_p_val:.3e}" + ("*" if up_p_val < 0.001 else "")
    p_val_down_str = f"{down_p_val:.3e}" + ("*" if down_p_val < 0.001 else "")
    
    if dicodon in ["TAA", "TAG", "TGA"]:
        dicodon += "*"
    
    up_freq_sum = upregulated_summed_frequencies.get(dicodon, 0)
    down_freq_sum = downregulated_summed_frequencies.get(dicodon, 0)
    all_genes_freq_sum = up_freq_sum + down_freq_sum
    total_freq_sum = up_freq_sum + down_freq_sum
    up_avg_freq = upregulated_avg_frequencies.get(dicodon, 0)
    down_avg_freq = downregulated_avg_frequencies.get(dicodon, 0)
    all_genes_avg_freq = all_genes_avg_frequencies.get(dicodon, 0)
    results.append([dicodon, up_a, up_b, down_a, down_b, p_val_up_str, p_val_down_str, log_odds_ratio, up_freq_sum, down_freq_sum, all_genes_freq_sum, total_freq_sum, p_value_log_odds_ratio, up_avg_freq, down_avg_freq, all_genes_avg_freq])

# Print results in a table
print(tabulate(results, headers=["Dicodon", "Up_A", "Up_B", "Down_A", "Down_B", "Up_P_Val", "Down_P_Val", "Log_Odds_Ratio", "Up_Freq_Sum", "Down_Freq_Sum", "All_Genes_Freq_Sum", "Total_Freq_Sum", "P_Value_Log_Odds_Ratio", "Up_Avg_Freq", "Down_Avg_Freq", "All_Genes_Avg_Freq"]))
