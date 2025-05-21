from collections import Counter
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

CODON_TABLE = {
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
    'UUA': 'Leucine', 'UUG': 'Leucine',
    'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
    'AUG': 'Methionine',
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
    'CAU': 'Histidine', 'CAC': 'Histidine',
    'CAA': 'Glutamine', 'CAG': 'Glutamine',
    'AAU': 'Asparagine', 'AAC': 'Asparagine',
    'AAA': 'Lysine', 'AAG': 'Lysine',
    'GAU': 'Aspartic acid', 'GAC': 'Aspartic acid',
    'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',
    'UGU': 'Cysteine', 'UGC': 'Cysteine',
    'UGG': 'Tryptophan',
    'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AGU': 'Serine', 'AGC': 'Serine',
    'AGA': 'Arginine', 'AGG': 'Arginine',
    'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine',
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

def plot_amino_acid_frequencies(mRNA_sequence):
    stop_codons = {'UAA', 'UAG', 'UGA'}
    freq = Counter()
    for i in range(0, len(mRNA_sequence), 3):
        codon = mRNA_sequence[i:i+3]
        if len(codon) < 3 or codon in stop_codons:
            break
        aa = CODON_TABLE.get(codon)
        if aa and aa != 'Stop':
            freq[aa] += 1

    if not freq:
        print("No amino acids to show — check your sequence?")
        return
    # Plotting
    sorted_items = sorted(freq.items(), key=lambda x: x[1], reverse=True) # Sort the amino acids by frequency
    labels, values = zip(*sorted_items) # Unzip the labels and values
    norm = Normalize(vmin=0, vmax=max(values) or 1) # Normalize the values for color mapping
    cmap = cm.get_cmap('GnBu') # Choose a colormap
    colors = [cmap(norm(v)) for v in values] 
    plt.figure(figsize=(10, 6))
    bars = plt.barh(labels, values, color=colors, edgecolor='black') # Horizontal bar chart
    for bar in bars:
        width = bar.get_width() 
        plt.text(width + 0.5, bar.get_y() + bar.get_height()/2, f"{int(width)}", va='center', fontsize=10)
    plt.xlabel("Count", fontsize=12)
    plt.title("Amino Acid Frequency Distribution", fontsize=16, fontweight='bold')
    plt.gca().invert_yaxis() # Invert y-axis to have the most frequent amino acid on top
    plt.tight_layout() # Adjust layout to fit labels
    plt.xlim(0, max(values) + 2) # Add some space to the right of the bars
    plt.show()


test_seq ="AUGGCUACCUUUGGGAUCAAGGCGUUGGAGCUGCUUCCUGACUGGAAUUGGACCGGAACCCUUGGUAC"
print("Plotting amino acid frequencies for test_seq...")
plot_amino_acid_frequencies(test_seq)
