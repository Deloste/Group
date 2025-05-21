from collections import Counter

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

def most_frequent_amino_acid(mRNA_sequence):
    stop_codons = {'UAA', 'UAG', 'UGA'}
    codon_freq = Counter()
    for i in range(0, len(mRNA_sequence), 3):
        codon = mRNA_sequence[i:i+3]
        if len(codon) < 3 or codon in stop_codons:
            break
        codon_freq[codon] += 1
    if not codon_freq:
        return None, 0
    # Find the most frequent codon(s)
    max_count = max(codon_freq.values())
    most_freq_codons = [codon for codon, count in codon_freq.items() if count == max_count]

    # Translate the most frequent codon(s) to amino acids
    translated_amino_acids = [CODON_TABLE[codon] for codon in most_freq_codons]
    # Remove duplicates
    amino_acids = list(set(translated_amino_acids))
    return {
        "most_frequent_amino_acids": amino_acids,
        "count": max_count
    }


test_seq = "AUGGUUGUUCUUCUUUAG"

result = most_frequent_amino_acid(test_seq)
if result["most_frequent_amino_acids"]:
    print("Most frequent amino acid(s) and frequency:")
    for amino_acid in result["most_frequent_amino_acids"]:
        print(f"{amino_acid}: {result['count']}")
else:
    print("No valid codons found before stop codon.")