
from collections import Counter

def most_frequent_nucleotide(mRNA_sequence):
    stop_codons = {'UAA', 'UAG', 'UGA'}
    codon_freq = Counter()
    for i in range(0, len(mRNA_sequence), 3):
        codon = mRNA_sequence[i:i+3]
        if len(codon) < 3 or codon in stop_codons:
            break
        codon_freq[codon] += 1
    if not codon_freq:
        return None, 0
    max_count = max(codon_freq.values())
    most_freq_codons = [codon for codon, count in codon_freq.items() if count == max_count]
    DNA_trinucleotides = [codon.replace('U', 'T') for codon in most_freq_codons]

    return { "most_frequent_trinucleotides": DNA_trinucleotides, 
            "count": max_count}

test_seq = "AUGAUGUUCUUCGUAA"

result = most_frequent_nucleotide(test_seq)
if result["most_frequent_trinucleotides"]:
    print("Most frequent DNA trinucleotide(s) and frequency:")
    for trinucleotide in result["most_frequent_trinucleotides"]:
        print(f"{trinucleotide}: {result['count']}")
else:
    print("No valid codons found before stop codon.")
