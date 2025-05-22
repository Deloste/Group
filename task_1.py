"""
pseudocode：
Function most_frequent_nucleotide with input mRNA_sequence:
    Initialize stop_codons as a set containing 'UAA', 'UAG', 'UGA'
    Create an empty dictionary called codon_freq to count occurrences
    Initialize index to 0
    
    While index is less than length of mRNA_sequence:
        Extract codon from mRNA_sequence starting at index, 3 characters long
        If codon length is less than 3 OR codon is in stop_codons:
            Exit loop
        Else:
            Increment count for this codon in codon_freq
            Increase index by 3
    
    If codon_freq is empty:
        Return structure with null trinucleotides and 0 count
    Else:
        Find maximum value in codon_freq counts
        Create list of codons with this maximum count
        Convert these codons to DNA format by replacing 'U' with 'T'
        Return structure with DNA trinucleotides list and max count

Test Case:
    test_seq = "AUGAUGUUCUUCGUAA"
    result = call most_frequent_nucleotide with test_seq
    If result contains trinucleotides:
        Print "Most frequent DNA trinucleotide(s) and frequency:"
        For each trinucleotide in result's list:
            Print trinucleotide and its count
    Else:
        Print "No valid codons found before stop codon."

Function read_sequence_txt with input filepath:
    Open file at filepath for reading
    Read file content, remove whitespace and newline characters
    Return cleaned sequence as single continuous string
"""
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

def read_sequence_txt(filepath):
    with open(filepath, 'r') as file:
        sequence = file.read().strip().replace('\n', '')
    return sequence
