"""
Auxilliary routines for DNA/RNA/AA sequences processing.
"""

# from Bio.Seq import Seq


DNA_LETTERS = list('ACGT')
RNA_LETTERS = list('ACGU')
AA_LETTERS = list('ACDEFGHIKLMNPQRSTVWY')


# translatable RNA codons
rna_codon_table = {
	'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 
	'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 
	'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 
	'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 
	'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 
	'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 
	'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 
	'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'
}


rna_stop_codons = ['UAA', 'UAG', 'UGA']


# amino acids composition by atoms
aa_composition_table = {
	'A': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 0},
	'C': {'C': 3, 'H': 7, 'N': 1, 'O': 2, 'S': 1},
	'D': {'C': 4, 'H': 7, 'N': 1, 'O': 4, 'S': 0},
	'E': {'C': 5, 'H': 9, 'N': 1, 'O': 4, 'S': 0},
	'F': {'C': 9, 'H': 11, 'N': 1, 'O': 2, 'S': 0},
	'G': {'C': 2, 'H': 5, 'N': 1, 'O': 2, 'S': 0},
	'H': {'C': 6, 'H': 9, 'N': 3, 'O': 2, 'S': 0},
	'I': {'C': 6, 'H': 13, 'N': 1, 'O': 2, 'S': 0},
	'K': {'C': 6, 'H': 14, 'N': 2, 'O': 2, 'S': 0},
	'L': {'C': 6, 'H': 13, 'N': 1, 'O': 2, 'S': 0},
	'M': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 1},
	'N': {'C': 4, 'H': 8, 'N': 2, 'O': 3, 'S': 0},
	'P': {'C': 5, 'H': 9, 'N': 1, 'O': 2, 'S': 0},
	'Q': {'C': 5, 'H': 10, 'N': 2, 'O': 3, 'S': 0},
	'R': {'C': 6, 'H': 14, 'N': 4, 'O': 2, 'S': 0},
	'S': {'C': 3, 'H': 7, 'N': 1, 'O': 3, 'S': 0},
	'T': {'C': 4, 'H': 9, 'N': 1, 'O': 3, 'S': 0},
	'V': {'C': 5, 'H': 11, 'N': 1, 'O': 2, 'S': 0},
	'W': {'C': 11, 'H': 12, 'N': 2, 'O': 2, 'S': 0},
	'Y': {'C': 9, 'H': 11, 'N': 1, 'O': 3, 'S': 0}
}


# Kyte & Doolittle index of hydrophobicity
aa_kd_table = {
	'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
	'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
	'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
	'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}


aa_codes_table = {
	'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
	'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
	'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
	'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'
}


aa_half_life = {
    "A": ["4.4 hour", ">20 hour", ">10 hour"], 
    "C": ["1.2 hour", ">20 hour", ">10 hour"], 
    "E": ["1 hour", "30 min", ">10 hour"], 
    "D": ["1.1 hour", "3 min", ">10 hour"], 
    "G": ["30 hour", ">20 hour", ">10 hour"], 
    "F": ["1.1 hour", "3 min", "2 min"], 
    "I": ["20 hour", "30 min", ">10 hour"], 
    "H": ["3.5 hour", "10 min", ">10 hour"], 
    "K": ["1.3 hour", "3 min", "2 min"], 
    "M": ["30 hour", ">20 hour", ">10 hour"], 
    "L": ["5.5 hour", "3 min", "2 min"], 
    "N": ["1.4 hour", "3 min", ">10 hour"], 
    "Q": ["0.8 hour", "10 min", ">10 hour"], 
    "P": [">20 hour", ">20 hour", "?"], 
    "S": ["1.9 hour", ">20 hour", ">10 hour"], 
    "R": ["1 hour", "2 min", "2 min"], 
    "T": ["7.2 hour", ">20 hour", ">10 hour"], 
    "W": ["2.8 hour", "3 min", "2 min"], 
    "V": ["100 hour", ">20 hour", ">10 hour"], 
    "Y": ["2.8 hour", "10 min", "2 min"]
}


atom_codes_table = {'C': 'Carbon', 'H': 'Hydrogen', 'N': 'Nitrogen', 'O': 'Oxygen', 'S': 'Sulfur'}


def rna_to_aa(seq):
	assert len(seq) > 0, 'Empty sequence.'
	assert len(seq) % 3 == 0, 'Sequence length is not product of 3.'
	assert all(map(lambda x: x in RNA_LETTERS, seq)), 'Not an RNA sequence.'
	assert all(map(
		lambda codon: codon not in rna_stop_codons,
		map(lambda i: seq[i:i+3], xrange(0, len(seq), 3))
	)), 'Sequence should not contain stop codons.'

	return ''.join(map(
		lambda codon: rna_codon_table[codon], 
		map(lambda i: seq[i:i+3], xrange(0, len(seq), 3))
	))


def dna_to_rna(seq):
	assert all(map(lambda x: x in DNA_LETTERS, seq)), 'Not an DNA sequence.'
	
	return seq.replace('T', 'U')


def dna_to_aa(seq):
	assert all(map(lambda x: x in DNA_LETTERS, seq)), 'Not an DNA sequence.'
	
	return rna_to_aa(dna_to_rna(seq))
