"""
Auxilliary routines for DNA/RNA/AA sequences processing.
"""

# from Bio.Seq import Seq


DNA_LETTERS = list('ACGT')
RNA_LETTERS = list('ACGU')
AA_LETTERS = list('ACDEFGHIKLMNPQRSTVWY')


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
