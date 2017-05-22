"""
The module tests implemented calculation versus Biopython/ExPASy results.
"""

import unittest
import os
from json import loads
from random import choice, sample

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sequence_property import (
	GRAVY,
	InstabilityIndex,
	IsoelectricPoint,
	MolecularWeight,
	SequenceComposition,

	ExtinctionCoefficient,
	AbsorptionCoefficient,
	AliphaticIndex
)
from seq_utils import rna_codon_table, rna_to_aa, RNA_LETTERS, DNA_LETTERS, AA_LETTERS
from main import process


class BioTests(unittest.TestCase):
	def test_random(self):
		"""
		Compare properties on random sequences (using Biopython).
		"""
		test_status = True
		num_tries = 20
		for peptide_len in xrange(10, 101):
			for _ in xrange(num_tries):
				seq = BioTests._gen_random_mrna(peptide_len)
				test_status &= BioTests._check_props(seq)

		self.assertTrue(test_status)

	def test_expasy(self):
		tests_path = 'expasy_tests'

		if not os.path.exists(tests_path):
			raise IOError('No ExPASy tests found.')
		for fn in filter(lambda fn: fn.split('.')[-1] == 'json', os.listdir(tests_path)):
			with open(tests_path + os.sep + fn) as f:
				item = loads(f.read())
				expasy_result = item['results']
				aa_seq = item['sequence']
				self.assertTrue(BioTests._check_expasy_result(expasy_result, aa_seq))

	@staticmethod
	def _check_expasy_result(result, seq):
		extract_value = lambda name: filter(lambda x: x['property'] == name, result)[0]['value']

		composition = SequenceComposition(seq).calculate_prop()
		ref_atom_count = extract_value('atom_count')
		atom_count = composition.atom_count
		return (
			all(BioTests._is_value_valid(ref_atom_count[let], atom_count[let]) for let in ref_atom_count) and
			BioTests._is_value_valid(extract_value('atom_num'), composition.atom_num) and
			BioTests._is_value_valid(extract_value('neg_charged_res'), composition.aa_count['D'] + composition.aa_count['E']) and
			BioTests._is_value_valid(extract_value('pos_charged_res'), composition.aa_count['R'] + composition.aa_count['K']) and
			BioTests._is_value_valid(extract_value('extinction_coefficient'), ExtinctionCoefficient(seq).calculate_prop()) and
			BioTests._is_value_valid(extract_value('absorption_coefficient'), AbsorptionCoefficient(seq).calculate_prop())
		)


	@staticmethod
	def _check_props(seq):
		return (
			BioTests._check_ii_prop(seq) 
			and BioTests._check_pi_prop(seq) 
			and BioTests._check_mw_prop(seq)
			and BioTests._check_composition_prop(seq)
			and BioTests._check_gravy_prop(seq)
		)

	@staticmethod
	def _check_ii_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).instability_index(),
			InstabilityIndex(rna_to_aa(seq)).calculate_prop()
		)

	@staticmethod
	def _check_pi_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).isoelectric_point(),
			IsoelectricPoint(rna_to_aa(seq)).calculate_prop()
		)

	@staticmethod
	def _check_mw_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).molecular_weight(),
			MolecularWeight(rna_to_aa(seq)).calculate_prop()
		)

	@staticmethod
	def _check_composition_prop(seq):
		ref_count = ProteinAnalysis(rna_to_aa(seq)).count_amino_acids()
		aa_count = SequenceComposition(rna_to_aa(seq)).calculate_prop().aa_count
		return all(
			BioTests._is_value_valid(ref_count[let], aa_count[let])
			for let in ref_count
		)

	@staticmethod
	def _check_gravy_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).gravy(),
			GRAVY(rna_to_aa(seq)).calculate_prop()
		)

	@staticmethod
	def _gen_random_mrna(peptide_len=30):
		return ''.join(map(lambda _: choice(rna_codon_table.keys()), xrange(peptide_len)))

	@staticmethod
	def _is_value_valid(ref_value, value):
		if isinstance(ref_value, int) and isinstance(value, int):
			return ref_value == value
		elif isinstance(ref_value, float) or isinstance(value, float):
			rtol = 0.01
			return abs(ref_value - value) <= rtol * abs(ref_value)


class Tests(unittest.TestCase):
	def test_correct(self):
		for sz in xrange(10, 101):
			self.assertTrue('error' not in process({'sequence': Tests._get_random_seq(AA_LETTERS), 'type': 'Protein'}))
			self.assertTrue('error' not in process({'sequence': Tests._get_random_seq(RNA_LETTERS), 'type': 'mRNA'}))
			self.assertTrue('error' not in process({'sequence': Tests._get_random_seq(DNA_LETTERS), 'type': 'mRNA'}))

	def test_incorrect(self):
		for sz in xrange(10, 101):
			seq = Tests._get_random_seq(RNA_LETTERS + DNA_LETTERS)
			res = process({'sequence': seq, 'type': 'mRNA'})
			if 'T' in seq and 'U' in seq:
				self.assertTrue('error' in res)

	@staticmethod
	def _get_random_seq(alphabet, sz=30):
		return ''.join(map(lambda _: choice(alphabet), xrange(sz)))


if __name__ == '__main__':
	unittest.main()
