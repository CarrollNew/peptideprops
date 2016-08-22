"""
The module tests implemented calculation versus Biopython results.
"""

import unittest
import urllib2

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sequence_property import InstabilityIndex, IsoelectricPoint, MolecularWeight
from seq_utils import rna_codon_table, rna_to_aa, RNA_LETTERS


class BioTests(unittest.TestCase):
	def test_random(self):
		"""
		Compare properties on random sequences.
		"""
		test_status = True
		num_tries = 20
		for peptide_len in xrange(10, 101):
			for _ in xrange(num_tries):
				seq = BioTests._gen_random_mrna(peptide_len)
				test_status &= BioTests._check_props(seq)

		self.assertTrue(test_status)

	def test_uniprot(self):
		"""
		Compare properties on UniProt sequences.
		"""
		# TODO: add some random uniprot entries calculation & comparison here
		pass

	@staticmethod
	def _check_props(seq):
		return (
			BioTests._check_ii_prop(seq) 
			and BioTests._check_pi_prop(seq) 
			and BioTests._check_mw_prop(seq)
		)

	@staticmethod
	def _check_ii_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).instability_index(),
			InstabilityIndex(seq).calculate_prop()
		)

	@staticmethod
	def _check_pi_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).isoelectric_point(),
			IsoelectricPoint(seq).calculate_prop()
		)

	@staticmethod
	def _check_mw_prop(seq):
		return BioTests._is_value_valid(
			ProteinAnalysis(rna_to_aa(seq)).molecular_weight(),
			MolecularWeight(seq).calculate_prop()
		)

	@staticmethod
	def _gen_random_mrna(peptide_len=30):
		return ''.join(np.random.choice(rna_codon_table.keys(), size=peptide_len, replace=True))

	@staticmethod
	def _is_value_valid(ref_value, value):
		if isinstance(ref_value, int) and isinstance(value, int):
			return ref_value == value
		elif isinstance(ref_value, float) or isinstance(value, float):
			return np.isclose(ref_value, value, rtol=0.05)

	@staticmethod
	def _read_random_uniprot_entry():
		entry_url_template = 'http://www.uniprot.org/uniprot/'
		random_entry_url = entry_url_template + '&'.join([
			'?query=reviewed:yes+AND+organism:9606', 'random=yes', 'format=txt'])
		data = urllib2.urlopen(random_entry_url).read()


if __name__ == '__main__':
	unittest.main()
