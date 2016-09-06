import os
import subprocess
import tempfile
import string
import random
import unittest
from json import dumps, loads

from seq_utils import RNA_LETTERS


class RNAPlotter(object):
	possible_rnaplot_paths = [
		'C://Program Files (x86)//ViennaRNA Package//RNAplot.exe',
		'/usr/local/bin/ViennaRNA/RNAplot',
		'/usr/bin/ViennaRNA/RNAplot',
		'/usr/bin/RNAplot',
		'/usr/local/bin/RNAplot'
	]

	def __init__(self, rna_seq, rna_ss, params=None):
		RNAPlotter._validate(rna_seq, rna_ss)
		self.rna_seq = rna_seq
		self.rna_ss = rna_ss

	def create_plot(self, fmt='svg'):
		seq_name = RNAPlotter._create_random_seq_name()
		in_file_name = tempfile.gettempdir() + os.sep + 'rnaplot_input.txt'
		out_file_name = '{}_ss.{}'.format(seq_name, fmt)
		self._write_input(in_file_name, seq_name)
		command = 'cat {} | {} -o {}'.format(in_file_name, RNAPlotter._get_rnaplot_app(), fmt)
		
		try:
			subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
			assert os.path.exists(out_file_name), 'Could not find result of RNAPlot.'
		except subprocess.CalledProcessError as e:
			print e.message
			out_file_name = ''

		os.remove(in_file_name)

		return out_file_name

	def _write_input(self, fn, name):
		with open(fn, 'w') as f:
			f.write('>{}\n'.format(name))
			f.write(self.rna_seq + '\n')
			f.write(self.rna_ss + '\n')

	@staticmethod
	def from_json(json_input):
		rna_seq = json_input['rna_seq']
		rna_ss = json_input['rna_ss']
		fmt = json_input['format'] if 'format' in json_input else 'svg'
		instance = RNAPlotter(rna_seq, rna_ss)
		return instance.create_plot(fmt)

	@staticmethod
	def _create_random_seq_name(seq_len=10):
		assert seq_len > 0, 'Can not create empty sequence name.'

		return ''.join(random.choice(string.ascii_letters + string.digits) for _ in xrange(seq_len))

	@staticmethod
	def _validate(rna_seq, rna_ss):
		assert all(map(lambda x: x in RNA_LETTERS, rna_seq)), \
			'Not an RNA sequence.'
		assert all(map(lambda x: x in ['(', ')', '.'], rna_ss)), \
			'Secondary structure should be in dot format.'
		assert len(rna_seq) == len(rna_ss), \
			'Sequence and secondary structure should have same length.'

	@staticmethod
	def _get_rnaplot_app():
		for app in RNAPlotter.possible_rnaplot_paths:
			if os.path.exists(app):
				return app
		raise ValueError('Could not find suitable RNAPlot app directory.')


class RNAPlotterTests(unittest.TestCase):
	def test_no_ss(self):
		num_tries = 100
		for _ in xrange(num_tries):
			for seq_len in xrange(10, 100):
				seq = ''.join(random.choice(RNA_LETTERS) for _ in xrange(seq_len))
				ss = ''.join(['.'] * seq_len)
				instance = RNAPlotter(seq, ss)
				fn_name = instance.create_plot()
				self.assertTrue(os.path.exists(fn_name))
				os.remove(fn_name)

	def test_formats(self):
		num_tries, seq_len = 100, 10
		for _ in xrange(num_tries):
			for fmt in ['ps', 'gml', 'svg']:
				seq = ''.join(random.choice(RNA_LETTERS) for _ in xrange(seq_len))
				ss = ''.join(['.'] * seq_len)
				instance = RNAPlotter(seq, ss)
				fn_name = instance.create_plot(fmt)
				self.assertTrue(os.path.exists(fn_name))
				os.remove(fn_name)

	def test_with_ss(self):
		def _create_valid_ss(seq):
			assert len(seq) > 5, 'I would not fold 5 or less b.p. RNAs.'

			split_point = len(seq) // 2
			num_brackets = random.randint(1, split_point)
			beg_inds = random.sample(xrange(split_point), num_brackets)
			end_inds = random.sample(xrange(split_point, 2 * split_point), num_brackets)
			return ''.join(['(' if i in beg_inds else ')' if i in end_inds else '.'
							for i in xrange(len(seq))])

		num_tries = 100
		for _ in xrange(num_tries):
			for seq_len in xrange(50, 100):
				seq = ''.join(random.choice(RNA_LETTERS) for _ in xrange(seq_len))
				ss = _create_valid_ss(seq)
				instance = RNAPlotter(seq, ss)
				fn_name = instance.create_plot()
				self.assertTrue(os.path.exists(fn_name))
				os.remove(fn_name)


if __name__ == '__main__':
	unittest.main()
