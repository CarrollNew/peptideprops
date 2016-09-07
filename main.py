import os
import sys

import urllib2
from json import dumps, loads

from peptide_props_factory import PeptidePropsFactory
from seq_utils import rna_to_aa


def _read_params(fn):
	with open(fn) as f:
		return loads(f.read())


def write_result(out_fn, result):
	with open(out_fn, 'w') as f:
		f.write(dumps(result, indent=4))


def process(params):
	if 'sequence' not in params:
		return {'error': 'Invalid parameters', 'message': 'No mRNA sequence provided.'}
	if 'properties' not in params:
		return {'error': 'Invalid parameters', 'message': 'No settings provided.'}

	rna_seq = params['sequence']
	try:
		aa_seq = rna_to_aa(rna_seq)
	except StandardError as e:
		return {'error': 'Invalid RNA sequence.', 'message': e.message}

	results = []
	for prop in params['properties']:
		try:
			value = PeptidePropsFactory.from_prop_code(prop['name'])(rna_seq).calculate_prop()
		except KeyError as e:
			value = PeptidePropsFactory.from_long_name(prop['name'])(rna_seq).calculate_prop()
		except StandardError as e:
			print e.message
			print 'Could not calculate property: {}'.format(prop['name'])
			value = None
		results.append({'property': prop['name'], 'value': value})
	result_dict = {'sequence': aa_seq, 'results': results}
	return result_dict


def main(args):
	in_fn = args[0]
	out_fn = args[1]
	result = process(_read_params(in_fn))
	write_result(out_fn, result)


if __name__ == '__main__':
	main(sys.argv[1:])
