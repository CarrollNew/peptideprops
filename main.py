import os
import sys

import urllib2
from json import dumps, loads

from peptide_props_factory import PeptidePropsFactory
from seq_utils import rna_to_aa


def _read_params(fn):
	with open(fn) as f:
		return loads(f.read())


def _write_result(out_fn, rna_seq, results):
	result_dict = {'sequence': rna_to_aa(rna_seq), 'results': results}
	with open(out_fn, 'w') as f:
		f.write(dumps(result_dict, indent=4))


def main(args):
	in_fn = args[0]
	out_fn = args[1]
	params = _read_params(in_fn)
	rna_seq = params['sequence']
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
	_write_result(out_fn, rna_seq, results)


if __name__ == '__main__':
	main(sys.argv[1:])
