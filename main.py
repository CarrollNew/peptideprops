import sys

from json import dumps, loads

from sequence_property_wrapper import ProteinPropertyWrapper, RNAPropertyWrapper
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
	if 'folding_temp' not in params:
		return {'error': 'Invalid parameters', 'message': 'No settings provided.'}
	rna_seq = params['sequence']
	results = []
	try:
		aa_seq = rna_to_aa(rna_seq)
	except StandardError as e:
		return {'error': 'Invalid mRNA sequence.', 'message': e.message}

	protein_analysis = ProteinPropertyWrapper(aa_seq)
	results.extend(protein_analysis.molecular_weight())
	results.extend(protein_analysis.composition())
	results.extend(protein_analysis.extinction_coefficient())
	results.extend(protein_analysis.absorption_coefficient())
	results.extend(protein_analysis.half_life())
	results.extend(protein_analysis.instability_index())
	results.extend(protein_analysis.isoelectric_point())
	results.extend(protein_analysis.aliphatic_index())
	results.extend(protein_analysis.gravy())
	rna_analysis = RNAPropertyWrapper([rna_seq])
	results.extend(rna_analysis.folding(temp=params['folding_temp']))
	result_dict = {'sequence': aa_seq, 'results': results}
	return result_dict


def process_batch_folding(params):
	if 'sequences' not in params:
		return {'error': 'Invalid parameters', 'message': 'No mRNA sequences provided.'}
	if 'folding_temp' not in params:
		return {'error': 'Invalid parameters', 'message': 'No settings provided.'}
	sequences = params['sequences']
	temp = params['folding_temp']
	return RNAPropertyWrapper(sequences).folding(temp=temp, be_concise=True)


def main(args):
	in_fn = args[0]
	out_fn = args[1]
	result = process(_read_params(in_fn))
	write_result(out_fn, result)


if __name__ == '__main__':
	main(sys.argv[1:])
