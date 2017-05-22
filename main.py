import sys

from json import dumps, loads

from sequence_property_wrapper import ProteinPropertyWrapper, RNAPropertyWrapper
from seq_utils import dna_to_rna, rna_to_aa, AA_LETTERS, RNA_LETTERS, DNA_LETTERS


def _read_params(fn):
	with open(fn) as f:
		return loads(f.read())


def _write_result(out_fn, result):
	with open(out_fn, 'w') as f:
		f.write(dumps(result, indent=4))


def process(params):
	sequence = params.get('sequence', None)
	sequence_type = params.get('type', None)
	folding_temp = params.get('folding_temp', 37)
	if not sequence:
		return {'error': 'Invalid parameters', 'message': 'No mRNA sequence provided.'}
	if not sequence_type:
		return {'error': 'Invalid parameters', 'message': 'No sequence type provided.'}
	results = {}
	if sequence_type == 'mRNA':
		if all(map(lambda let: let in DNA_LETTERS, sequence)):
			sequence = dna_to_rna(sequence)
		elif any(map(lambda let: let not in RNA_LETTERS, sequence)):
			return {'error': 'Not an RNA sequence.'}
		rna_results = []
		rna_results = RNAPropertyWrapper([sequence]).folding(temp=folding_temp)
		results[sequence_type] = rna_results
		try:
			sequence = rna_to_aa(sequence)
			sequence_type = 'Protein'
		except StandardError as e:
			return results
	if sequence_type == 'Protein':
		if not all(map(lambda let: let in AA_LETTERS, sequence)):
			return {'error': 'Not a protein sequence.'}
		protein_results = []
		protein_analysis = ProteinPropertyWrapper(sequence)
		protein_results.extend(protein_analysis.molecular_weight())
		protein_results.extend(protein_analysis.composition())
		protein_results.extend(protein_analysis.extinction_coefficient())
		protein_results.extend(protein_analysis.absorption_coefficient())
		protein_results.extend(protein_analysis.half_life())
		protein_results.extend(protein_analysis.instability_index())
		protein_results.extend(protein_analysis.isoelectric_point())
		protein_results.extend(protein_analysis.aliphatic_index())
		protein_results.extend(protein_analysis.gravy())
		protein_results.extend(protein_analysis.sequence())
		results[sequence_type] = protein_results
		return results
	return {'error': 'Invalid parameters', 'message': 'Invalid sequence type.'}


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
	_write_result(out_fn, result)


if __name__ == '__main__':
	main(sys.argv[1:])
