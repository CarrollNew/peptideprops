import abc

from sequence_property import (
	AbsorptionCoefficient,
	AliphaticIndex,
	ExtinctionCoefficient,
	FoldingMFE,
	GRAVY,
	HalfLife,
	InstabilityIndex,
	IsoelectricPoint,
	MolecularWeight,
	SequenceComposition
)


class SequencePropertyWrapper(object):
	__metaclass__ = abc.ABCMeta


class ProteinPropertyWrapper(SequencePropertyWrapper):
	def __init__(self, seq):
		self.seq = seq

	def absorption_coefficient(self):
		return [
			{
				'property': 'absorption_coefficient',
				'description': 'assuming all pairs of Cys residues form cystines',
				'name': 'Abs 0.1%',
				'unit': '',
				'value': AbsorptionCoefficient(self.seq).calculate_prop()
			},
			{
				'property': 'absorption_coefficient',
				'description': 'assuming all pairs of Cys residues are reduced',
				'name': 'Abs 0.1%',
				'unit': '',
				'value': AbsorptionCoefficient(self.seq).calculate_prop(cys_reduced=True)
			}
		]

	def aliphatic_index(self):
		return [{
			'property': 'aliphatic_index',
			'description': '',
			'name': 'Aliphatic index',
			'unit': '',
			'value': AliphaticIndex(self.seq).calculate_prop()
		}]

	def composition(self):
		res = SequenceComposition(self.seq).calculate_prop()
		return [
			{
				'property': 'aa_count',
				'description': '',
				'name': 'Amino acid composition',
				'unit': '',
				'value': res.aa_count
			},
			{
				'property': 'aa_len',
				'description': '',
				'name': 'Number of amino acids',
				'unit': '',
				'value': res.aa_len
			},
			{
				'property': 'atom_count',
				'description': '',
				'name': 'Atomic composition',
				'unit': '',
				'value': res.atom_count
			},
			{
				'property': 'atom_num',
				'description': '',
				'name': 'Total number of atoms',
				'unit': '',
				'value': res.atom_num
			},
			{
				'property': 'neg_charged_res',
				'description': '',
				'name': 'Number of negatively charged residues (Asp + Glu)',
				'unit': '',
				'value': res.aa_count['D'] + res.aa_count['E']
			},
			{
				'property': 'pos_charged_res',
				'description': '',
				'name': 'Number of positively charged residues (Arg + Lys)',
				'unit': '',
				'value': res.aa_count['R'] + res.aa_count['K']
			}
		]

	def extinction_coefficient(self):
		return [
			{
				'property': 'extinction_coefficient',
				'description': 'assuming all pairs of Cys residues form cystines',
				'name': 'Extinction coefficient',
				'unit': '1 / M / cm',
				'value': ExtinctionCoefficient(self.seq).calculate_prop()
			},
			{
				'property': 'extinction_coefficient',
				'description': 'assuming all pairs of Cys residues are reduced',
				'name': 'Extinction coefficient',
				'unit': '1 / M / cm',
				'value': ExtinctionCoefficient(self.seq).calculate_prop(cys_reduced=True)
			}
		]

	def gravy(self):
		return [{
			'property': 'gravy',
			'description': '',
			'name': 'Grand average of hydropathicity (GRAVY)',
			'unit': '',
			'value': GRAVY(self.seq).calculate_prop()
		}]

	def half_life(self):
		headers = ['Mammalian reticulocytes, in vitro', 'Yeast, in vivo', 'E. coli, in vivo']
		return [{
			'property': 'half_life',
			'description': '',
			'name': 'Estimated half-life',
			'unit': '',
			'value': {headers[i]: HalfLife(self.seq).calculate_prop(i) for i in xrange(len(headers))}
		}]

	def instability_index(self):
		return [{
			'property': 'instability_index',
			'description': '',
			'name': 'Instability index',
			'unit': '',
			'value': InstabilityIndex(self.seq).calculate_prop()
		}]

	def isoelectric_point(self):
		return [{
			'property': 'isoelectric_point',
			'description': '',
			'name': 'Isoelectric point',
			'unit': '',
			'value': IsoelectricPoint(self.seq).calculate_prop()
		}]

	def molecular_weight(self):
		return [
			{
				'property': 'molecular_weight',
				'description': 'at average resolution',
				'name': 'Molecular weight',
				'unit': 'g / mol',
				'value': MolecularWeight(self.seq).calculate_prop()
			},
			{
				'property': 'molecular_weight',
				'description': 'at monoisotopic resolution',
				'name': 'Molecular weight',
				'unit': 'g / mol',
				'value': MolecularWeight(self.seq).calculate_prop(monoisotopic=True)
			}
		]

	def sequence(self):
		return [{
			'property': 'sequence',
			'value': self.seq
		}]


class RNAPropertyWrapper(SequencePropertyWrapper):
	def __init__(self, sequences):
		self.sequences = sequences

	def folding(self, temp=37, be_concise=False):
		res = FoldingMFE(self.sequences).calculate_prop(temp)
		if be_concise:
			return [{'folding': item.folding, 'energy': item.energy} for item in res]
		else:
			assert len(res) == 1
			return [
				{
					'property': 'folding',
					'description': '',
					'name': 'RNA structure in bracket format',
					'unit': '',
					'value': res[0].folding
				},
				{
					'property': 'energy',
					'description': '',
					'name': 'RNA minimum free energy',
					'unit': 'kcal / mol',
					'value': res[0].energy
				}
			]
