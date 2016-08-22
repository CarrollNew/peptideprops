from sequence_property import (
	ExtinctionCoefficient, IsoelectricPoint, InstabilityIndex, MolecularWeight
	)
from sequence_property import FoldingMFE


class PeptidePropsFactory(object):
	param_objects = {
		'extinction_coefficient': ExtinctionCoefficient,
		'isoelectric_point': IsoelectricPoint,
		'instability_index': InstabilityIndex,
		'molecular_weight': MolecularWeight,
		'folding': FoldingMFE
	}

	code_objects = {
		'ext_coef': 'extinction_coefficient',
		'pHI': 'isoelectric_point',
		'II': 'instability_index',
		'mw': 'molecular_weight',
		'folding': 'folding'
	}

	@staticmethod
	def from_prop_code(code):
		for code_object in PeptidePropsFactory.code_objects:
			if code == code_object:
				return PeptidePropsFactory.param_objects[PeptidePropsFactory.code_objects[code]]
		raise KeyError('Unknown property code: {}.'.format(code))

	@staticmethod
	def from_long_name(name):
		for param_object in PeptidePropsFactory.param_objects:
			if name == param_object:
				return PeptidePropsFactory.param_objects[name]
		raise KeyError('Unknown property name: {}.'.format(name))
