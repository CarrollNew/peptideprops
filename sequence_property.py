import os
import abc
import subprocess
import tempfile
from collections import namedtuple

from seq_utils import rna_to_aa, AA_LETTERS, RNA_LETTERS
from seq_utils import aa_composition_table, aa_codes_table, atom_codes_table, aa_kd_table
from seq_utils import aa_half_life


FoldingResult = namedtuple('FoldingResult', ['folding', 'energy'])
CompositionResult = namedtuple('CompositionResult', ['aa_count', 'aa_len', 'atom_count', 'atom_num'])


class SequenceProperty(object):
	__metaclass__ = abc.ABCMeta

	def calculate_prop():
		pass


class PeptideProperty(SequenceProperty):
	__metaclass__ = abc.ABCMeta

	def calculate_prop():
		pass


class RNAProperty(SequenceProperty):
	__metaclass__ = abc.ABCMeta

	def calculate_prop():
		pass


class AbsorptionCoefficient(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq

	def calculate_prop(self, monoisotopic=False, cys_reduced=False):
		mw = MolecularWeight(self.seq).calculate_prop(monoisotopic)
		ext_coef = ExtinctionCoefficient(self.seq).calculate_prop(cys_reduced)
		return ext_coef / mw


class MolecularWeight(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq

	def calculate_prop(self, monoisotopic=False):
		weights = MolecularWeight._get_weights_table(monoisotopic)
		water_weight = MolecularWeight._get_water_weight(monoisotopic)
		return (sum(map(lambda letter: self.seq.count(letter) * weights[letter], AA_LETTERS))
			- (len(self.seq) - 1) * water_weight)

	@staticmethod
	def _get_weights_table(monoisotopic=False):
		def _get_monoisotopic_weights():
			return {
			    "A": 89.047678, 
			    "C": 121.019749, 
			    "E": 147.053158, 
			    "D": 133.037508, 
			    "G": 75.032028, 
			    "F": 165.078979, 
			    "I": 131.094629, 
			    "H": 155.069477, 
			    "K": 146.105528, 
			    "M": 149.051049, 
			    "L": 131.094629, 
			    "O": 255.158292, 
			    "N": 132.053492, 
			    "Q": 146.069142, 
			    "P": 115.063329, 
			    "S": 105.042593, 
			    "R": 174.111676, 
			    "U": 168.964203, 
			    "T": 119.058243, 
			    "W": 204.089878, 
			    "V": 117.078979, 
			    "Y": 181.073893
			}

		def _get_weights():
			return {
			    "A": 89.0932, 
			    "C": 121.1582, 
			    "E": 147.1293, 
			    "D": 133.1027, 
			    "G": 75.0666, 
			    "F": 165.1891, 
			    "I": 131.1729, 
			    "H": 155.1546, 
			    "K": 146.1876, 
			    "M": 149.2113, 
			    "L": 131.1729, 
			    "O": 255.3134, 
			    "N": 132.1179, 
			    "Q": 146.1445, 
			    "P": 115.1305, 
			    "S": 105.0926, 
			    "R": 174.201, 
			    "U": 168.0532, 
			    "T": 119.1192, 
			    "W": 204.2252, 
			    "V": 117.1463, 
			    "Y": 181.1885
			}

		return _get_monoisotopic_weights() if monoisotopic else _get_weights()

	@staticmethod
	def _get_water_weight(monoisotopic=False):
		return 18.010565 if monoisotopic else 18.0153


class InstabilityIndex(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq

	def calculate_prop(self):
		pairwise_weights = InstabilityIndex._get_dipeptide_weights()
		sum_weights = sum(map(lambda i: pairwise_weights[self.seq[i:i+2]], xrange(len(self.seq) - 1)))
		return 10.0 * sum_weights / len(self.seq)

	@staticmethod
	def _get_dipeptide_weights():
		return {
		    "GW": 13.34, 
		    "GV": 1.0, 
		    "GT": -7.49, 
		    "GS": 1.0, 
		    "GR": 1.0, 
		    "GQ": 1.0, 
		    "GP": 1.0, 
		    "GY": -7.49, 
		    "GG": 13.34, 
		    "GF": 1.0, 
		    "GE": -6.540000000000001, 
		    "GD": 1.0, 
		    "GC": 1.0, 
		    "GA": -7.49, 
		    "GN": -7.49, 
		    "GM": 1.0, 
		    "GL": 1.0, 
		    "GK": -7.49, 
		    "GI": -7.49, 
		    "GH": 1.0, 
		    "ME": 1.0, 
		    "MD": 1.0, 
		    "MG": 1.0, 
		    "MF": 1.0, 
		    "MA": 13.34, 
		    "MC": 1.0, 
		    "MM": -1.8799999999999997, 
		    "ML": 1.0, 
		    "MN": 1.0, 
		    "MI": 1.0, 
		    "MH": 58.279999999999994, 
		    "MK": 1.0, 
		    "MT": -1.8799999999999997, 
		    "MW": 1.0, 
		    "MV": 1.0, 
		    "MQ": -6.540000000000001, 
		    "MP": 44.94, 
		    "MS": 44.94, 
		    "MR": -6.540000000000001, 
		    "MY": 24.68, 
		    "FP": 20.26, 
		    "FQ": 1.0, 
		    "FR": 1.0, 
		    "FS": 1.0, 
		    "FT": 1.0, 
		    "FV": 1.0, 
		    "FW": 1.0, 
		    "FY": 33.601, 
		    "FA": 1.0, 
		    "FC": 1.0, 
		    "FD": 13.34, 
		    "FE": 1.0, 
		    "FF": 1.0, 
		    "FG": 1.0, 
		    "FH": 1.0, 
		    "FI": 1.0, 
		    "FK": -14.029999999999998, 
		    "FL": 1.0, 
		    "FM": 1.0, 
		    "FN": 1.0, 
		    "SY": 1.0, 
		    "SS": 20.26, 
		    "SR": 20.26, 
		    "SQ": 20.26, 
		    "SP": 44.94, 
		    "SW": 1.0, 
		    "SV": 1.0, 
		    "ST": 1.0, 
		    "SK": 1.0, 
		    "SI": 1.0, 
		    "SH": 1.0, 
		    "SN": 1.0, 
		    "SM": 1.0, 
		    "SL": 1.0, 
		    "SC": 33.6, 
		    "SA": 1.0, 
		    "SG": 1.0, 
		    "SF": 1.0, 
		    "SE": 20.26, 
		    "SD": 1.0, 
		    "YI": 1.0, 
		    "YH": 13.34, 
		    "YK": 1.0, 
		    "YM": 44.94, 
		    "YL": 1.0, 
		    "YN": 1.0, 
		    "YA": 24.68, 
		    "YC": 1.0, 
		    "YE": -6.540000000000001, 
		    "YD": 24.68, 
		    "YG": -7.49, 
		    "YF": 1.0, 
		    "YY": 13.34, 
		    "YQ": 1.0, 
		    "YP": 13.34, 
		    "YS": 1.0, 
		    "YR": -15.91, 
		    "YT": -7.49, 
		    "YW": -9.37, 
		    "YV": 1.0, 
		    "LF": 1.0, 
		    "LG": 1.0, 
		    "LD": 1.0, 
		    "LE": 1.0, 
		    "LC": 1.0, 
		    "LA": 1.0, 
		    "LN": 1.0, 
		    "LL": 1.0, 
		    "LM": 1.0, 
		    "LK": -7.49, 
		    "LH": 1.0, 
		    "LI": 1.0, 
		    "LV": 1.0, 
		    "LW": 24.68, 
		    "LT": 1.0, 
		    "LR": 20.26, 
		    "LS": 1.0, 
		    "LP": 20.26, 
		    "LQ": 33.6, 
		    "LY": 1.0, 
		    "RT": 1.0, 
		    "RV": 1.0, 
		    "RW": 58.279999999999994, 
		    "RP": 20.26, 
		    "RQ": 20.26, 
		    "RR": 58.279999999999994, 
		    "RS": 44.94, 
		    "RY": -6.540000000000001, 
		    "RD": 1.0, 
		    "RE": 1.0, 
		    "RF": 1.0, 
		    "RG": -7.49, 
		    "RA": 1.0, 
		    "RC": 1.0, 
		    "RL": 1.0, 
		    "RM": 1.0, 
		    "RN": 13.34, 
		    "RH": 20.26, 
		    "RI": 1.0, 
		    "RK": 1.0, 
		    "VH": 1.0, 
		    "VI": 1.0, 
		    "EM": 1.0, 
		    "EL": 1.0, 
		    "EN": 1.0, 
		    "EI": 20.26, 
		    "EH": -6.540000000000001, 
		    "EK": 1.0, 
		    "EE": 33.6, 
		    "ED": 20.26, 
		    "EG": 1.0, 
		    "EF": 1.0, 
		    "EA": 1.0, 
		    "EC": 44.94, 
		    "VM": 1.0, 
		    "EY": 1.0, 
		    "VN": 1.0, 
		    "ET": 1.0, 
		    "EW": -14.029999999999998, 
		    "EV": 1.0, 
		    "EQ": 20.26, 
		    "EP": 20.26, 
		    "ES": 20.26, 
		    "ER": 1.0, 
		    "VP": 20.26, 
		    "VQ": 1.0, 
		    "VR": 1.0, 
		    "VT": -7.49, 
		    "VW": 1.0, 
		    "KC": 1.0, 
		    "KA": 1.0, 
		    "KG": -7.49, 
		    "KF": 1.0, 
		    "KE": 1.0, 
		    "KD": 1.0, 
		    "KK": 1.0, 
		    "KI": -7.49, 
		    "KH": 1.0, 
		    "KN": 1.0, 
		    "KM": 33.6, 
		    "KL": -7.49, 
		    "KS": 1.0, 
		    "KR": 33.6, 
		    "KQ": 24.64, 
		    "KP": -6.540000000000001, 
		    "KW": 1.0, 
		    "KV": -7.49, 
		    "KT": 1.0, 
		    "KY": 1.0, 
		    "DN": 1.0, 
		    "DL": 1.0, 
		    "DM": 1.0, 
		    "DK": -7.49, 
		    "DH": 1.0, 
		    "DI": 1.0, 
		    "DF": -6.540000000000001, 
		    "DG": 1.0, 
		    "DD": 1.0, 
		    "DE": 1.0, 
		    "DC": 1.0, 
		    "DA": 1.0, 
		    "DY": 1.0, 
		    "DV": 1.0, 
		    "DW": 1.0, 
		    "DT": -14.029999999999998, 
		    "DR": -6.540000000000001, 
		    "DS": 20.26, 
		    "DP": 1.0, 
		    "DQ": 1.0, 
		    "QQ": 20.26, 
		    "QP": 20.26, 
		    "QS": 44.94, 
		    "QR": 1.0, 
		    "QT": 1.0, 
		    "QW": 1.0, 
		    "QV": -6.540000000000001, 
		    "QY": -6.540000000000001, 
		    "QA": 1.0, 
		    "QC": -6.540000000000001, 
		    "QE": 20.26, 
		    "QD": 20.26, 
		    "QG": 1.0, 
		    "QF": -6.540000000000001, 
		    "QI": 1.0, 
		    "QH": 1.0, 
		    "QK": 1.0, 
		    "QM": 1.0, 
		    "QL": 1.0, 
		    "QN": 1.0, 
		    "WG": -9.37, 
		    "WF": 1.0, 
		    "WE": 1.0, 
		    "WD": 1.0, 
		    "WC": 1.0, 
		    "WA": -14.029999999999998, 
		    "WN": 13.34, 
		    "WM": 24.68, 
		    "WL": 13.34, 
		    "WK": 1.0, 
		    "WI": 1.0, 
		    "WH": 24.68, 
		    "WW": 1.0, 
		    "WV": -7.49, 
		    "WT": -14.029999999999998, 
		    "WS": 1.0, 
		    "WR": 1.0, 
		    "WQ": 1.0, 
		    "WP": 1.0, 
		    "WY": 1.0, 
		    "PR": -6.540000000000001, 
		    "PS": 20.26, 
		    "PP": 20.26, 
		    "PQ": 20.26, 
		    "PV": 20.26, 
		    "PW": -1.8799999999999997, 
		    "PT": 1.0, 
		    "PY": 1.0, 
		    "PC": -6.540000000000001, 
		    "PA": 20.26, 
		    "PF": 20.26, 
		    "PG": 1.0, 
		    "PD": -6.540000000000001, 
		    "PE": 18.38, 
		    "PK": 1.0, 
		    "PH": 1.0, 
		    "PI": 1.0, 
		    "PN": 1.0, 
		    "PL": 1.0, 
		    "PM": -6.540000000000001, 
		    "CK": 1.0, 
		    "CI": 1.0, 
		    "CH": 33.6, 
		    "CN": 1.0, 
		    "CM": 33.6, 
		    "CL": 20.26, 
		    "CC": 1.0, 
		    "CA": 1.0, 
		    "CG": 1.0, 
		    "CF": 1.0, 
		    "CE": 1.0, 
		    "CD": 20.26, 
		    "CY": 1.0, 
		    "CS": 1.0, 
		    "CR": 1.0, 
		    "CQ": -6.540000000000001, 
		    "CP": 20.26, 
		    "CW": 24.68, 
		    "CV": -6.540000000000001, 
		    "CT": 33.6, 
		    "IY": 1.0, 
		    "VA": 1.0, 
		    "VC": 1.0, 
		    "VD": -14.029999999999998, 
		    "VE": 1.0, 
		    "VF": 1.0, 
		    "VG": -7.49, 
		    "IQ": 1.0, 
		    "IP": -1.8799999999999997, 
		    "IS": 1.0, 
		    "IR": 1.0, 
		    "VL": 1.0, 
		    "IT": 1.0, 
		    "IW": 1.0, 
		    "IV": -7.49, 
		    "II": 1.0, 
		    "IH": 13.34, 
		    "IK": -7.49, 
		    "VS": 1.0, 
		    "IM": 1.0, 
		    "IL": 20.26, 
		    "VV": 1.0, 
		    "IN": 1.0, 
		    "IA": 1.0, 
		    "VY": -6.540000000000001, 
		    "IC": 1.0, 
		    "IE": 44.94, 
		    "ID": 1.0, 
		    "IG": 1.0, 
		    "IF": 1.0, 
		    "HY": 44.94, 
		    "HR": 1.0, 
		    "HS": 1.0, 
		    "HP": -1.8799999999999997, 
		    "HQ": 1.0, 
		    "HV": 1.0, 
		    "HW": -1.8799999999999997, 
		    "HT": -6.540000000000001, 
		    "HK": 24.68, 
		    "HH": 1.0, 
		    "HI": 44.94, 
		    "HN": 24.68, 
		    "HL": 1.0, 
		    "HM": 1.0, 
		    "HC": 1.0, 
		    "HA": 1.0, 
		    "HF": -9.37, 
		    "HG": -9.37, 
		    "HD": 1.0, 
		    "HE": 1.0, 
		    "NH": 1.0, 
		    "NI": 44.94, 
		    "NK": 24.68, 
		    "NL": 1.0, 
		    "NM": 1.0, 
		    "NN": 1.0, 
		    "NA": 1.0, 
		    "NC": -1.8799999999999997, 
		    "ND": 1.0, 
		    "NE": 1.0, 
		    "NF": -14.029999999999998, 
		    "NG": -14.029999999999998, 
		    "NY": 1.0, 
		    "NP": -1.8799999999999997, 
		    "NQ": -6.540000000000001, 
		    "NR": 1.0, 
		    "NS": 1.0, 
		    "NT": -7.49, 
		    "NV": 1.0, 
		    "NW": -9.37, 
		    "TY": 1.0, 
		    "TV": 1.0, 
		    "TW": -14.029999999999998, 
		    "TT": 1.0, 
		    "TR": 1.0, 
		    "TS": 1.0, 
		    "TP": 1.0, 
		    "TQ": -6.540000000000001, 
		    "TN": -14.029999999999998, 
		    "TL": 1.0, 
		    "TM": 1.0, 
		    "TK": 1.0, 
		    "TH": 1.0, 
		    "TI": 1.0, 
		    "TF": 13.34, 
		    "TG": -7.49, 
		    "TD": 1.0, 
		    "TE": 20.26, 
		    "TC": 1.0, 
		    "TA": 1.0, 
		    "AA": 1.0, 
		    "AC": 44.94, 
		    "AE": 1.0, 
		    "AD": -7.49, 
		    "AG": 1.0, 
		    "AF": 1.0, 
		    "AI": 1.0, 
		    "AH": -7.49, 
		    "AK": 1.0, 
		    "AM": 1.0, 
		    "AL": 1.0, 
		    "AN": 1.0, 
		    "AQ": 1.0, 
		    "AP": 20.26, 
		    "AS": 1.0, 
		    "AR": 1.0, 
		    "AT": 1.0, 
		    "AW": 1.0, 
		    "AV": 1.0, 
		    "AY": 1.0, 
		    "VK": -1.8799999999999997
		}


class IsoelectricPoint(PeptideProperty):
	positive_pKs = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98}
	negative_pKs = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0}
	pKcterminal = {'D': 4.55, 'E': 4.75}
	pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.7}
	charged_aas = ('K', 'R', 'H', 'D', 'E', 'C', 'Y')

	def __init__(self, seq):
		self.seq = seq
		self.aa_count = SequenceComposition(seq).calculate_prop().aa_count
		self.charged_aas_content = self._select_charged(self.aa_count)

	def calculate_prop(self):
		pos_pKs = dict(IsoelectricPoint.positive_pKs)
		neg_pKs = dict(IsoelectricPoint.negative_pKs)
		nterm = self.seq[0]
		cterm = self.seq[-1]
		if nterm in IsoelectricPoint.pKnterminal:
			pos_pKs['Nterm'] = IsoelectricPoint.pKnterminal[nterm]
		if cterm in IsoelectricPoint.pKcterminal:
			neg_pKs['Cterm'] = IsoelectricPoint.pKcterminal[cterm]

        # Bracket between pH1 and pH2
		pH = 7.0
		Charge = self._chargeR(pH, pos_pKs, neg_pKs)
		if Charge > 0.0:
			pH1 = pH
			Charge1 = Charge
			while Charge1 > 0.0:
				pH = pH1 + 1.0
				Charge = self._chargeR(pH, pos_pKs, neg_pKs)
				if Charge > 0.0:
					pH1 = pH
					Charge1 = Charge
				else:
					pH2 = pH
					Charge2 = Charge
					break
		else:
			pH2 = pH
			Charge2 = Charge
			while Charge2 < 0.0:
				pH = pH2 - 1.0
				Charge = self._chargeR(pH, pos_pKs, neg_pKs)
				if Charge < 0.0:
					pH2 = pH
					Charge2 = Charge
				else:
					pH1 = pH
					Charge1 = Charge
					break

		# Bisection
		while pH2 - pH1 > 0.0001 and Charge != 0.0:
			pH = (pH1 + pH2) / 2.0
			Charge = self._chargeR(pH, pos_pKs, neg_pKs)
			if Charge > 0.0:
				pH1 = pH
				Charge1 = Charge
			else:
				pH2 = pH
				Charge2 = Charge
		return pH

	def _chargeR(self, pH, pos_pKs, neg_pKs):
		PositiveCharge = 0.0
		for aa, pK in pos_pKs.items():
			CR = 10 ** (pK - pH)
			partial_charge = CR / (CR + 1.0)
			PositiveCharge += self.charged_aas_content[aa] * partial_charge

		NegativeCharge = 0.0
		for aa, pK in neg_pKs.items():
			CR = 10 ** (pH - pK)
			partial_charge = CR / (CR + 1.0)
			NegativeCharge += self.charged_aas_content[aa] * partial_charge

		return PositiveCharge - NegativeCharge

	def _select_charged(self, AminoAcidsContent):
		charged = {}
		for aa in IsoelectricPoint.charged_aas:
			charged[aa] = float(AminoAcidsContent[aa])
		charged['Nterm'] = 1.0
		charged['Cterm'] = 1.0
		return charged


class ExtinctionCoefficient(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq
	
	def calculate_prop(self, cys_reduced=False):
		ext_params = ExtinctionCoefficient._get_ext_params()
		return (
			self.seq.count('Y') * ext_params[0] + self.seq.count('W') * ext_params[1]
			+ int(not cys_reduced) * self.seq.count('C') / 2 * ext_params[2]
		)

	@staticmethod
	def _get_ext_params(nm=280):
		"""
		[Tyr (Y), Trp (W), Cys (C)]
		"""
		params_table = {
			276: [1450, 5400, 145],
			278: [1400, 5600, 127],
			279: [1345, 5660, 120],
			280: [1490, 5500, 125],
			282: [1200, 5600, 100]
		}
		if nm in params_table:
			return params_table[nm]
		raise KeyError('Specified wavelength not found.')


class SequenceComposition(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq
	
	def calculate_prop(self):
		aa_count = {k: self.seq.count(k) for k in AA_LETTERS}
		aa_len = sum(aa_count.values())
		atom_count = {
			atom: sum(aa_composition_table[aa][atom] * aa_count[aa] for aa in aa_count)
			for atom in atom_codes_table
		}
		atom_count['H'] -= 2 * (aa_len - 1)
		atom_count['O'] -= (aa_len - 1)
		atom_num = sum(atom_count.values())
		return CompositionResult(aa_count, aa_len, atom_count, atom_num)


class GRAVY(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq

	def calculate_prop(self):
		return 1.0 * sum(aa_kd_table[aa] for aa in self.seq) / len(self.seq)


class AliphaticIndex(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq

	def calculate_prop(self):
		a = 2.9
		b = 3.9
		return 100.0 * (
			self.seq.count('A') + a * self.seq.count('V') + b * (self.seq.count('I') + self.seq.count('L'))
			) / len(self.seq)


class HalfLife(PeptideProperty):
	def __init__(self, seq):
		self.seq = seq

	def calculate_prop(self, cell_type):
		nterm = self.seq[0]
		assert cell_type in [0, 1, 2]
		return aa_half_life[nterm][cell_type]
		

class FoldingMFE(RNAProperty):
	possible_rnafold_paths = [
		'/usr/local/bin/RNAfold',
		'C://Program Files (x86)//ViennaRNA Package//RNAfold.exe',
		'/usr/local/bin/ViennaRNA/RNAfold',
		'/usr/bin/ViennaRNA/RNAfold',
		'/usr/bin/RNAfold'
	]
	
	def __init__(self, sequences):
		self.sequences = sequences

	def calculate_prop(self, temp):
		in_fasta = tempfile.gettempdir() + os.sep + 'rnafold_input.fasta'
		out_fasta = tempfile.gettempdir() + os.sep + 'rnafold_output'
		seq_names = self._write_input(in_fasta)
		command = [
			FoldingMFE._get_rnafold_app(), 
			'-i', in_fasta, 
			'-o', out_fasta, 
			'--noPS', 
			'-T', str(temp)
		]
		try:
			subprocess.check_output(command, stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError as e:
			print e.message
		return self._read_results(out_fasta, seq_names)

	def _write_input(self, fn):
		names = map(lambda i: 'sequence_{}'.format(i), xrange(len(self.sequences)))
		with open(fn, 'w') as f:
			for i in xrange(len(self.sequences)):
				f.write('>{}\n'.format(names[i]))
				f.write('{}\n'.format(self.sequences[i]))
		return names

	def _read_results(self, out_fasta, seq_names):
		return map(lambda name: self._parse_output(out_fasta + '_{}.fold'.format(name)), seq_names)

	def _parse_output(self, fn):
		with open(fn) as f:
			contents = f.read().splitlines()
		strs = contents[-1].split()
		folding = strs[0]
		energy = float(strs[-1][:-1].replace('(', '').replace(')', ''))
		return FoldingResult(folding, energy)

	@staticmethod
	def _get_rnafold_app():
		for app in FoldingMFE.possible_rnafold_paths:
			if os.path.exists(app):
				return app
		raise ValueError('Could not find suitable RNAfold app directory.')
