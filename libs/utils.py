import argparse

from rdkit import Chem


def sanitize_and_smiles(mol):
	mol = Chem.RemoveHs(mol)
	smi = Chem.MolToSmiles(mol)
	mol = Chem.MolFromSmiles(smi)
	smi = Chem.MolToSmiles(mol)
	return smi


def canonicalized_mol(smi):
	mol = Chem.MolFromSmiles(smi)
	smi = Chem.MolToSmiles(mol)
	mol = Chem.MolFromSmiles(smi)
	return mol


def str2bool(v):
	v = str(v)
	if v.lower() in ['yes', 'true', 't', 'y', '1']:
		return True
	elif v.lower() in ['no', 'false', 'f', 'n', '0']:
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected')
