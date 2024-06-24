#!/usr/bin/env python3

import os
import argparse

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw

from libs.hopping_utils import check_mol
from libs.hopping_utils import HeteroShuffle
 
from libs.utils import canonicalized_mol


def main(args):
	dir_ = os.path.abspath(args.dir_)

	mol = canonicalized_mol(args.smiles)
	core = canonicalized_mol(args.core)

	hetero_shuffle = HeteroShuffle(mol, core)
	products = hetero_shuffle.generate_mols()
	print ("Total number of generated products:", len(products))

	if len(products) > 0:
		id_list = [args.prefix+'_'+str(i).rjust(3, '0') for i in range(len(products))]
		smi_list = [Chem.MolToSmiles(product) for product in products]

		df = pd.DataFrame({})
		df['ID'] = id_list
		df['SMILES'] = smi_list

		csv_path = os.path.join(dir_, args.output+'.csv')
		df.to_csv(csv_path, index=False)
	
		img_path = os.path.join(dir_, args.output+'.png')
		img = Draw.MolsToGridImage(products, legends=id_list, molsPerRow=4)
		img.save(img_path)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--smiles', type=str, required=True, 
						help='SMILES of query molecule')
	parser.add_argument('-c', '--core', type=str, required=True, 
						help='SMILES of core scaffold')
	parser.add_argument('-o', '--output', type=str, required=True, 
						help='Output CSV File')
	parser.add_argument('-p', '--prefix', type=str, default='HOP', 
						help='Prefix for ID generation')
	parser.add_argument('-d', '--dir_', type=str, default=os.getcwd(), 
						help='Directory to save the files')
	args = parser.parse_args()

	main(args)
