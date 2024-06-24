#!/usr/bin/env python3

import os

import pandas as pd
import argparse

from rdkit import Chem
from rdkit.Chem import Draw

from libs.decoration_utils import add_groups
from libs.decoration_utils import edit_initial_smiles
from libs.decoration_utils import get_num_attachment
from libs.decoration_utils import get_branch_list

from libs.utils import sanitize_and_smiles
from libs.utils import str2bool


def main(args):
	dir_ = os.path.abspath(args.dir_)

	smi_init = edit_initial_smiles(args.smiles)
	num_attachment = get_num_attachment(smi_init)
	
	df_rgroup = pd.read_csv(args.rgroup_data)
	df_recipe = pd.read_csv(args.recipe)
	if len(df_recipe) != num_attachment:
		print ("Number of attachment point is different to number of R-groups in the recipe")
		exit(-1)
	
	branch_list = get_branch_list(
		df_rgroup=df_rgroup,
		df_recipe=df_recipe,
		combi=args.combi
	)

	mol_list = []
	print ("Number of branches to add:", len(branch_list))
	for branch in branch_list:
		try:
			mol = add_groups(Chem.MolFromSmiles(smi_init), branch)
			mol_list.append(mol)
		except:
			pass

	if len(mol_list) > 0:
		id_list = [args.prefix+'_'+str(i).rjust(3, '0') for i in range(len(mol_list))]
		smi_list = [sanitize_and_smiles(mol) for mol in mol_list]

		csv_path = os.path.join(dir_, args.output+'.csv')
		df_output = pd.DataFrame({})
		df_output['ID'] = id_list
		df_output['SMILES'] = smi_list
		df_output.to_csv(csv_path, index=False) 

		img_path = os.path.join(dir_, args.output+'.png')
		img = Draw.MolsToGridImage(mol_list, legends=id_list, molsPerRow=4)
		img.save(img_path)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--smiles', type=str, required=True, 
						help='SMILES of query molecule')
	parser.add_argument('-r', '--recipe', type=str, required=True, 
						help='Recipe')
	parser.add_argument('-c', '--combi', type=str2bool, default=False, 
						help='Whether to make all combinations')
	parser.add_argument('-o', '--output', type=str, required=True, 
						help='Output CSV File')
	parser.add_argument('-p', '--prefix', type=str, default='HOP', 
						help='Prefix for ID generation')
	parser.add_argument('-d', '--dir_', type=str, default=os.getcwd(), 
						help='Directory to save the files')

	data_path = '/home/seongok/works/design_modules/data/rgroup_list.csv'
	parser.add_argument('--rgroup_data', type=str, default=data_path, 
						help='')
	args = parser.parse_args()

	main(args)
