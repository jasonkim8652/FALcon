#!/usr/bin/env python3

import os
import time
import argparse

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDeprotect import Deprotect
from rdkit.Chem import Draw

from libs.rxn_utils import get_template
from libs.rxn_utils import run_substitution

from libs.utils import canonicalized_mol
from libs.utils import str2bool


def main(
		args,
		df_template,
		df_bb,
		reaction_idx,
		revert_,
	):
	dir_ = os.path.abspath(args.dir_)

	mol_init = canonicalized_mol(args.smiles)

	rxn_name, template = get_template(
		df_template=df_template,
		reaction_idx=reaction_idx,
	)
	print (rxn_name, "is chosen for proceeding reactions")
	print ("Number of building blocks:", len(df_bb))
	rxn = AllChem.ReactionFromSmarts(template)

	product_list = []
	id_bb_list = []
	for idx in range(len(df_bb)):
		row = df_bb.iloc[idx]
		id_bb = row['ID']
		smi_bb = row['SMILES']

		product = []
		if '.' not in smi_bb:
			mol_bb = Chem.MolFromSmiles(smi_bb)
			reactants = (mol_init, mol_bb)
			if str2bool(revert_):
				reactants = (mol_bb, mol_init)
			products = run_substitution(rxn, reactants)

			for product in products:
				mol_ = Chem.MolFromSmiles(product)
				mol_ = Deprotect(mol_)
				smi_ = Chem.MolToSmiles(mol_)
				product_list.append(smi_)
				id_bb_list.append(id_bb)

	id_product_list = [args.prefix + '_' + str(i).rjust(5, '0') for i in range(len(product_list))]
	df_product = pd.DataFrame({})
	df_product['ID'] = id_product_list
	df_product['SMILES'] = product_list
	df_product['ID_BB'] = id_bb_list

	df_product = df_product.drop_duplicates(subset=['SMILES'])

	csv_path = os.path.join(dir_, args.output+'.csv')
	df_product.to_csv(csv_path, index=False)
	print ("Finished")


if __name__ == '__main__':
	template_path = '/home/seongok/works/design_modules/data/rxn_gen_template.csv'
	df_template = pd.read_csv(template_path)
	print ("List of reaction templates implemented:")
	print (df_template)
	reaction_idx = int(input("Index of reaction to use: "))
	print ("-"*95)
	print ("\n")

	bb_list_path = '/home/seongok/works/design_modules/data/rxn_building_block.csv'
	df_bb = pd.read_csv(bb_list_path)
	print ("List of Building block files:")
	print (df_bb)
	bb_idx = input("Index of building block file to use: ")
	print ("-"*95)
	print ("\n")

	bb_name = str(df_bb.iloc[int(bb_idx)]['Name'])
	bb_path = '/home/seongok/works/design_modules/data/'+bb_name+'.csv'
	df_bb = pd.read_csv(bb_path)

	revert_ = input("Whether to revert the reaction template: ")

	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--smiles', type=str, required=True, 
						help='SMILES of query molecule')
	parser.add_argument('-o', '--output', type=str, default='test', 
						help='Output CSV File')
	parser.add_argument('-p', '--prefix', type=str, default='RXGN', 
						help='Prefix for ID generation')
	parser.add_argument('-d', '--dir_', type=str, default=os.getcwd(), 
						help='Directory to save the files')
	'''
	parser.add_argument('-r', '--revert', type=str2bool, default=False, 
						help='Whether to revert the chosen template')
	parser.add_argument('-b', '--bb_path', type=str, required=True, 
						help='SMILES of core scaffold')
	parser.add_argument('-i', '--reaction_idx', type=int, default=reaction_idx, 
						help='Reaction idx in the template file')
	parser.add_argument('--group1', type=str, required=True, 
						help='Functional group to be reacted in the anchor')
	parser.add_argument('--group2', type=str, required=True, 
						help='Functional group to be reacted in other BBs')
	'''
	args = parser.parse_args()

	main(
		args=args,
		df_template=df_template,
		df_bb=df_bb,
		reaction_idx=reaction_idx,
		revert_=revert_,
	)
