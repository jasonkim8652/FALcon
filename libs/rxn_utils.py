from rdkit import Chem
from rdkit.Chem.rdDeprotect import Deprotect


def get_template(
		df_template,
		reaction_idx,
	):
	if reaction_idx < len(df_template):
		row = df_template.iloc[reaction_idx]
		rxn = str(row['Rxn'])
		template = str(row['Template'])
	else:
		print ("Wrong reaction index for choosing reaction template")
		exit(-1)
	return rxn, template
		

def run_substitution(
		rxn,
		reactants,
	):
	product_list = []
	try:
		products = rxn.RunReactants(reactants)
		for product in products:
			smi_ = Chem.MolToSmiles(product[0])
			if smi_ not in product_list:
				product_list.append(smi_)
				#product_list.append(Chem.MolToSmiles(Deprotect(product[0])))
	except:
		print ("Error")
	return product_list
