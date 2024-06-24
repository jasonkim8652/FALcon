from rdkit import Chem
from rdkit import DataStructs

from itertools import product


def add_groups(m, mapping):
    while True:
        # Find wildcard atom if available, otherwise exit
        a = None
        for a_ in m.GetAtoms():
            if a_.GetAtomicNum() == 0:
                a = a_
                break
        if a is None:
            break

        # Get appropriate group to substitute in from mapping
        group = mapping[a.GetAtomMapNum()]
        if group.GetNumAtoms() == 1:

            # Directly substitute atom in, if single atom group
            a.SetAtomicNum(group.GetAtomWithIdx(0).GetAtomicNum())
            a.SetAtomMapNum(0)
        else:

            # Set wildcard atoms to having AtomMapNum 1000 for tracking
            a.SetAtomMapNum(1000)

            for a_ in group.GetAtoms():
                if a_.GetAtomicNum() == 0:
                    a_.SetAtomMapNum(1000)

            # Put group and base molecule together and make it editable
            m = Chem.CombineMols(m, group)
            m = Chem.RWMol(m)

            # Find using tracking number the atoms to merge in new molecule
            a1 = None
            a2 = None
            bonds = []
            for a in m.GetAtoms():
                if a.GetAtomMapNum() == 1000:
                    if a1 is None:
                        a1 = a
                    else:
                        a2 = a

            # Find atoms to bind together based on atoms to merge
            b1 = a1.GetBonds()[0]
            start = (b1.GetBeginAtomIdx() if b1.GetEndAtomIdx() == a1.GetIdx()
                else b1.GetEndAtomIdx())

            b2 = a2.GetBonds()[0]
            end = (b2.GetBeginAtomIdx() if b2.GetEndAtomIdx() == a2.GetIdx()
                else b2.GetEndAtomIdx())

            # Add the connection and remove original wildcard atoms
            m.AddBond(start, end, order=Chem.rdchem.BondType.SINGLE)
            m.RemoveAtom(a1.GetIdx())
            m.RemoveAtom(a2.GetIdx())
    return m


def edit_initial_smiles(smi):
	smi = smi.replace('[A1]', '[*:1]')
	smi = smi.replace('[A2]', '[*:2]')
	smi = smi.replace('[A3]', '[*:3]')
	smi = smi.replace('[A4]', '[*:4]')
	mol = Chem.MolFromSmiles(smi)
	smi = Chem.MolToSmiles(mol)
	return smi


def get_num_attachment(smi):
	mol = Chem.MolFromSmiles(smi)
	n_ = 0
	for atom in mol.GetAtoms():
		if atom.GetAtomicNum() == 0:
			n_ += 1
	return n_


def get_branch_list(
		df_rgroup,
		df_recipe,
		combi,
	):

	rgroup_list = []
	for i in range(len(df_recipe)):
		row_ = df_recipe.iloc[i]
		class_ = row_['Class']

		rgroups = []
		if class_.lower() == 'all':
			rgroups = list(df_rgroup['Motif'])
		else:
			tmp_list = class_.split(':')
			class_list = [item.lower() for item in tmp_list]
			rgroups = list(df_rgroup[df_rgroup['Class'].isin(class_list)]['Motif'])
		rgroup_list.append(rgroups)
	
	branch_list = []
	if combi:
		all_combi = list(product(*rgroup_list))
		for item in all_combi:
			dict_ = {}
			for idx, rgroup in enumerate(item):
				dict_[idx+1] = Chem.MolFromSmiles(rgroup)
			branch_list.append(dict_)
	else:
		len_ = len(rgroup_list)
		for idx, rgroups in enumerate(rgroup_list):
			tmp_list = []
			for i in range(len_):
				tmp_list.append(['[*][H]'])
			tmp_list[idx] = rgroups
			all_combi = list(product(*tmp_list))
			for item in all_combi:
				dict_ = {}
				for idx, rgroup in enumerate(item):
					dict_[idx+1] = Chem.MolFromSmiles(rgroup)
				branch_list.append(dict_)

	return branch_list
