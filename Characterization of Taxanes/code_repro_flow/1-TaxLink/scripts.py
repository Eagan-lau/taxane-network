# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 08:10:33 2024

@author: DELL
"""


import os
import json
import requests
import subprocess
import numpy as np
import pandas as pd

from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, inchi


with open('save_folder/chemical_templetes/templates_general.json') as js:
    chemical_templetes = list(json.load(js).keys())
    chemical_rxns = [AllChem.ReactionFromSmarts(temp) for temp in chemical_templetes]

biological_templetes = pd.read_csv('save_folder/biological_templetes/BioTemplates.txt', sep=',')['template']
'''
with open('save_folder/biological_templetes/BioTemplates.txt') as txt:
    biological_templetes = txt.readlines()
'''
biological_rxns = [AllChem.ReactionFromSmarts(temp) for temp in biological_templetes]

    
def flatten_list(lst):
    flattened = []
    for item in lst:
        if isinstance(item, tuple):
            flattened.extend(flatten_list(item))
        else:
            flattened.append(item)
    return flattened


def predict_products_biological(smi):
    mol = Chem.MolFromSmiles(smi)
    products = [rxn.RunReactants([mol]) for rxn in biological_rxns]
    products = [s for s in products if len(s) > 0]
    products = flatten_list(products)
    products_smiles = np.unique([Chem.MolToSmiles(p) for p in products])
    return list(products_smiles)


tax_list = pd.read_excel('data/TaxList1.xlsx')
smi_list = tax_list['Isomeric SMILES']
mol_list = [Chem.MolFromSmiles(s) for s in smi_list]
mol_list = [m for m in mol_list if m is not None]
smi_list = [Chem.MolToSmiles(m) for m in mol_list]
short_keys = [inchi.MolToInchiKey(m)[:14] for m in mol_list if m is not None]

link_matrix = np.zeros((len(smi_list), len(smi_list)))
from_to_list = []

# biological transformation
for i, smi in enumerate(tqdm(smi_list)):
    prods = predict_products_biological(smi)
    prods_mols = [Chem.MolFromSmiles(m) for m in prods]
    prods_shortkeys = [inchi.MolToInchiKey(m)[:14] for m in prods_mols if m is not None]
    j = np.where([k in prods_shortkeys for k in short_keys])[0]
    if len(j) > 0:
        for jj in j:
            link_matrix[i,jj] = 1
            link_matrix[jj,i] = 1
            from_to_list.append([jj, i])

with open('smiles_list.txt', 'w') as txt:
    txt.writelines(smi_list)
pd.DataFrame(link_matrix).to_csv('linkage_new.csv')
pd.DataFrame(from_to_list).to_csv('from_to_list.csv')
