from rdkit import Chem
from rdkit.Chem import Draw
from accfg import (AccFG, draw_mol_with_fgs, molimg, 
                   img_grid,  compare_mols, draw_compare_mols, set_atom_idx, canonical_smiles,
                   draw_RascalMCES, print_fg_tree)

from IPython.display import Image
import networkx as nx
import argparse

afg = AccFG(print_load_info=True)

def run(smi, show_atoms, show_graph):
    smi = canonical_smiles(smi)
    print(f'Input SMILES (canonical): {Chem.MolToSmiles(set_atom_idx(smi))}')
    if show_graph:
        fgs,fg_graph = afg.run(smi, show_atoms=show_atoms, show_graph=show_graph)
        print_fg_tree(fg_graph, fgs.keys(), show_atom_idx=True)
    else:
        fgs = afg.run(smi, show_atoms=show_atoms, show_graph=show_graph)
        print(fgs)
        
def run_compare(smi1, smi2, similarityThreshold):
    smi1 = canonical_smiles(smi1)
    smi2 = canonical_smiles(smi2)
    mol1_fg_alkane_df, mol2_fg_alkane_df = compare_mols(smi1, smi2, afg=afg, similarityThreshold=similarityThreshold)
    print(f'***Functional group comparison***')
    print(f'Similarity threshold: {similarityThreshold}')
    print('---'*10)
    print(f'Target molecule: {Chem.MolToSmiles(set_atom_idx(smi1))}')
    print(f'Unique functional groups in target: {mol1_fg_alkane_df[0]}')
    print(f'Unique alkane in target: {mol1_fg_alkane_df[1]}')
    print('---'*10)
    print(f'Reference molecule: {Chem.MolToSmiles(set_atom_idx(smi2))}')
    print(f'Unique functional groups in reference: {mol2_fg_alkane_df[0]}')
    print(f'Unique alkane in reference: {mol2_fg_alkane_df[1]}')
    

def parse_args():
    parser = argparse.ArgumentParser(description='Run AccFG on a SMILES string')
    parser.add_argument('smi', type=str, help='SMILES string')
    parser.add_argument('--show_atoms', default=True, help='Show the atoms in the functional groups')
    parser.add_argument('--show_graph', default=True, help='Show the functional group graph')
    parser.add_argument('--compare_smi', type=str, default=None, help='The SMILES strings to be compared')
    parser.add_argument('--similarityThreshold', type=float, default=0.5, help='The similarity threshold for comparing')
    return parser

if __name__ == '__main__':
    '''
    # Get functional groups from SMILES
    python run_accfg.py 'CN(C)/N=N/C1=C(NC=N1)C(=O)N'

    # Compare two molecules
    python run_accfg.py 'CNC(=O)Cc1nc(-c2ccccc2)cs1' --compare_smi 'CCNCCc1nc2ccccc2s1'
    '''
    parser = parse_args()
    args = parser.parse_args()
    if not args.compare_smi:
        run(args.smi, args.show_atoms, args.show_graph)
    else:
        smi1, smi2 = args.smi, args.compare_smi
        run_compare(smi1, smi2, args.similarityThreshold)