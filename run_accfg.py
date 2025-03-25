from rdkit import Chem
from rdkit.Chem import Draw
from accfg import (AccFG, draw_mol_with_fgs, molimg, 
                   img_grid,  compare_mols, draw_compare_mols,
                   draw_RascalMCES, print_fg_tree)

from IPython.display import Image
import networkx as nx
import argparse

afg = AccFG(print_load_info=True)

def run(smi, show_atoms, show_graph):
    if show_graph:
        fgs,fg_graph = afg.run(smi, show_atoms=show_atoms, show_graph=show_graph)
        print_fg_tree(fg_graph, fgs.keys(), show_atom_idx=True)
    else:
        fgs = afg.run(smi, show_atoms=show_atoms, show_graph=show_graph)
        print(fgs)

def parse_args():
    parser = argparse.ArgumentParser(description='Run AccFG on a SMILES string')
    parser.add_argument('smi', type=str, help='SMILES string')
    parser.add_argument('--show_atoms', default=True, help='Show the atoms in the functional groups')
    parser.add_argument('--show_graph', default=True, help='Show the functional group graph')
    return parser

if __name__ == '__main__':
    parser = parse_args()
    args = parser.parse_args()
    run(args.smi, args.show_atoms, args.show_graph)