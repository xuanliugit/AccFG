from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdRascalMCES
from rdkit.Chem import Draw
from rdkit import DataStructs
import heapq
import re
import csv
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np
import os
from collections import Counter
from tqdm import tqdm
from IPython.display import display
import swifter
from IPython.display import SVG
import networkx as nx
import warnings
from rdkit.Chem import rdFMCS
from collections import defaultdict

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class AccFG():
    def __init__(self, common_fgs=True, heterocycle_fgs=True, user_defined_fgs={}, print_load_info=False):
        log_text = ""
        if common_fgs:
            self.dict_fgs_common_path = os.path.join(PROJECT_DIR, 'accfg/fgs_common.csv')
            
            self.dict_fgs_common = self.csv_to_dict(self.dict_fgs_common_path)
            log_text += f"Loaded {len(self.dict_fgs_common)} common functional groups. "
        else:
            self.dict_fgs_common = {}
        if heterocycle_fgs:
            self.dict_fg_heterocycle_path =  os.path.join(PROJECT_DIR,'accfg/fgs_heterocycle.csv')
            
            self.dict_fg_heterocycle = self.csv_to_dict(self.dict_fg_heterocycle_path)
            log_text += f"Loaded {len(self.dict_fg_heterocycle)} heterocycle groups. "
        else:
            self.dict_fg_heterocycle = {}
        if user_defined_fgs:
            self.dict_fgs_user_defined = self.process_user_defined_fgs(user_defined_fgs)
            log_text += f"Loaded {len(user_defined_fgs)} user-defined functional groups. "
        else:
            self.dict_fgs_user_defined = {}
        self.dict_fgs = {**self.dict_fgs_common, **self.dict_fg_heterocycle, **self.dict_fgs_user_defined}
        
        if print_load_info:
            print(f'{log_text}Total {len(self.dict_fgs)} functional groups loaded.')
        
    def _is_fg_in_mol(self, mol, fg):
        fgmol = Chem.MolFromSmarts(fg)
        mol = Chem.MolFromSmiles(mol.strip())
        mapped_atoms = Chem.Mol.GetSubstructMatches(mol, fgmol, uniquify=True)
        if_mapped = len(mapped_atoms) > 0
        return if_mapped, mapped_atoms
    
    def _get_bonds_from_match(self, query_mol, mol, atom_match):
        """
        Args:
        query_mol: query molecule used to match
        mol: molecule matched
        atom_match: result of GetSubstructMatch (i.e. matched atom idx list)
        Returns:
        list of matched bond indices, or None
        """
        bonds = []
        if isinstance(atom_match, (list, tuple)):
            pass
        else:
            atom_match = list(atom_match)
        for bond in query_mol.GetBonds():
            idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonds.append(mol.GetBondBetweenAtoms(atom_match[idx1], atom_match[idx2]).GetIdx())
        return bonds
    
    def process_user_defined_fgs(self, user_defined_fgs):
        user_defined_fgs_edit = {}
        for fg_name, fg_smi in user_defined_fgs.items():
            fg_smi_edit = Chem.MolToSmiles(Chem.MolFromSmiles(fg_smi))
            fg_smi_edit = fg_smi_edit.replace('[nH]','[n]')
            user_defined_fgs_edit[fg_name] = fg_smi_edit
        return user_defined_fgs_edit
    
    def run(self, smiles: str, show_atoms=True, show_graph=False, canonical=True) -> dict:
        """
        Input a molecule SMILES or name.
        Returns a list of functional groups identified by their common name (in natural language).
        """
        if canonical:
            smiles = canonical_smiles(smiles)
        
        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = {
                executor.submit(self._is_fg_in_mol, smiles, fg): name
                for name, fg in self.dict_fgs.items()
            }
            fgs_in_molec = {futures[future]: future.result()[1] for future in futures if future.result()[0]}
            
            # Check if the functional groups are subgroups of other functional groups
            # Build FG graph
            fg_graph = nx.DiGraph()
            fg_graph.add_nodes_from(list(fgs_in_molec.keys()))
            
            for name, mapped_atoms in list(fgs_in_molec.items()):
                fg_graph.nodes[name]['mapped_atoms'] = mapped_atoms
                
                remained_mapped_atoms_tuple_list = list(mapped_atoms) # a list of tuples containing atoms for the functional group
                for ref_name, ref_mapped_atoms in list(fgs_in_molec.items()):
                    if name != ref_name:
                        for target_atoms in mapped_atoms: # check if the target atoms are a subset of the reference atoms
                            for ref_atoms in ref_mapped_atoms:
                                if (set(target_atoms) < set(ref_atoms)) and ('derivative' not in ref_name):#and (target_atoms in remained_mapped_atoms_tuple_list) 
                                    if target_atoms in remained_mapped_atoms_tuple_list: remained_mapped_atoms_tuple_list.remove(target_atoms)
                                    fg_graph.add_edge(ref_name, name)
                                    
                                elif (set(target_atoms) == set(ref_atoms)) and ('derivative' not in ref_name):#and (target_atoms in remained_mapped_atoms_tuple_list) 
                                    # If mapping the same set of atoms Check if the number of bonds is smaller than the reference
                                    mol = Chem.MolFromSmiles(smiles)
                                    query_mol_ref = Chem.MolFromSmarts(self.dict_fgs[ref_name])
                                    query_mol_target = Chem.MolFromSmarts(self.dict_fgs[name])
                                    
                                    ref_bonds = self._get_bonds_from_match(query_mol_ref, mol, ref_atoms)
                                    target_bonds = self._get_bonds_from_match(query_mol_target, mol, target_atoms)
                                    if len(target_bonds) < len(ref_bonds):
                                        if target_atoms in remained_mapped_atoms_tuple_list: remained_mapped_atoms_tuple_list.remove(target_atoms)
                                        fg_graph.add_edge(ref_name, name)
                                    if len(target_bonds) == len(ref_bonds): # only remove atoms
                                        if target_atoms in remained_mapped_atoms_tuple_list: remained_mapped_atoms_tuple_list.remove(target_atoms)
                        if len(remained_mapped_atoms_tuple_list) == 0:
                            fgs_in_molec.pop(name,None)
                            if show_graph:
                                continue
                            else:
                                break
                            # break
                if len(remained_mapped_atoms_tuple_list) > 0:
                    fgs_in_molec[name] = remained_mapped_atoms_tuple_list
        
        if show_atoms and not show_graph:
            return fgs_in_molec
        elif show_atoms and show_graph:
            return fgs_in_molec, fg_graph
        else:
            return list(fgs_in_molec.keys())
        #except:
        #    return None
            
    def run_freq(self, smiles: str) -> str:
        """
        Input a molecule SMILES or name.
        Returns a list of functional groups identified by their common name (in natural language).
        """
        try:
            fgs = self.run(smiles, show_atoms=True)
            fgs_freq = [(fg, len(mapped_atoms)) for fg, mapped_atoms in fgs.items()]
            return fgs_freq
        except:
            return None
            # return "Wrong argument. Please input a valid molecular SMILES."
    def csv_to_dict(self, csv_file):
        data = {}
        with open(csv_file, 'r') as file:
            #df = pd.read_csv(file, quotechar='"', comment='#')
            reader = csv.DictReader(filter(lambda row: row[0]!='#', file))
            for row in reader:
                key = row.pop('Functional Group')
                data[key] = row.pop('SMARTS Pattern')
        return data
    

# Comparison
def canonical_smiles(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

def get_RascalMCES(smiles1, smiles2, similarityThreshold=0.7):
    if isinstance(smiles1, str) and isinstance(smiles2, str):
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
    else:
        mol1 = smiles1
        mol2 = smiles2
    opts = rdRascalMCES.RascalOptions()
    opts.ringMatchesRingOnly = False#True
    # opts.completeRingsOnly = True
    opts.ignoreAtomAromaticity = True#True #
    if similarityThreshold:
        opts.similarityThreshold = similarityThreshold
        
    res = rdRascalMCES.FindMCES(mol1, mol2,opts)
    return res


def remove_atoms_from_mol(mol, atom_set):
    ed_mol = Chem.RWMol(mol)
    ed_mol.BeginBatchEdit()
    for atom in atom_set:
        ed_mol.RemoveAtom(atom)
    ed_mol.CommitBatchEdit()
    return ed_mol.GetMol()


def get_unique_fgs_with_all_atoms(target_fgs, ref_fgs):
    if isinstance(target_fgs, str):
        target_fgs = eval(target_fgs)
    if isinstance(ref_fgs, str):
        ref_fgs = eval(ref_fgs)
    unique_target_fgs = []
    for fg in target_fgs:
        if fg not in ref_fgs:
            unique_target_fgs.append((fg,len(target_fgs[fg]),target_fgs[fg]))
        elif fg in ref_fgs and len(target_fgs[fg]) > len(ref_fgs[fg]):
            unique_target_fgs.append((fg,len(target_fgs[fg])-len(ref_fgs[fg]),target_fgs[fg]))
    unique_ref_fgs = []
    for fg in ref_fgs:
        if fg not in target_fgs:
            unique_ref_fgs.append((fg,len(ref_fgs[fg]),ref_fgs[fg]))
        elif len(ref_fgs[fg]) > len(target_fgs[fg]):
            unique_ref_fgs.append((fg,len(ref_fgs[fg])-len(target_fgs[fg]),ref_fgs[fg]))
    return unique_target_fgs, unique_ref_fgs

def process_unique_fgs_atoms(unique_fgs, mapped_atoms):
    '''
    Only keep the atoms that are not in the mapped atoms
    '''
    unique_fgs_atoms = []
    for fg_name, number, atom_list in unique_fgs:
        unique_atom_list = []
        if number == len(atom_list):
            unique_fgs_atoms.append((fg_name, number, atom_list))
            continue
        for atom_set in atom_list:
            if set(atom_set).issubset(set(mapped_atoms)):
                continue
            unique_atom_list.append(atom_set)

        assert len(unique_atom_list) == number
        unique_fgs_atoms.append((fg_name, number, unique_atom_list))
    return unique_fgs_atoms

def flatten_fg_diff_atoms(fg_diff_atoms):
    return [atom for fgs in fg_diff_atoms for atoms in fgs for atom in atoms]

def get_alkane_and_atom_from_remain_mol(remain_mol_alkane):
    alkane_frags = Chem.GetMolFrags(remain_mol_alkane)
    alkane_list = []
    for alkane_frag in alkane_frags:
        atom_list = []
        atom_idx_list = []
        for atom_index in alkane_frag:
            atom = remain_mol_alkane.GetAtomWithIdx(atom_index)
            atom_list.append(atom.GetSymbol())
            atom_idx_list.append(int(atom.GetProp('atomNote')))
        atom_count = Counter(atom_list)
        if 'C' in atom_count and len(atom_count)==1:
            number = atom_count['C']
            alkane_list.append((f'C{number} alkane',atom_idx_list))
        elif len(atom_count)==0:
            return []
        else:
            raise ValueError(f'Error on {Chem.MolToSmiles(remain_mol_alkane)}')         
    alkane_list_dict = dict()
    for alkane, atom_num_list in alkane_list:
        alkane_list_dict.setdefault(alkane, []).append(atom_num_list)
    alkane_list_with_len = [(alkane, len(atom_list), atom_list) for alkane, atom_list in alkane_list_dict.items()]
    return alkane_list_with_len 
def set_atom_idx(smi, label = 'molAtomMapNumber'):
    #https://chemicbook.com/2021/03/01/how-to-show-atom-numbers-in-rdkit-molecule.html
    if isinstance(smi, str):
        mol  = Chem.MolFromSmiles(smi)
    else:
        mol = smi
    for atom in mol.GetAtoms():
        atom.SetProp(label,str(atom.GetIdx()))
    return mol

def merge_alkane_synonyms(fg_list):
    merged_dict = defaultdict(list)
    for fg_name, count, atom_list in fg_list:
        merged_dict[fg_name].extend(atom_list)
    merged_list = [(fg_name, len(atom_list), atom_list) for fg_name, atom_list in merged_dict.items()]
    return merged_list

def get_alkane_diff_split(target_remain_mol_frags, ref_remain_mol_frags):
    '''
    Split the remaining molecules into smaller fragments and compare them with the reference remaining molecules.
    '''
    target_remain_alkane = []
    ref_remain_alkane = []
     
    for i in range(len(target_remain_mol_frags)):
        target_remain_mol_frag = target_remain_mol_frags[i]
        ref_remain_mol_frag = ref_remain_mol_frags[i]
        res = rdFMCS.FindMCS([target_remain_mol_frag, ref_remain_mol_frag])
        mcs_smarts = res.smartsString
        mcs_mol = Chem.MolFromSmarts(res.smartsString)
        
        target_remain_mol_frag_match_atoms = target_remain_mol_frag.GetSubstructMatch(mcs_mol)
        ref_remain_mol_frag_match_atoms = ref_remain_mol_frag.GetSubstructMatch(mcs_mol)
        
        target_remain_mol_frag_match_atoms = remove_atoms_from_mol(target_remain_mol_frag, set(target_remain_mol_frag_match_atoms))
        ref_remain_mol_frag_match_atoms = remove_atoms_from_mol(ref_remain_mol_frag, set(ref_remain_mol_frag_match_atoms))
        target_remain_frag_alkane = get_alkane_and_atom_from_remain_mol(target_remain_mol_frag_match_atoms)
        ref_remain_frag_alkane = get_alkane_and_atom_from_remain_mol(ref_remain_mol_frag_match_atoms)
        if target_remain_frag_alkane != []:
            target_remain_alkane.extend(target_remain_frag_alkane)
        if ref_remain_frag_alkane != []:
            
            ref_remain_alkane.extend(ref_remain_frag_alkane)
            
    for i in range(len(target_remain_mol_frags), len(ref_remain_mol_frags)):
        ref_remain_mol_frag = ref_remain_mol_frags[i]
        ref_remain_frag_alkane = get_alkane_and_atom_from_remain_mol(ref_remain_mol_frag)
        if ref_remain_frag_alkane != []:
            ref_remain_alkane.extend(ref_remain_frag_alkane)
    return merge_alkane_synonyms(target_remain_alkane), merge_alkane_synonyms(ref_remain_alkane)

def get_alkane_diff(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms):
    target_fg_diff_atoms = [unique_atom_list for _,_,unique_atom_list in unique_target_fgs_atoms]
    ref_fg_diff_atoms = [unique_atom_list for _,_,unique_atom_list in unique_ref_fgs_atoms]

    target_fg_diff_atoms = flatten_fg_diff_atoms(target_fg_diff_atoms)
    ref_fg_diff_atoms = flatten_fg_diff_atoms(ref_fg_diff_atoms)
    
    target_mol = Chem.MolFromSmiles(target_smiles)
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    target_mol = set_atom_idx(target_mol,'atomNote')
    ref_mol = set_atom_idx(ref_mol,'atomNote')
    
    target_remain_mol = remove_atoms_from_mol(target_mol, set(target_fg_diff_atoms))
    ref_remain_mol = remove_atoms_from_mol(ref_mol, set(ref_fg_diff_atoms))
    Chem.SanitizeMol(target_remain_mol)
    Chem.SanitizeMol(ref_remain_mol)
    
    mces_result_on_remain = get_RascalMCES(target_remain_mol, ref_remain_mol, similarityThreshold=0.01)
    if len(mces_result_on_remain) != 0:
        target_mapped_atoms = [atom_pair[0] for atom_pair in mces_result_on_remain[0].atomMatches()]
        ref_mapped_atoms = [atom_pair[1] for atom_pair in mces_result_on_remain[0].atomMatches()]
        target_remain_mol_alkane = remove_atoms_from_mol(target_remain_mol, set(target_mapped_atoms))
        ref_remain_mol_alkane = remove_atoms_from_mol(ref_remain_mol, set(ref_mapped_atoms))
        
        target_remain_alkane = get_alkane_and_atom_from_remain_mol(target_remain_mol_alkane)
        ref_remain_alkane = get_alkane_and_atom_from_remain_mol(ref_remain_mol_alkane)
        return target_remain_alkane, ref_remain_alkane
    else: # If the MCES result is empty, try to split the remaining molecules into smaller fragments
        target_remain_mol_frags = Chem.GetMolFrags(target_remain_mol, asMols=True)
        ref_remain_mol_frags = Chem.GetMolFrags(ref_remain_mol, asMols=True)
        
        target_remain_alkane = []
        ref_remain_alkane = []
        
        if len(target_remain_mol_frags) <= len(ref_remain_mol_frags):
            target_remain_alkane,ref_remain_alkane = get_alkane_diff_split(target_remain_mol_frags, ref_remain_mol_frags)
            return target_remain_alkane, ref_remain_alkane
        else:
            ref_remain_alkane, target_remain_alkane = get_alkane_diff_split(ref_remain_mol_frags, target_remain_mol_frags)
            return target_remain_alkane, ref_remain_alkane
        
    
def get_alkane_diff_loose(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms, target_mapped_atoms, ref_mapped_atoms):
    '''
    Use this method when the MCES result is empty. This method is not as accurate as the get_alkane_diff method.
    '''
    target_fg_diff_atoms = [unique_atom_list for _,_,unique_atom_list in unique_target_fgs_atoms]
    ref_fg_diff_atoms = [unique_atom_list for _,_,unique_atom_list in unique_ref_fgs_atoms]

    target_fg_diff_atoms = flatten_fg_diff_atoms(target_fg_diff_atoms)
    ref_fg_diff_atoms = flatten_fg_diff_atoms(ref_fg_diff_atoms)
    
    target_mol = Chem.MolFromSmiles(target_smiles)
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    target_mol = set_atom_idx(target_mol,'atomNote')
    ref_mol = set_atom_idx(ref_mol,'atomNote')
    
    target_atom_to_remove = set(target_fg_diff_atoms) | set(target_mapped_atoms)
    ref_atom_to_remove = set(ref_fg_diff_atoms) | set(ref_mapped_atoms)
    
    target_remain_mol = remove_atoms_from_mol(target_mol, set(target_atom_to_remove))
    ref_remain_mol = remove_atoms_from_mol(ref_mol, set(ref_atom_to_remove))
    Chem.SanitizeMol(target_remain_mol)
    Chem.SanitizeMol(ref_remain_mol)
    
    target_remain_alkane = get_alkane_and_atom_from_remain_mol(target_remain_mol)
    ref_remain_alkane = get_alkane_and_atom_from_remain_mol(ref_remain_mol)
    return target_remain_alkane, ref_remain_alkane


def compare_mols(target_smiles, ref_smiles, afg = AccFG(), similarityThreshold=0.7, canonical=True):
    if canonical:
        target_smiles = canonical_smiles(target_smiles)
        ref_smiles = canonical_smiles(ref_smiles)
    
    mces_result = get_RascalMCES(target_smiles, ref_smiles, similarityThreshold)
    if len(mces_result) == 0:
        warnings.warn(f'target_smiles: {target_smiles} and ref_smiles: {ref_smiles} has low similarity. MCES result is empty. Try to lower the similarityThreshold.')
        target_mapped_atoms = []
        ref_mapped_atoms = []
    else:
        target_mapped_atoms = [atom_pair[0] for atom_pair in mces_result[0].atomMatches()]
        ref_mapped_atoms = [atom_pair[1] for atom_pair in mces_result[0].atomMatches()]

    target_fg = afg.run(target_smiles)
    ref_fg = afg.run(ref_smiles)
    
    unique_target_fgs, unique_ref_fgs = get_unique_fgs_with_all_atoms(target_fg, ref_fg)
    
    unique_target_fgs_atoms = process_unique_fgs_atoms(unique_target_fgs, target_mapped_atoms)
    unique_ref_fgs_atoms = process_unique_fgs_atoms(unique_ref_fgs, ref_mapped_atoms)
    try:
        target_remain_alkane, ref_remain_alkane = get_alkane_diff(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms)
    except:
        try:
            target_remain_alkane, ref_remain_alkane = get_alkane_diff_loose(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms, target_mapped_atoms, ref_mapped_atoms)
            warnings.warn("Using loose method to get the remaining alkanes.")
        except:
            warnings.warn("Cannot get the remaining alkanes.")
            return (unique_target_fgs_atoms, []), (unique_ref_fgs_atoms, [])
    return (unique_target_fgs_atoms, target_remain_alkane), (unique_ref_fgs_atoms, ref_remain_alkane)

