from rdkit import Chem
from rdkit.Chem import rdRascalMCES
from collections import Counter
import warnings
from rdkit.Chem import rdFMCS
from collections import defaultdict
import warnings
from .main import AccFG, canonical_smiles

afg = AccFG()

def get_RascalMCES(smiles1, smiles2, similarityThreshold=0.7):
    if isinstance(smiles1, str) and isinstance(smiles2, str):
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
    else:
        mol1 = smiles1
        mol2 = smiles2
    opts = rdRascalMCES.RascalOptions()
    opts.ringMatchesRingOnly = True#True
    # opts.completeRingsOnly = True
    opts.ignoreAtomAromaticity = False#False#True #
    opts.maxBondMatchPairs = 1000
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

        # assert len(unique_atom_list) == number
        if len(unique_atom_list) != number:
            raise ValueError(f'Error on {unique_fgs} and {mapped_atoms}')
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
            alkane_list.append((f'C{number} alkane',tuple(atom_idx_list)))
        elif len(atom_count)==0:
            return []
        else:
            raise ValueError(f'Error on {Chem.MolToSmiles(remain_mol_alkane)}')         
    alkane_list_dict = dict()
    for alkane, atom_num_list in alkane_list:
        alkane_list_dict.setdefault(alkane, []).append(atom_num_list)
    alkane_list_with_len = [(alkane, len(atom_list), atom_list) for alkane, atom_list in alkane_list_dict.items()]
    return alkane_list_with_len 

def project_atom_num_to_atom_note(mol):
    for atom in mol.GetAtoms():
        atom.SetProp('atomNote', str(atom.GetProp('molAtomMapNumber')))
    return mol

def get_alkane_and_atom_from_remain_mol_with_remap(remain_mol_alkane, original_mol_with_atom_num):
    if remain_mol_alkane.GetNumAtoms() == 0:
        return []
    
    # 1. Get atomNote list from remain_mol_alkane
    remain_mol_alkane_atom_note_list = []
    for atom in remain_mol_alkane.GetAtoms():
        remain_mol_alkane_atom_note_list.append(int(atom.GetProp('atomNote')))
    # 2. Get fragment smiles from original_mol and remain_mol_alkane_atom_note_list
    
    fragment = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(original_mol_with_atom_num, remain_mol_alkane_atom_note_list))
    fragment_atom_note = project_atom_num_to_atom_note(fragment)
    alkane_frags = Chem.GetMolFrags(fragment_atom_note)

    alkane_list = []
    for alkane_frag in alkane_frags:
        atom_list = []
        atom_idx_list = []
        for atom_index in alkane_frag:
            atom = fragment.GetAtomWithIdx(atom_index)
            atom_list.append(atom.GetSymbol())
            atom_idx_list.append(int(atom.GetProp('atomNote')))
        atom_count = Counter(atom_list)
        if 'C' in atom_count and len(atom_count)==1:
            number = atom_count['C']
            alkane_list.append((f'C{number} alkane',tuple(atom_idx_list)))
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
    # MCS method
    Split the remaining molecules into smaller fragments and compare them with the reference remaining molecules.
    '''
    target_remain_alkane = []
    ref_remain_alkane = []
    
    for i in range(len(target_remain_mol_frags)):
        target_remain_mol_frag = target_remain_mol_frags[i]
        ref_remain_mol_frag = ref_remain_mol_frags[i]
        res = rdFMCS.FindMCS([target_remain_mol_frag, ref_remain_mol_frag])
        #res = rdFMCS.FindMCS([target_remain_mol_frag, ref_remain_mol_frag],ringMatchesRingOnly=True,completeRingsOnly=True)
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
def get_alkane_diff_from_mol_MCS(target_remain_mol, ref_remain_mol):
    # MCS method directly on the remaining molecules if target/ref only has one fragments
    target_remain_alkane = []
    ref_remain_alkane = []
    
    res = rdFMCS.FindMCS([target_remain_mol, ref_remain_mol])
    mcs_smarts = res.smartsString
    mcs_mol = Chem.MolFromSmarts(res.smartsString)

    target_remain_mol_match_atoms = target_remain_mol.GetSubstructMatch(mcs_mol)
    ref_remain_mol_match_atoms = ref_remain_mol.GetSubstructMatch(mcs_mol)
    target_remain_mol_frag_match_atoms = remove_atoms_from_mol(target_remain_mol, set(target_remain_mol_match_atoms))
    ref_remain_mol_frag_match_atoms = remove_atoms_from_mol(ref_remain_mol, set(ref_remain_mol_match_atoms))
    target_remain_frag_alkane = get_alkane_and_atom_from_remain_mol(target_remain_mol_frag_match_atoms)
    ref_remain_frag_alkane = get_alkane_and_atom_from_remain_mol(ref_remain_mol_frag_match_atoms)
    if target_remain_frag_alkane != []:
        target_remain_alkane.extend(target_remain_frag_alkane)
    if ref_remain_frag_alkane != []:
        ref_remain_alkane.extend(ref_remain_frag_alkane)
    return merge_alkane_synonyms(target_remain_alkane), merge_alkane_synonyms(ref_remain_alkane)



def get_alkane_diff_legacy(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms):
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
    mces_result_on_remain = get_RascalMCES(target_remain_mol, ref_remain_mol, similarityThreshold=0.1)
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
def get_atoms_from_diff(diff_tuple):
    atoms_list = []
    if diff_tuple == ([],[]): return []
    for diff_list in diff_tuple:
        if diff_list:
            for fg_list in diff_list:
                atoms = [atom for atom_list in fg_list[2] for atom in atom_list]
                atoms_list.extend(atoms)
    return list(set(atoms))

def get_atoms_list_from_diff(diff_tuple):
    atoms_list = []
    if diff_tuple == ([],[]): return []
    for diff_list in diff_tuple:
        if diff_list:
            for fg_list in diff_list:
                atoms = fg_list[2]
                atoms_list.extend(atoms)
    return atoms_list

def get_atoms_list_from_fg_list(fg_list):
    atoms_list = []
    if fg_list == []: return []
    for fgs in fg_list:
        atoms = fgs[2]
        atoms_list.extend(atoms)
    return atoms_list

def get_outer_bond_from_atoms(mol, atoms):
    '''
    fg_list: [(fg, number, atoms_list),...]
    outer_bond_list: [(atom1_idx, atom2_idx),...]
    '''
    outer_bond_list = []
    for atom_idx in atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in atoms:
                outer_bond_list.append((atom.GetIdx(), nbr.GetIdx()))
    return outer_bond_list

def get_outer_atoms_from_atoms(mol, atoms):
    '''
    fg_list: [(fg, number, atoms_list),...]
    outer_atoms_list: [atom_idx,...]
    '''
    outer_atoms_list = []
    for atom_idx in atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in atoms:
                outer_atoms_list.append(nbr.GetIdx())
    return outer_atoms_list
 
def get_atom_idx_from_atom_note(mol, atom_note_list):
    atom_note_idx_dict = {}
    for atom in mol.GetAtoms():
        if int(atom.GetProp('atomNote')) in atom_note_list:
            atom_note_idx_dict[int(atom.GetProp('atomNote'))] = atom.GetIdx()
    return atom_note_idx_dict

def add_hs_from_idx(mol, idx_list):
    for idx in idx_list:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs()+1)
    return mol

def remove_atoms_build_bond(mol, atom_idx, outer_atoms):
    ed_mol = Chem.RWMol(mol)
    ed_mol.BeginBatchEdit()
    for idx in atom_idx:
        ed_mol.RemoveAtom(idx)
    ed_mol.AddBond(outer_atoms[0], outer_atoms[1], order=Chem.rdchem.BondType.SINGLE)
    ed_mol.CommitBatchEdit()
    return ed_mol.GetMol()

def remove_atoms_add_hs(mol, atom_idx, outer_atoms):
    ed_mol = Chem.RWMol(mol)
    ed_mol.BeginBatchEdit()
    for idx in atom_idx:
        ed_mol.RemoveAtom(idx)
    for idx in outer_atoms:
        atom = ed_mol.GetAtomWithIdx(idx)
        
        atom.SetNumExplicitHs(atom.GetNumExplicitHs()+1)
        
    ed_mol.CommitBatchEdit()
    return ed_mol.GetMol()

def remove_fg_list_from_mol(mol, fg_list):
    '''
    fg_list: [(fg, number, atoms_list),...]
    '''
    atoms_list_to_remove = get_atoms_list_from_fg_list(fg_list)

    for atoms in atoms_list_to_remove:
        atom_note_idx_dict = get_atom_idx_from_atom_note(mol, atoms)
        atom_idx = list(atom_note_idx_dict.values())
        outer_atoms = get_outer_atoms_from_atoms(mol, atom_idx)
        outer_atoms_unique = list(set(outer_atoms))

        if len(outer_atoms_unique) == 2 and not mol.GetBondBetweenAtoms(outer_atoms_unique[0], outer_atoms_unique[1]):
            # Two outer atoms are not bonded
            mol = remove_atoms_build_bond(mol, atom_idx, outer_atoms_unique)
        else:
            mol = remove_atoms_add_hs(mol, atom_idx, outer_atoms)

    return mol

def get_outer_bond_from_fg_list(mol, fg_list):
    '''
    fg_list: [(fg, number, atoms_list),...]
    fg_list_outer_bonds: [(fg, number, atoms_list, outer_bonds_list),...]
    '''
    fg_list_outer_bonds = []
    for fg, number, atoms_list in fg_list:
        outer_bonds_list = []
        for atoms in atoms_list:
            outer_bonds = get_outer_bond_from_atoms(mol, atoms)
            outer_bonds_list.append(outer_bonds)
        fg_list_outer_bonds.append((fg, number, atoms_list, outer_bonds_list))
    return fg_list_outer_bonds

def get_alkane_diff_MCES(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms):

    target_mol = Chem.MolFromSmiles(target_smiles)
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    target_mol = set_atom_idx(target_mol,'atomNote')
    ref_mol = set_atom_idx(ref_mol,'atomNote')

    target_remain_mol = remove_fg_list_from_mol(target_mol, unique_target_fgs_atoms)
    ref_remain_mol = remove_fg_list_from_mol(ref_mol, unique_ref_fgs_atoms)

    Chem.SanitizeMol(target_remain_mol)
    Chem.SanitizeMol(ref_remain_mol)
    
    mces_result_on_remain = get_RascalMCES(target_remain_mol, ref_remain_mol, similarityThreshold=0.1)
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

def get_alkane_diff(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms):
    
    target_mol = Chem.MolFromSmiles(target_smiles)
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    target_mol = set_atom_idx(target_mol,'atomNote')
    ref_mol = set_atom_idx(ref_mol,'atomNote')
    
    target_mol_with_atom_num = Chem.MolFromSmiles(target_smiles)
    target_mol_with_atom_num = set_atom_idx(target_mol_with_atom_num,'molAtomMapNumber')
    ref_mol_with_atom_num = Chem.MolFromSmiles(ref_smiles)
    ref_mol_with_atom_num = set_atom_idx(ref_mol_with_atom_num,'molAtomMapNumber')

    target_remain_mol = remove_fg_list_from_mol(target_mol, unique_target_fgs_atoms)
    ref_remain_mol = remove_fg_list_from_mol(ref_mol, unique_ref_fgs_atoms)

    Chem.SanitizeMol(target_remain_mol)
    Chem.SanitizeMol(ref_remain_mol)

    target_remain_mol_frags = Chem.GetMolFrags(target_remain_mol, asMols=True)
    ref_remain_mol_frags = Chem.GetMolFrags(ref_remain_mol, asMols=True)
    # Edit
    if len(target_remain_mol_frags) == 1 and len(ref_remain_mol_frags) == 1:
        mces_result_on_remain = get_RascalMCES(target_remain_mol, ref_remain_mol, similarityThreshold=0.1)
        if len(mces_result_on_remain) != 0:
            target_mapped_atoms = [atom_pair[0] for atom_pair in mces_result_on_remain[0].atomMatches()]
            ref_mapped_atoms = [atom_pair[1] for atom_pair in mces_result_on_remain[0].atomMatches()]

            target_remain_mol_alkane = remove_atoms_from_mol(target_remain_mol, set(target_mapped_atoms))
            ref_remain_mol_alkane = remove_atoms_from_mol(ref_remain_mol, set(ref_mapped_atoms))
            
            target_remain_alkane = get_alkane_and_atom_from_remain_mol_with_remap(target_remain_mol_alkane, target_mol_with_atom_num)
            
            ref_remain_alkane = get_alkane_and_atom_from_remain_mol_with_remap(ref_remain_mol_alkane, ref_mol_with_atom_num)
            return target_remain_alkane, ref_remain_alkane
        
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
        if len(unique_target_fgs_atoms) == 0 and len(unique_ref_fgs_atoms) == 0 and len(mces_result) == 0:
            warnings.warn(f'Check the MCES setting! Error on {target_smiles} and {ref_smiles}.')
            return None
        elif len(unique_target_fgs_atoms) == 0 and len(unique_ref_fgs_atoms) == 0 and len(mces_result) != 0:
            target_remain_alkane, ref_remain_alkane = get_alkane_diff_MCES(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms)
        else:
            target_remain_alkane, ref_remain_alkane = get_alkane_diff(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms)
    except:
        try:
            warnings.warn("Using loose method to get the remaining alkanes.")
            target_remain_alkane, ref_remain_alkane = get_alkane_diff_loose(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms, target_mapped_atoms, ref_mapped_atoms)
        except:
            warnings.warn("Cannot get the remaining alkanes.")
            return None
    return (unique_target_fgs_atoms, target_remain_alkane), (unique_ref_fgs_atoms, ref_remain_alkane)
