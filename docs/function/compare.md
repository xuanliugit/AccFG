# accfg.compare

Comparison utilities built on top of AccFG that highlight structural differences between two molecules. Includes wrappers around RDKit’s Rascal MCES implementation plus extensive helpers for handling alkane residues.

## High-Level Workflow

### `compare_mols(target_smiles, ref_smiles, afg=AccFG(), similarityThreshold=0.7, canonical=True)`

End-to-end comparison between a target and reference molecule.

- **Parameters**
  - `target_smiles`, `ref_smiles` (`str`): Input SMILES strings.
  - `afg` (`accfg.main.AccFG`): Detector instance used to extract functional groups.
  - `similarityThreshold` (`float`): Passed to `get_RascalMCES` to control MCES sensitivity.
  - `canonical` (`bool`): Canonicalise both SMILES before further processing.
- **Returns**: Tuple `((unique_target_fgs_atoms, target_remain_alkane), (unique_ref_fgs_atoms, ref_remain_alkane))`.
  - `unique_*_fgs_atoms`: `list[tuple[str, int, list[tuple[int, ...]]]]` highlighting functional groups unique to each molecule.
  - `*_remain_alkane`: Similar list describing leftover alkane fragments.
- **Behaviour**
  1. Compute the Rascal MCES between the two molecules.
  2. Run AccFG on both molecules and extract functional group matches.
  3. Derive groups unique to each side (`get_unique_fgs_with_all_atoms` and `process_unique_fgs_atoms`).
  4. Analyse remaining structure to classify unshared alkane fragments via `get_alkane_diff` or fallback strategies.
- **Warnings**: Emits warnings when similarity is low, MCES results are empty, or when falling back to looser heuristics.
- **Usage**

```python
from accfg.compare import compare_mols

target, ref = compare_mols("CCO", "CC(=O)O")
target_unique, target_alkanes = target
```

## MCES Utilities

### `get_RascalMCES(smiles1, smiles2, similarityThreshold=0.7)`

Wrapper around `rdRascalMCES.FindMCES`. Accepts SMILES strings or RDKit molecules.

- **Parameters**: Two molecules plus an optional similarity threshold (default `0.7`).
- **Returns**: Rascal MCES result object (list of match records) from RDKit.

### `remove_atoms_from_mol(mol, atom_set)`

Produce a copy of a molecule with a set of atom indices removed using `Chem.RWMol`.

- **Parameters**: `mol` (`rdkit.Chem.Mol`), `atom_set` (`Iterable[int]`).
- **Returns**: New RDKit molecule with the selected atoms deleted.

### `get_unique_fgs_with_all_atoms(target_fgs, ref_fgs)`

Compute functional groups present in one molecule but not the other and track the raw atom lists.

- **Parameters**: `target_fgs`, `ref_fgs` (dicts as returned by `AccFG.run`).
- **Returns**: `(unique_target, unique_ref)` where each entry is a list of `(fg_name, count, atom_lists)`.
- Handles stringified dict input by calling `eval` (legacy compatibility).

### `process_unique_fgs_atoms(unique_fgs, mapped_atoms)`

Filter atom lists so only atoms outside the MCES-mapped region remain.

- **Parameters**
  - `unique_fgs`: Output from `get_unique_fgs_with_all_atoms`.
  - `mapped_atoms` (`list[int]`): Atoms covered by the MCES match.
- **Returns**: List of `(fg_name, count, filtered_atom_lists)`. Raises `ValueError` if counts disagree after filtering.

### `flatten_fg_diff_atoms(fg_diff_atoms)`

Flatten a list of functional group atom lists into a simple one-dimensional list of atom indices.

## Alkane Fragment Helpers

### `get_alkane_and_atom_from_remain_mol(remain_mol_alkane)`

Analyse fragments in a residual molecule to identify pure carbon chains (`Cn alkane`) based on `atomNote` properties.

- **Parameters**: Residual molecule with `atomNote` assigned to each atom.
- **Returns**: List of tuples `(name, count, atom_idx_tuples)` where `name` looks like `C3 alkane`.

### `project_atom_num_to_atom_note(mol)`

Copy each atom’s `molAtomMapNumber` property into `atomNote` and return the mutated molecule.

### `get_alkane_and_atom_from_remain_mol_with_remap(remain_mol_alkane, original_mol_with_atom_num)`

Variant of `get_alkane_and_atom_from_remain_mol` that reprojects fragment atoms back onto the original molecule using atom notes from `original_mol_with_atom_num`.

### `set_atom_idx(smi, label='molAtomMapNumber')`

Tag atoms with their index under the provided property key. Accepts SMILES or RDKit molecule.

### `merge_alkane_synonyms(fg_list)`

Combine multiple entries referring to the same alkane label by concatenating their atom lists.

### `get_alkane_diff_split(target_remain_mol_frags, ref_remain_mol_frags)`

Compare multiple residual fragments by iterating over fragment pairs and finding their MCS.

- **Parameters**: Two sequences of RDKit fragment molecules (as produced by `Chem.GetMolFrags(..., asMols=True)`).
- **Returns**: Tuple `(target_alkanes, ref_alkanes)` after merging synonyms.

### `get_alkane_diff_from_mol_MCS(target_remain_mol, ref_remain_mol)`

Convenience wrapper when both residual molecules contain a single fragment. Performs MCES once and classifies the leftover atoms as alkanes.

### `get_alkane_diff_legacy(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms)`

Legacy approach that removes unique functional group atoms and then searches for remaining alkanes via MCES, with fragment splitting fallback.

### `get_alkane_diff_MCES(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms)`

Modern approach that removes unique functional groups (after converting to RDKit molecules with atom mapping notes) and runs MCES on the remaining structures. Falls back to fragment splitting if no MCES is found.

### `get_alkane_diff(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms)`

Main entry point for computing alkane differences. Handles remapping of atom notes and chooses between direct MCES, fragment splitting, or fallback heuristics.

### `get_alkane_diff_loose(target_smiles, unique_target_fgs_atoms, ref_smiles, unique_ref_fgs_atoms, target_mapped_atoms, ref_mapped_atoms)`

Fallback used when standard MCES approaches fail. Removes both unique FG atoms and mapped atoms, then labels the residue as alkanes directly.

## Atom Index & Bond Utilities

These helpers manipulate atom index collections during FG comparisons.

- `get_atoms_from_diff(diff_tuple)`: Flatten the atom indices from a pair of FG diff lists, returning unique indices.
- `get_atoms_list_from_diff(diff_tuple)`: Similar to above but preserves grouping.
- `get_atoms_list_from_fg_list(fg_list)`: Convenience for extracting atom lists from `(fg, count, atoms)` structures.
- `get_outer_bond_from_atoms(mol, atoms)`: Enumerate bonds that connect a given atom set to external atoms.
- `get_outer_atoms_from_atoms(mol, atoms)`: Return indices of neighbouring atoms outside the provided set.
- `get_atom_idx_from_atom_note(mol, atom_note_list)`: Map `atomNote` values back to RDKit indices.
- `add_hs_from_idx(mol, idx_list)`: Increment explicit hydrogen counts for specific atoms.
- `remove_atoms_build_bond(mol, atom_idx, outer_atoms)`: Remove a set of atoms and connect two outer atoms with a new single bond.
- `remove_atoms_add_hs(mol, atom_idx, outer_atoms)`: Remove atoms and compensate neighbours by adding hydrogens.
- `remove_atoms_list_from_mol(mol, atoms_list)`: Iteratively apply the above removal routines for multiple atom groups.
- `remove_fg_list_from_mol(mol, fg_list)`: Remove all atoms associated with functional groups from a molecule by delegating to `remove_atoms_list_from_mol`.
- `get_outer_bond_from_fg_list(mol, fg_list)`: Annotate functional group entries with their associated outer bonds.

## Legacy & Support Constants

- `afg = AccFG()`: Module-level default detector used when callers omit an explicit instance. Reuse it for convenience but create your own if you need custom functional group definitions.
