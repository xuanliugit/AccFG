# accfg.main

Core functional group detection logic, including the `AccFG` class and helper utilities for loading SMARTS patterns.

## `AccFG` Class

### Constructor

`AccFG(common_fgs=True, heterocycle_fgs=True, user_defined_fgs=None, print_load_info=False, lite=False)`

- **Parameters**
  - `common_fgs` (`bool`): Load bundled common functional groups from `accfg/fgs_common.csv`.
  - `heterocycle_fgs` (`bool`): Load heterocycle definitions from `accfg/fgs_heterocycle.csv`.
  - `user_defined_fgs` (`dict[str, str] | None`): Optional SMARTS overrides/additions supplied by the caller.
  - `print_load_info` (`bool`): Emit a summary of how many groups were loaded.
  - `lite` (`bool`): When true, the CSV reader skips rows marked with `%` or `#`, yielding a lighter rule set.
- **Behaviour**
  - Populates dictionaries for each source and merges them into `self.dict_fgs`.
  - User-defined patterns are canonicalised and `[nH]` is normalised to `[n]` via `process_user_defined_fgs`.
  - Stores whether the lite mode was requested for later use.

### `_is_fg_in_mol(self, mol, fg)`

Internal helper that returns `(is_present, matches)` for a SMARTS query on an RDKit molecule.

- `is_present` (`bool`): Whether any substructure matches were found.
- `matches` (`tuple[tuple[int, ...], ...]`): Atom index tuples for each match (RDKit order).

### `_freq_fg_in_mol(self, mol, fg)`

Count the occurrences of a functional group SMARTS in a molecule.

- Returns the integer `count` or `False` when no match is found.
- Useful when tallying frequency rather than raw match coordinates.

### `_get_bonds_from_match(self, query_mol, mol, atom_match)`

Translate an atom match from the query SMARTS into the corresponding bond indices on the target molecule. Used for pruning nested functional groups.

- **Parameters**: RDKit molecules for the query and target, plus a list/tuple of atom indices returned by RDKit.
- **Returns**: `list[int]` of bond indices in `mol`.

### `process_user_defined_fgs(self, user_defined_fgs)`

Sanitise custom functional group definitions.

- If a SMARTS string already uses `;` (indicating an explicit SMARTS pattern), it is preserved.
- Otherwise, the SMILES is canonicalised and aromatic `[nH]` handling is normalised to `[n]`.
- Returns a new dictionary with processed SMARTS strings.

### `run(self, smiles, show_atoms=True, show_graph=False, canonical=True)`

High-level entry point that accepts a SMILES string, optionally canonicalises it, converts it into an RDKit molecule, and defers to `run_mol`.

- **Parameters**
  - `smiles` (`str`): Molecule SMILES.
  - `show_atoms` (`bool`): When true, include atom index matches in the output.
  - `show_graph` (`bool`): When true, also include the functional group hierarchy graph.
  - `canonical` (`bool`): Canonicalise SMILES before processing.
- **Returns**
  - If `show_atoms` and not `show_graph`: `dict[str, list[tuple[int, ...]]]`.
  - If `show_atoms` and `show_graph`: tuple `(fg_dict, networkx.DiGraph)`.
  - Otherwise: `list[str]` of detected functional group names.
- **Usage**

```python
from accfg.main import AccFG

afg = AccFG(print_load_info=True)
fgs = afg.run("CC(=O)O")
```

### `run_mol(self, mol, show_atoms=True, show_graph=False, use_atom_map_num=False)`

Operate directly on an RDKit molecule and locate all functional groups using a `ProcessPoolExecutor`.

- **Parameters**
  - `mol` (`rdkit.Chem.Mol`): Pre-built molecule.
  - `show_atoms`, `show_graph`: Same semantics as `run`.
  - `use_atom_map_num` (`bool`): Remap atom indices to the `molAtomMapNumber` property when present (useful for aligned molecules).
- **Behaviour**
  - Launches up to four parallel worker processes to evaluate SMARTS patterns.
  - Builds a directed graph of parent/child relationships to remove redundant nested matches (e.g., derivatives).
- **Returns**: Same as `run`.
- **Notes**
  - When `use_atom_map_num=True`, the resulting atom lists are remapped from indices to mapping numbers.
  - For `show_graph=True`, the returned `DiGraph` stores `mapped_atoms` on each node for downstream visualisation.

### `run_freq(self, smiles)`

Return per-functional-group occurrence counts for a SMILES string.

- **Returns**: `list[tuple[str, int]]` on success or `None` if detection fails.
- **Usage**

```python
counts = afg.run_freq("c1ccccc1O")
```

### `csv_to_dict(self, csv_file, lite=False)`

Load functional group SMARTS definitions from a CSV file.

- Uses `csv.DictReader`, filtering comment lines. In lite mode, both `%` and `#` prefixed rows are skipped.
- Returns `dict[str, str]` mapping the functional group name to its SMARTS pattern.

## Internal Utilities

The following helpers support post-processing inside `run_mol`:

- `_get_bonds_from_match` and `_freq_fg_in_mol` are used to compare structural context between candidate groups.
- The functional group directed graph leverages `networkx` to reason about parent/child relationships before emitting results.

## Module-Level Helpers

### `canonical_smiles(smi)`

Convert any SMILES string into its RDKit canonical form using `Chem.MolToSmiles(Chem.MolFromSmiles(smi))`. Useful for normalising inputs before comparison.