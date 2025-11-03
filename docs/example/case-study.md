# Case Study: Functional Group Analysis With AccFG

This case study walks through the workflows demonstrated in `example.ipynb`. Follow along to learn how to extract functional groups, visualise annotations, compare molecules, and extend AccFG with custom SMARTS definitions.

## 1. Environment Setup

```python
from rdkit import Chem
from rdkit.Chem import Draw

from accfg import (
    AccFG,
    draw_mol_with_fgs,
    molimg,
    img_grid,
    compare_mols,
    draw_compare_mols,
    draw_RascalMCES,
    print_fg_tree,
)

afg = AccFG(print_load_info=True)
```

The autoreload extension ensures edits to `accfg` are picked up automatically during interactive sessions. Instantiating `AccFG` loads the bundled functional group dictionaries and prints a summary when `print_load_info=True`.

## 2. Single-Molecule Functional Group Extraction

```python
smi = "CN(C)/N=N/C1=C(NC=N1)C(=O)N"
fgs, fg_graph = afg.run(smi, show_atoms=True, show_graph=True, canonical=False)

print("Functional groups:\n", fgs)
print("----------------")
print_fg_tree(fg_graph, fgs.keys(), show_atom_idx=True)
```

- `AccFG.run` returns both the mapping of functional group names to atom indices and, when `show_graph=True`, the dependency tree showing how groups nest.
- `print_fg_tree` renders the `networkx.DiGraph` as an ASCII hierarchy, optionally appending mapped atom indices. This is useful for inspecting derivative relationships (e.g., a carbonyl contained in an amide).

### Visualising FG Highlights

```python
mol_fg_img = molimg(
    draw_mol_with_fgs(
        smi,
        canonical=False,
        alpha=1,
    )
)
mol_fg_img
```

`draw_mol_with_fgs` runs AccFG internally, highlights each match with a distinct colour, and returns PNG bytes generated via RDKit’s Cairo drawer. Wrapping the bytes with `molimg` (a Pillow helper) gives you a display-ready image inside notebooks or scripts.

## 3. Batch Highlighting: Ertl Dataset

```python
ertl_smis = [
    "Cc1nc(NS(=O)(=O)c2ccc(N)cc2)nc(C)c1",
    "NC(=N)c1ccc(C=Cc2ccc(cc2O)C(=N)N)cc1",
    # ...
]

img_list = [molimg(draw_mol_with_fgs(smi, alpha=1)) for smi in ertl_smis]
mols_img = img_grid(img_list)
mols_img
```

- Iterate over the SMILES list, generating one highlighted image per molecule.
- `img_grid` composes the individual images into a tiled canvas (the notebook saves a high-resolution PNG in `results/result_on_ertl_mols.png` for reuse).
- Resizing the final grid (`mols_img.resize(...)`) produces presentation-friendly dimensions.

## 4. Comparing Two Molecules

```python
smi_1, smi_2 = (
    "CNC(=O)Cc1nc(-c2ccccc2)cs1",
    "CCNCCc1nc2ccccc2s1",
)
diff = compare_mols(smi_1, smi_2)
diff
```

`compare_mols` returns the functional group differences and any unique alkane fragments for both molecules. This object is ideal for summarising edits before/after a design change.

### Visual Comparison

```python
img = img_grid(draw_compare_mols(smi_1, smi_2), num_columns=2)
img
```

- `draw_compare_mols` wraps `compare_mols`, converts the results into highlight dictionaries, and supplies ready-to-display PIL images.
- For a structural alignment perspective, call `draw_RascalMCES(smi_1, smi_2)` to draw the Rascal MCES matches directly.

### Extending to Multiple Pairs

```python
smiles_pairs = [
    (
        "CC(=O)NC[C@H]1CN(c2ccc(-n3ccc(C#N)c3)c(F)c2)C(=O)O1",
        "CC(=O)NC[C@H]1CN(c2ccc(N3CCN(c4ccc(C#N)cn4)CC3)c(F)c2)C(=O)O1",
    ),
    # additional pairs...
]

for smi_a, smi_b in smiles_pairs:
    print(compare_mols(smi_a, smi_b))
    img = draw_RascalMCES(smi_a, smi_b, subImgSize=(600, 600))
    compare_img = img_grid(
        draw_compare_mols(smi_a, smi_b, img_size=(1000, 800)),
        num_columns=2,
        cell_height=800,
        cell_width=1000,
    )
```

The loop prints differences, generates MCES imagery, and saves side-by-side comparisons for later review. Adjust `img_size` or `subImgSize` to balance resolution and processing time.

## 5. Injecting User-Defined Functional Groups

To extend the detector with custom SMARTS:

```python
my_fgs_dict = {
    "Cephem": "O=C(O)C1=CCS[C@@H]2CC(=O)N12",
    "Thioguanine": "Nc1nc(=S)c2[nH]cnc2[nH]1",
}
my_afg = AccFG(user_defined_fgs=my_fgs_dict, print_load_info=True)
```

When you pass `user_defined_fgs`, AccFG normalises the SMARTS (canonicalises SMILES and handles `[nH]` vs `[n]`) before merging them into the active dictionary. Existing labels are overwritten, which allows you to tweak built-in definitions without editing CSVs.

### Visualise the New Patterns

```python
cephalosporin_C = "CC(=O)OCC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)CCC[C@H](C(=O)O)N)SC1)C(=O)O"
fgs, fg_graph = my_afg.run(cephalosporin_C, show_atoms=True, show_graph=True)
print_fg_tree(fg_graph, fgs.keys(), show_atom_idx=True)

molimg(
    draw_mol_with_fgs(
        cephalosporin_C,
        afg=my_afg,
        img_size=(900, 900),
    )
)
```

Repeat the pattern for additional molecules (`mol_6_Thioguanosine`, etc.) to verify the custom groups fire as expected. The notebook saves rendered PNGs in `results/` for documentation.


## 6. Next Steps

- Swap in your own SMILES lists to analyse project-specific compounds.
- Extend the user-defined dictionary with domain-specific SMARTS and validate them using the visualisation routines.
- Use the CHEMBL aggregation template to profile other collections (e.g., internal libraries or screening hits).

By following these examples you can replicate the notebook outputs while understanding the reasoning behind each step. The combination of AccFG’s detection logic, RDKit visualisation helpers, and Pandas analytics offers a full pipeline for functional group centric studies.
