# Edit Functional Group Lists

AccFG loads functional group definitions from CSV files packaged with the library. You can extend or customise these dictionaries by editing the CSV files directly or by supplying user-defined SMARTS at runtime.

## 1. Built-in CSV dictionaries

The default dictionaries live under `accfg/`:

- `accfg/fgs_common.csv`: General functional groups.
- `accfg/fgs_heterocycle.csv`: Heterocycles and related ring systems.

Both files follow a custom CSV format tailored for `AccFG.csv_to_dict()`:

- Lines beginning with `%` are treated as comments and ignored.
- In non-lite mode, rows whose `Functional Group` column begins with `#` are considered active entries (the `#` is stripped during import).
- Each row must define two columns: `Functional Group` (label) and `SMARTS Pattern` (SMARTS string describing the substructure).

Before editing, create a backup copy so you can revert if a new entry causes parsing errors.

## 2. Adding new entries

1. Open the target CSV in a text editor or spreadsheet tool that preserves plain text formatting (UTF-8, no BOM).
2. Append a new row using the existing headers:

   | Functional Group | SMARTS Pattern        |
   |------------------|-----------------------|
   | Example group    | `C(=O)O[H]`           |

3. Save the file and rerun your script. AccFG will pick up the changes on the next import.

If you run AccFG in lite mode (`AccFG(lite=True)`), rows prefixed with `#` are skipped, so ensure the entries you need are not commented out.

## 3. Supplying user-defined functional groups in code

You can inject additional definitions without touching the CSV files by passing a dictionary to the constructor:

```python
from accfg import AccFG

custom = {
    "Cephem": "O=C(O)C1=CCS[C@@H]2CC(=O)N12",
    "Thioguanine": "Nc1nc(=S)c2[nH]cnc2[nH]1",
}

afg = AccFG(user_defined_fgs=custom, print_load_info=True)
```

The constructor normalises SMILES inputs (converts `[nH]` to `[n]` where required) and merges them with the built-in dictionaries. Named entries override existing labels, making it easy to tweak behaviour temporarily.

## 4. Testing your changes

After modifying the functional group lists, run a simple extraction and inspect the results:

```python
fgs, graph = afg.run("CC(=O)O", show_atoms=True, show_graph=True)
print(fgs)
```

If AccFG raises an error, double-check the SMARTS syntax and the CSV formatting around your new row. Invalid SMARTS patterns are the most common source of failures.
