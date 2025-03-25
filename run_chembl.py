from rdkit import Chem
from rdkit.Chem import Draw
from accfg import AccFG
import pandas as pd
from collections import Counter
from chembl_structure_pipeline import standardizer as sdz
from tqdm.auto import tqdm
tqdm.pandas()

afg = AccFG(print_load_info=True)

d = pd.read_csv("../efgs/ch33_inchis.csv",sep='\t', low_memory = False)
d = d[pd.notna(d.standard_inchi)].reset_index(drop = True) # make sure there are no missing inchis
d = d[d.standard_inchi.apply(lambda x: "." not in x.split("/")[1])].reset_index(drop = True) # make sure we have no complex molecules
print(d.shape) # (839063, 1)

d = d.sample(frac=1, random_state=4210).reset_index(drop = True) # scramble the df
n = d.shape[0] # Here one can put an alternative much smaller n 
d = d[d.index.isin(range(n))].reset_index(drop = True) # keep only n
def sdz_standardize_inchi(inchi):
    try:
        return sdz.standardize_mol(Chem.MolFromInchi(inchi))
    except:
        return None

d["mol"] = d.standard_inchi.progress_apply(lambda x: sdz_standardize_inchi(x))
d = d[pd.notna(d.mol)].reset_index(drop = True) # remove entries without molecule object
print(d.shape[0])
d['smiles'] = d['mol'].progress_apply(lambda x: Chem.MolToSmiles(x))
d['fgs'] = d['smiles'].progress_apply(lambda x: afg.run(x, canonical=False))
result_df = d[['smiles', 'fgs']]
result_df.to_csv("molecule_data/chembl_33_fgs.csv", index = False)
def get_count(fg_dict):
    return {name:len(fgs) for name, fgs in fg_dict.items()}
d['fgs_count'] = d['fgs'].apply(lambda x: get_count(x))
# Combine all dictionaries in 'fgs_count' into a single Counter
total_counter = Counter()
for fg_count in d['fgs_count']:
    total_counter.update(fg_count)

# Convert the Counter to a DataFrame and sort by count in descending order
fgs_count_df = pd.DataFrame(total_counter.items(), columns=['Functional Group', 'Count'])
fgs_count_df = fgs_count_df.sort_values(by='Count', ascending=False).reset_index(drop=True)

# Export the DataFrame to a CSV file
fgs_count_df.to_csv("molecule_data/chembl_33_fgs_count.csv", index=False)
