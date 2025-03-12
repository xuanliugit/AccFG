from rdkit import Chem
from rdkit.Chem import Draw
from accfg import AccFG
import pandas as pd
from collections import Counter
from tqdm.auto import tqdm
tqdm.pandas()

def canonical_smiles(smi):
    try:
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        return smi
    except:
        return None
     

def run_fgs_csv(path, smiles_col='smiles', save_dir=True, afg=AccFG()):
    data_name = path.split('/')[-1].split('.')[0]
    try:
        df = pd.read_csv(path)
    except:
        df = pd.read_csv(path, sep=';', on_bad_lines='skip')
    smiles_fgs_df = df[[smiles_col]].rename(columns={smiles_col: 'smiles'})
    print(f'Processing {data_name} dataset...')
    smiles_fgs_df['smiles'] = smiles_fgs_df['smiles'].progress_apply(lambda x: canonical_smiles(x))
    smiles_fgs_df = smiles_fgs_df.dropna()
    print('Calculating functional groups...')
    smiles_fgs_df['fgs'] = smiles_fgs_df['smiles'].progress_apply(lambda x: afg.run(x))
    if save_dir:
        smiles_fgs_df.to_csv(f'molecule_data/{data_name}_fgs.csv', index=False)
    return smiles_fgs_df

if __name__ == '__main__':
    path_list = ['molecule_data/bace.csv',
                 'molecule_data/BBBP.csv',
                 'molecule_data/Lipophilicity.csv',
                 'molecule_data/chembl_approved_sm.csv']
    smiles_col_list = ['mol', 
                       'smiles', 
                       'smiles',
                       'Smiles']
    for path, smiles_col in zip(path_list, smiles_col_list):
        run_fgs_csv(path, smiles_col=smiles_col)