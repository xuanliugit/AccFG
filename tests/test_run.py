from accfg import AccFG

def test_lite():
    afg = AccFG(print_load_info=False, lite=True)
    smi = 'CCO'
    fgs = afg.run(smi, show_atoms=True, show_graph=False)
    assert {'hydroxy': [(2,)]} == fgs
    
def test_full_0():
    afg = AccFG(print_load_info=False, lite=False)
    smi = 'CCO'
    fgs = afg.run(smi, show_atoms=True, show_graph=False)
    assert {'primary hydroxyl': [(2,)]} == fgs

def test_full_1():
    afg = AccFG(print_load_info=False, lite=False)
    smi = 'O=C(O)C1=CCS[C@@H]2CC(=O)N12'
    fgs = afg.run(smi, show_atoms=True, show_graph=False)
    assert {'alkene': [(3, 4)],
            'azetidin-2-one': [(10, 9, 8, 7, 11)],
            'carboxylic acid': [(1, 0, 2)],
            'dialkyl thioether': [(6,)]} == fgs

def test_show_atoms_false():
    afg = AccFG(print_load_info=False, lite=True)
    smi = 'OCCCCO'
    fgs = afg.run(smi, show_atoms=False, show_graph=False)
    assert ['hydroxy'] == fgs
    
def test_stress_run_time():
    afg = AccFG(print_load_info=False, lite=True)
    smi = ['O=C(O)C1=CCS[C@@H]2CC(=O)N12']*100
    for s in smi:
        fgs = afg.run(s, show_atoms=True, show_graph=False)
    assert True
    
def test_user_defined_fgs():
    my_fgs_dict = {'Cephem': 'O=C(O)C1=CCS[C@@H]2CC(=O)N12', 'Thioguanine': 'Nc1nc(=S)c2[nH]cnc2[nH]1'}
    my_afg = AccFG(user_defined_fgs=my_fgs_dict,print_load_info=False)
    cephalosporin_C = 'CC(=O)OCC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)CCC[C@H](C(=O)O)N)SC1)C(=O)O'
    fgs = my_afg.run(cephalosporin_C, show_atoms=True, show_graph=False)
    assert {'primary aliphatic amine': [(21,)], 
            'carboxylic acid': [(22, 23, 24)], 
            'carboxylic ester': [(1, 2, 3)], 
            'secondary amide': [(15, 16, 14)], 
            'Cephem': [(8, 7, 9, 6, 5, 27, 26, 25, 13, 11, 12, 10)]} == fgs
   