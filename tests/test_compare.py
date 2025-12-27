from accfg import compare_mols

def test_compare():
    smi_1,smi_2 = ('CNC(=O)Cc1nc(-c2ccccc2)cs1','CCNCCc1nc2ccccc2s1')
    diff = compare_mols(smi_1, smi_2)
    assert diff == (([('secondary amide', 1, [(2, 3, 1)]), 
                      ('benzene', 1, [(8, 9, 10, 11, 12, 13)]), 
                      ('thiazole', 1, [(7, 14, 15, 5, 6)])], 
                     []), 
                    ([('secondary aliphatic amine', 1, [(2,)]), 
                      ('benzo[d]thiazole', 1, [(9, 10, 11, 12, 13, 5, 6, 7, 8)])], 
                     [('C2 alkane', 1, [(3, 4)])]))


    