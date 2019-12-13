#!/usr/bin/env python3
"""The tree that we output is currently labelled as SRRXXX
This script will produce another tree *_sp_names_bh that will contain
the actual species names in the tree
"""
import sys

class AnnotateTree:
    def __init__(self):
        self.dict_of_names = {'SRR1795735':'SRR1795735_b_pseudominutum', 'SRR1793320':'SRR1793320_b_min_Mf1.05b',
                 'SRR1793321':'SRR1793321_b_min_rt002', 'SRR1793322':'SRR1793322_b_min_mac703',
                 'SRR1793323':'SRR1793323_b_min_rt351', 'SRR1793324':'SRR1793324_b_psygm_rt141',
                 'SRR1793325':'SRR1793325_b_psygm_bMf10_14b.02', 'SRR1793326':'SRR1793326_b_psygm_PurPflex',
                 'SRR1793327':'SRR1793327_b_psygm_HIAp', 'SRR1795737':'SRR1795737_b_aenigmatum',
                 'SRR8292140':'durusdinium_in_hospite_28_3', 'SRR8292146':'durusdinium_free_28_3',
                 'SRR8292139':'durusdinium_in_hospite_28_4', 'SRR8292145':'durusdinium_free_28_4',
                 'SRR2298893':'cladocopium_c1_mi_27_13_8', 'SRR2298892':'cladocopium_c1_mi_27_13_7',
                 'SRR2298891':'cladocopium_c1_mi_32_13_6', 'SRR2298890':'cladocopium_c1_mi_32_13_5',
                 'SRR1300302':'fugasium_1', 'SRR1300303':'fugasium_2', 
                 'SRR1300304':'fugasium_3', 'SRR1300305':'fugasium_4'}
        # A list of the paths to the trees that we want to annotate
        self.tree_paths = [path for path in sys.argv[1:]]

    def annotate_trees(self):
        for tree_path in self.tree_paths:
            with open(tree_path, 'r') as f:
                tree_file = f.readlines()
            if len(tree_file) != 1:
                raise RuntimeError(f'The tree file is over multiple lines ({len(tree_file)})')
            else:
                tree_file = self._replace_using_dict(string=tree_file[0], replacement_dict=self.dict_of_names)
                with open(f'{tree_path}_named', 'w') as f:
                    f.write(f'{tree_file}')

    def _replace_using_dict(self, string, replacement_dict):
        new_string = string
        for k, v in replacement_dict.items():
            new_string = new_string.replace(k, v)
        return new_string


at = AnnotateTree()
at.annotate_trees()