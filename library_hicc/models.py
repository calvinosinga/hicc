#!/bin/python3

"""
This file stores just simple get methods that return the names of the models that were used
in the HIH2 catalogue.
"""

def get_hiptl_models():
    return ['GD14','GK11','K13','S14']

def get_hisubhalo_models():
    models = ['GD14','GK11','K13','S14']
    proj = ['map','vol']
    res = []
    for m in models:
        for p in proj:
            res.append('m_hi_%s_%s'%(m,p))
    res.append('m_hi_L08_map')
    return res