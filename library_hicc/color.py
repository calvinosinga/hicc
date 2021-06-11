#!/usr/bin/env python3
import numpy as np
"""
Stores the color/resolution definitions.
"""

def is_red_nelson(gr, stmass, run):
    """
    Returns a boolean array of the subhalos that are red using
    Dylan's and Benedikt's definition of the color cuts. The run
    determines if these cuts are shifted vertically up or down in the
    gr-stmass plane to test the sensitivity of the results to the color
    definition. Note that this doesn't account for the resolution.
    """
    if run == 'high':
        return gr> 0.7 + 0.02*(np.log10(stmass)-10.28)
    elif run == 'low':
        return gr> 0.6 + 0.02*(np.log10(stmass)-10.28)
    elif run == 'mid':
        return gr> 0.65 + 0.02*(np.log10(stmass)-10.28)

def is_resolved_nelson(stmass, gasmass):
    """
    tests if the subhalo is well-resolved, using the mean baryonic mass that
    was taken from Pillepich et al. 2018 (check to make sure that this is the
    correct value). The subhalo is resolved if the gasmass OR the stellar mass
    is greater than the reference mass.
    """
    MEANBARYONICMASS=1.4e6 #solar masses  
    refmass = MEANBARYONICMASS*200
    resolved_in_stmass = stmass > refmass
    resolved_in_gasmass = gasmass > refmass
    return resolved_in_gasmass + resolved_in_stmass

def is_resolved_stmass(stmass):
    MEANBARYONICMASS=1.4e6 #solar masses  
    refmass = MEANBARYONICMASS*200
    return stmass > refmass

def is_resolved_gas_stmass(stmass, gasmass):
    MEANBARYONICMASS=1.4e6 #solar masses  
    refmass = MEANBARYONICMASS*200
    resolved_in_stmass = stmass > refmass
    resolved_in_gasmass = gasmass > refmass
    return resolved_in_gasmass * resolved_in_stmass
