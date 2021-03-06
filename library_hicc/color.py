#!/usr/bin/env python3

"""
Stores the color/resolution definitions.
"""

def is_red_nelson(gr, stmass, run):
    if run == 'high':
        return gr> 0.675 + 0.02*(np.log10(stmass)-10.28)
    elif run == 'low':
        return gr> 0.625 + 0.02*(np.log10(stmass)-10.28)
    elif run == 'mid':
        return gr> 0.65 + 0.02*(np.log10(stmass)-10.28)

def is_resolved_nelson(stmass, gasmass):
    """
    tests if the subhalo is well-resolved.
    """
    MEANBARYONICMASS=1.4e6 #solar masses  
    refmass = MEANBARYONICMASS*200
    resolved_in_stmass = stmass > refmass
    resolved_in_gasmass = gasmass > refmass
    return resolved_in_gasmass * resolved_in_stmass  