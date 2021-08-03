# does hydrotools depend on colossus?
# how much memory/space should I expect this to take?

from hydrotools.interface import interface as iface_run

#num_processes = 1
#randomize_order = True
#sim = 'tng75' # I think this is the right identifier
#snap_idx = 99
#verbose = False
#paranoid = True
#output_path = '/lustre/cosinga/final_fields/'
#output_compression = 'gzip'
# map_get = True
# n_max_extract = None, rank_by = None,
map_fields = ['f_mol_GD14_vol','f_mol_GK11_vol', 'f_mol_K13_vol','f_mol_S14_vol','f_neutral_H', 'gas_rho']
iface_run.extractGalaxyData(sim='tng75', snap_idx=99, paranoid=True, output_path='/lustre/cosinga/HI-color/results/', 
            file_suffix='_final', output_compression='gzip', map_get=True, map_fields=map_fields)
