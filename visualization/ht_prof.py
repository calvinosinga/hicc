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
xtf = ['gr_dust']
shf = ['gr_normal']
iface_run.extractGalaxyData(sim='tng75', snap_idx=99, n_max_extract=10, paranoid=True, output_path='/lustre/cosinga/HI-color/results/', 
            file_suffix='_try', output_compression='gzip',catxt_get=True,catxt_fields=xtf, catsh_get=True, catsh_fields=shf)
