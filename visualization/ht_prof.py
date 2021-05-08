# does hydrotools depend on colossus?
# how much memory/space should I expect this to take?

from hydrotools.interface import interface as iface_run

num_processes = 1
randomize_order = True
sim = 'tng75' # I think this is the right identifier
snap_idx = 99
verbose = False
paranoid = True
output_path = '/lustre/cosinga/final_fields/'
output_compression = 'gzip'
map_get = True
n_max_extract = None, rank_by = None,
# extractGalaxyData([above parameters])