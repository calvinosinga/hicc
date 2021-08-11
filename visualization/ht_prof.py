# does hydrotools depend on colossus?
# how much memory/space should I expect this to take?

from hydrotools.interface import interface as iface_run


xtf = ['gr_dust']
shf = ['gr_normal']
# fields from all of the particles
ptf = ['Coordinates','Velocities','ParticleIDs']
testmax = 10
iface_run.extractGalaxyData(sim='tng75', snap_idx=99, n_max_extract=testmax, paranoid=True, output_path='/lustre/cosinga/HI-color/results/', 
            file_suffix='_try', output_compression='gzip',catxt_get=True,catxt_fields=xtf, catsh_get=True, catsh_fields=shf,
            ptlgas_get = True, ptlgas_fields=ptf,
            ptldm_get=True, ptldm_fields=ptf,
            ptlstr_get=True, ptlstr_fields=ptf)

