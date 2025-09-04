import openmc
from openmc.deplete import d1s
from openmc_model import openmc_model


# Fetch model
model = openmc_model()

# Include regular mesh photon flux tally
mesh = openmc.RegularMesh().from_domain(
    model.geometry,
    dimension=[50, 50, 50],
    # 100 voxels in x and y axis directions and 1 voxel in z as we want a xy plot
)

rm_filter = openmc.MeshFilter(mesh)
photon_filter = openmc.ParticleFilter(['photon'])
tally = openmc.Tally(name='flux_rm')
tally.filters = [photon_filter, rm_filter]
tally.scores = ['flux']
model.tallies.append(tally)

# Set up d1s calculation
model.settings.photon_transport = True
openmc.config['chain_file'] = '/home/segantin/openmc_models/CROSS_SECTIONS/chain_endfb80_sfr.xml'
model.settings.use_decay_photons = True
d1s.prepare_tallies(model)

model.run(cwd='d1s_results/', threads=8)
