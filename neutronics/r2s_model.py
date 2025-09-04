import openmc
import openmc.deplete
from openmc_model import openmc_model
import numpy as np
from pathlib import Path

# Fetch model
model = openmc_model()

model.settings.particles = int(1e1)

# r2s depletion
openmc.config['chain_file'] = '/home/segantin/openmc_models/CROSS_SECTIONS/chain_endfb80_sfr.xml'

# deplete all materials & collect all nuclides
all_nuclides = set()
for i, m in enumerate(model.materials):
    model.materials[i].depletable = True
    all_nuclides.update(m.get_nuclides())

umesh = openmc.UnstructuredMesh("libra-one_nosrc.mesh.h5m", library="moab")

flux_in_each_voxel, micro_xs = openmc.deplete.get_microxs_and_flux(
    model=model,
    domains=umesh,
    energies=[0, 30e6],  # one energy bin from 0 to 30MeV
    chain_file=openmc.config['chain_file'],
    # needed otherwise the statepoint file is produced in an unknown temporary directory
    run_kwargs={'cwd': '.'},
    nuclides=list(all_nuclides)  # Convert set to list
)

print('\\DONE 1!\\')
# Read in the unstructured from the statepoint, this contains additional information (centroids and volumes) compared to the umesh object
sp_filename = f'statepoint.100.h5'
sp = openmc.StatePoint(sp_filename)
print('\\DONE 2!\\')
# normally with regular meshes I would get the mesh from the tally
# but with unstructured meshes the tally does not contain the mesh
# however we can get it from the statepoint file
umesh_from_sp = sp.meshes[umesh.id]
# reading a unstructured mesh from the statepoint trigger internal code in the mesh
#  object so that its centroids and volumes become known.
# centroids and volumes are needed for the get_values and write_data_to_vtk steps
centroids = umesh_from_sp.centroids
mesh_vols = umesh_from_sp.volumes

# Calcualte the material volumes for each mesh element from the unstructured mesh
mat_vols = umesh_from_sp.material_volumes(model=model, n_samples=1_000_000)

print('\\DONE 3!\\')

# Make a new fresh material for every tet in the unstructured mesh.
# Assign the material volume for each tet as the volume is needed to deplete the material.
# Get material IDs from my_materials object
material_ids = [mat.id for mat in model.materials]

materials_for_every_mesh_voxel = []
for i in range(len(mat_vols[material_ids[0]])):
    material_id = next(
        (mid for mid in material_ids if mat_vols[mid][i] > 0), None)
    if material_id is None:
        # handle this case if needed
        continue

    material = next(mat for mat in model.materials if mat.id == material_id)

    # Create a new material instance for this voxel
    new_mat = material.clone()
    new_mat.id = i
    # Use the volume of this material in this voxel from mat_vol
    new_mat.volume = mat_vols[material_id][i]
    materials_for_every_mesh_voxel.append(new_mat)

print('\\DONE 4!\\')
# Define irradiation and cooling time steps.
# Set source rates to zero during decay-only steps.
timesteps = np.diff(np.logspace(0.1, 7, num=10))
timesteps = np.insert(timesteps, 0, 10)  # [s] add 10s of shot time
source_rates = np.zeros(len(timesteps))
# n/s  # [n/s] initial source rate // Full power H-mode DD pulse is 7E18 n/s
source_rates[0] = 1e9  # [n/s]

# Perform the activation / depletion / transmutation of all the materials
# constructing the operator, note we pass in the flux and micro xs
operator = openmc.deplete.IndependentOperator(
    materials=openmc.Materials(materials_for_every_mesh_voxel),
    # Flux in each group in [n-cm/src] for each domain
    fluxes=[flux[0] for flux in flux_in_each_voxel],
    micros=micro_xs,
    reduce_chain_level=5,
    normalization_mode="source-rate"
)

integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=timesteps,
    # a 5 second pulse of neutrons followed by 120 seconds of decay
    source_rates=source_rates,
    timestep_units='s'
)

integrator.integrate()

print('\\DONE 5!\\')

# Make a decay photon flux tally on a regular mesh for the photon / gamma dose.
# The tallies will be used in each of the gamma simulations to see the shutdown dose
# Include regular mesh photon flux tally
mesh = openmc.RegularMesh().from_domain(
    model.geometry,
    dimension=[50, 50, 50],
    # 100 voxels in x and y axis directions and 1 voxel in z as we want a xy plot
)

rm_filter = openmc.MeshFilter(mesh)
photon_filter = openmc.ParticleFilter(['photon'])
tally = openmc.Tally(name='gflux_rm')
tally.filters = [photon_filter, rm_filter]
tally.scores = ['flux']
model.tallies.append(tally)

print('\\DONE 6!\\')
# We will collect the gamma source for all cooling time steps
# Extract all the materials and get their gamma emission spectrum
# Turn these gamma spectra into source terms for later use
results = openmc.deplete.Results.from_hdf5("depletion_results.h5")

all_mesh_sources = []
# skip the first time step as it is the irradiation step
for i_cool in range(1, len(timesteps)):
    all_sources = []
    for i, mesh_vol in enumerate(mesh_vols):
        material_id = str(i)

        activated_material = results[i_cool].get_material(material_id)
        activated_material.volume = mesh_vol
        energy = activated_material.get_decay_photon_energy(
            clip_tolerance=1e-6,
            units='Bq',
        )

        if energy:
            strength = energy.integral()
        # for the strength == None case
        else:
            strength = 0

        my_source = openmc.IndependentSource(
            energy=energy,
            particle="photon",
            strength=strength,
            # constraints={'domains':my_material}
        )

        all_sources.append(my_source)

    # Make a mesh source out of the IndependentSource just made in the inner loop
    mesh_source = openmc.MeshSource(
        mesh=umesh_from_sp,
        sources=all_sources,
    )

    all_mesh_sources.append(mesh_source)

# Makes and runs a simulation model for each time meshsource that has been made.
# Make simulation settings for the gamma transport simulation
my_gamma_settings = openmc.Settings()
my_gamma_settings.run_mode = "fixed source"
my_gamma_settings.batches = 100
my_gamma_settings.particles = int(1e3)
my_gamma_settings.output = {'summary': False}
my_gamma_settings.photon_transport = True

all_gamma_sp_filename = []

for mesh_source in all_mesh_sources:

    my_gamma_settings.source = mesh_source

    # here we use the same pristine materials from before neutron irradaiton as the burnup is low
    # and the materials have not changed much so they would not perterb the neutron spectrum significantly
    # you could also use the activated materials from the depletion results but this would significantly slow the simulation down
    model_gamma = openmc.Model(
        model.geometry, model.materials, my_gamma_settings, model.tallies)

    # Make the model for the gamma / photon transport and run the simulation
    # a folder will be made for each photon transposrt, XML files will be saves there as well as the plot
    # Create directories with parents=True to ensure parent directories are created
    output_dir = Path(f"photons/photon_at_time_{i_cool}")
    output_dir.mkdir(parents=True, exist_ok=True)

    gamma_sp_filename = model_gamma.run(cwd=str(output_dir))
    all_gamma_sp_filename.append(gamma_sp_filename)


print('\\DONE 7!\\')
