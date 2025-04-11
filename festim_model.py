import festim as F
import openmc2dolfinx

reader = openmc2dolfinx.StructuredGridReader("out.vtk")
# t_production = reader.create_dolfinx_function("mean")


from dolfinx.io import gmshio, XDMFFile
from mpi4py import MPI


from create_mesh import gmsh, mesh_comm, model_rank

mesh, ct, ft = gmshio.model_to_mesh(gmsh.model, comm=mesh_comm, rank=model_rank)

mesh.geometry.x[:,]


with XDMFFile(MPI.COMM_WORLD, "mt.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_meshtags(ct, mesh.geometry)

with XDMFFile(MPI.COMM_WORLD, "ft.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_meshtags(ft, mesh.geometry)


model = F.HydrogenTransportProblem()
model.volume_meshtags = ct
model.facet_meshtags = ft

salt = F.Material(D_0=100, E_D=0)

top_surface = F.SurfaceSubdomain(id=1)
volume = F.VolumeSubdomain(id=1, material=salt)
model.subdomains = [top_surface, volume]

model.mesh = F.Mesh(mesh)

T = F.Species("T")
model.species = [T]

model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=top_surface, value=0.0, species=T)
]

model.sources = [F.ParticleSource(value=1, volume=volume, species=T)]

model.temperature = 650.0

model.settings = F.Settings(atol=1e0, rtol=1e-10, transient=True, final_time=10)

model.settings.stepsize = 1

model.exports = [
    F.VTXSpeciesExport(filename="tritium_conc.bp", field=T),
]

# from dolfinx import log

# log.set_log_level(log.LogLevel.INFO)

model.initialise()
model.run()
