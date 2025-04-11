import festim as F
import openmc2dolfinx
from festim.helpers import nmm_interpolate
import dolfinx
from dolfinx.io import gmshio, XDMFFile
from mpi4py import MPI

from create_mesh import gmsh, mesh_comm, model_rank

irradiation_time = 10  # s
neutron_rate = 1e8  # n/s
percm3_to_perm3 = 1e6

# Read OpenMC tally results
reader = openmc2dolfinx.StructuredGridReader("out.vtk")
t_production = reader.create_dolfinx_function("mean")
mesh_source = t_production.function_space.mesh
mesh_source.geometry.x[:] *= 1e-2  # cm to m
mesh_source.geometry.x[:, 1] += -0.027
mesh_source.geometry.x[:, 2] += -0.45

t_production.x.array[:] *= neutron_rate  # T/n/cm3 to T/s/cm3
t_production.x.array[:] *= percm3_to_perm3  # T/s/cm3 to T/s/m3


mesh, ct, ft = gmshio.model_to_mesh(gmsh.model, comm=mesh_comm, rank=model_rank)

mesh.geometry.x[:] *= 1e-3  # mm to m

# rotate 90 degrees around x axis (switch y and z)
y_values = mesh.geometry.x[:, 1].copy()
z_values = mesh.geometry.x[:, 2].copy()
mesh.geometry.x[:, 1] = z_values
mesh.geometry.x[:, 2] = y_values


V = dolfinx.fem.functionspace(mesh, ("CG", 1))
t_production_on_wedge = dolfinx.fem.Function(V)
t_production_on_wedge.name = "T production"


nmm_interpolate(t_production_on_wedge, t_production)


with XDMFFile(MPI.COMM_WORLD, "mt.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_meshtags(ct, mesh.geometry)

with XDMFFile(MPI.COMM_WORLD, "ft.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_meshtags(ft, mesh.geometry)

with XDMFFile(MPI.COMM_WORLD, "t_production.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(t_production_on_wedge)


# NOTE need to override these methods in ParticleSource until a
# new version of festim is released
class ValueFromOpenMC(F.Value):
    explicit_time_dependent = True
    temperature_dependent = False

    def update(self, t):
        if t < irradiation_time:
            return
        else:
            self.fenics_object.x.array[:] = 0.0


# NOTE need to overrid this because ParticleSource won't accept a F.Value as value
# need to fix in FESTIM
class SourceFromOpenMC(F.ParticleSource):
    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        if isinstance(value, F.Value):
            self._value = value
        else:
            self._value = F.Value(value)


model = F.HydrogenTransportProblem()
model.volume_meshtags = ct
model.facet_meshtags = ft

salt = F.Material(D_0=0.5, E_D=0)

top_surface = F.SurfaceSubdomain(id=1)
volume = F.VolumeSubdomain(id=1, material=salt)
model.subdomains = [top_surface, volume]

model.mesh = F.Mesh(mesh)

T = F.Species("T")
model.species = [T]

model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=top_surface, value=0.0, species=T)
]

model.sources = [
    SourceFromOpenMC(
        value=ValueFromOpenMC(t_production_on_wedge), volume=volume, species=T
    )
]

model.temperature = 650.0

model.settings = F.Settings(atol=1e-10, rtol=1e-10, final_time=100)

model.settings.stepsize = F.Stepsize(
    0.1,
    growth_factor=1.1,
    cutback_factor=0.9,
    target_nb_iterations=4,
    milestones=[irradiation_time],
)

release_rate = F.SurfaceFlux(field=T, surface=top_surface)

model.exports = [
    F.VTXSpeciesExport(filename="tritium_conc.bp", field=T),
    release_rate,
]

# from dolfinx import log

# log.set_log_level(log.LogLevel.INFO)

model.initialise()
model.run()

import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import numpy as np

wedge_angle = 22.5  # degrees
cumulative_release = cumulative_trapezoid(release_rate.data, release_rate.t, initial=0)

release_rate.data = (
    np.array(release_rate.data) * 360 / wedge_angle
)  # convert to release rate in cm2/s

fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

axs[0].plot(release_rate.t, release_rate.data)
axs[0].set_ylabel("Release Rate (/s)")

axs[1].plot(release_rate.t, cumulative_release)
axs[1].set_xlabel("Time (s)")
axs[1].set_ylabel("Cumulative Release (#)")

plt.show()
