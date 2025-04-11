import festim as F
import openmc2dolfinx

reader = openmc2dolfinx.StructuredGridReader("out.vtk")
# t_production = reader.create_dolfinx_function("mean")


from dolfinx.io import gmshio, XDMFFile
from mpi4py import MPI

# mesh = gmshio.read_from_msh("wedge.msh", comm=MPI.COMM_WORLD)

from create_mesh import gmsh, mesh_comm, model_rank

mesh, ct, _ = gmshio.model_to_mesh(gmsh.model, comm=mesh_comm, rank=model_rank)


with XDMFFile(MPI.COMM_WORLD, "mt.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_meshtags(ct, mesh.geometry)
