import gmsh
from mpi4py import MPI

mesh_comm = MPI.COMM_WORLD
model_rank = 0

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("modelo_1")

gmsh.merge("wedge.step")

# Synchronize the geometry
gmsh.model.geo.synchronize()


# Add a physical group for the volume
volume_tags = gmsh.model.getEntities(dim=3)  # Get all 3D entities (volumes)
if volume_tags:
    gmsh.model.addPhysicalGroup(
        3, [v[1] for v in volume_tags], tag=1
    )  # Add all volumes to a physical group with tag 1
    gmsh.model.setPhysicalName(3, 1, "Volume")  # Set a name for the physical group

gmsh.model.mesh.generate(3)
gmsh.model.mesh.setOrder(2)

gmsh.model.geo.synchronize()
gmsh.fltk.run()


# save mesh to .msh file
gmsh.write("wedge.msh")

# gmsh.finalize()
