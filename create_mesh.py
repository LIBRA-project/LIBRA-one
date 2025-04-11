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
salt_volume_tag = 1

if volume_tags:
    gmsh.model.addPhysicalGroup(
        3, [v[1] for v in volume_tags], tag=salt_volume_tag
    )  # Add all volumes to a physical group with tag 1
    gmsh.model.setPhysicalName(
        3, salt_volume_tag, "Volume"
    )  # Set a name for the physical group

surface_tags = gmsh.model.getEntities(dim=2)  # Get all 2D entities (surfaces)
print(len(surface_tags))

tags = {
    "top": [1, [0]],
    "wall": [2, [1, 2, 4, 6, 7]],
    "heater": [3, [3, 5]],
}

for name, (tag, surfaces) in tags.items():
    gmsh.model.addPhysicalGroup(
        2, [surface_tags[i][1] for i in surfaces], tag=tag, name=name
    )

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 10)  # Maximum element size

gmsh.model.mesh.generate(3)

# gmsh.model.mesh.setOrder(2)  # uncomment to use quadratic elements

gmsh.model.geo.synchronize()
# gmsh.fltk.run()  # uncomment to show the GUI


# save mesh to .msh file
gmsh.write("wedge.msh")

# gmsh.finalize()
