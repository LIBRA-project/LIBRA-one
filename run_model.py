import openmc
import openmc.model
import numpy as np
import single_wall_model
import openmc_data_downloader as odd

reflector_thickness = 10

breeder_height = 80
breeder_thickness = 40


# Name: Portland concrete
# Density: 2.3 g/cm3
# Reference: PNNL Report 15870 (Rev. 1)
# Describes: facility foundation, floors, walls
Concrete = openmc.Material(name="Concrete")
Concrete.set_density("g/cm3", 2.3)
Concrete.add_nuclide("H1", 0.168759, "ao")
Concrete.add_element("C", 0.001416, "ao")
Concrete.add_nuclide("O16", 0.562524, "ao")
Concrete.add_nuclide("Na23", 0.011838, "ao")
Concrete.add_element("Mg", 0.0014, "ao")
Concrete.add_nuclide("Al27", 0.021354, "ao")
Concrete.add_element("Si", 0.204115, "ao")
Concrete.add_element("K", 0.005656, "ao")
Concrete.add_element("Ca", 0.018674, "ao")
Concrete.add_element("Fe", 0.004264, "ao")

air = openmc.Material(name="Air")
air.add_element("C", 0.00012399, "wo")
air.add_element("N", 0.75527, "wo")
air.add_element("O", 0.23178, "wo")
air.add_element("Ar", 0.012827, "wo")
air.set_density("g/cm3", 0.0012)

floor_rpp = openmc.model.RectangularParallelepiped(-300, 300, -300, 300, -150, -100)
floor_rpp.xmin.boundary_type = "vacuum"
floor_rpp.xmax.boundary_type = "vacuum"
floor_rpp.ymin.boundary_type = "vacuum"
floor_rpp.ymax.boundary_type = "vacuum"
floor_rpp.zmin.boundary_type = "vacuum"
boundary_top_plane = openmc.ZPlane(200, boundary_type="vacuum")

floor_cell = openmc.Cell(region=-floor_rpp, fill=Concrete, name="floor")

settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.batches = 100
settings.inactive = 0
settings.particles = int(1e5)
# settings.volume_calculations = [vol]
# settings.photon_transport = True
settings.photon_transport = False

libra_reg, libra_system_cell, materials, src, salt_cell, salt_material, salt_vol = (
    single_wall_model.build_libra_xl(
        salt_height=breeder_height,
        salt_thickness=breeder_thickness,
        reflector_thickness=reflector_thickness,
    )
)
materials.append(Concrete)
materials.append(air)

outside_air_reg = (
    +floor_rpp.xmin
    & -floor_rpp.xmax
    & +floor_rpp.ymin
    & -floor_rpp.ymax
    & +floor_rpp.zmax
    & -boundary_top_plane
    & ~libra_reg
)
outside_air_cell = openmc.Cell(region=outside_air_reg, fill=air, name="outside air")

universe = openmc.Universe(cells=[libra_system_cell, floor_cell, outside_air_cell])
geometry = openmc.Geometry(universe)

settings.source = src


t_tally = openmc.Tally(name="tritium tally")
salt_filter = openmc.MaterialFilter([salt_material])
salt_cell_filter = openmc.CellFilter([salt_cell])
t_tally.filters.append(salt_filter)
t_tally.filters.append(salt_cell_filter)
t_tally.scores = ["(n,Xt)"]

t_li_tally = openmc.Tally(name="tritium lithium tally")
t_li_tally.filters.append(salt_filter)
t_li_tally.nuclides = ["Li6", "Li7"]
t_li_tally.scores = ["(n,Xt)"]


tallies = openmc.Tallies([t_tally, t_li_tally])

plot1 = openmc.Plot.from_geometry(geometry)
plot1.pixels = (1000, 1500)
plot1.width = [200, 300]
plot1.basis = "xz"
plot1.origin = [0, 0, 25]
plot1.color_by = "material"

plots = openmc.Plots([plot1])

odd.download_cross_section_data(
    materials,
    libraries=["ENDFB-8.0-NNDC"],
    set_OPENMC_CROSS_SECTIONS=True,
    particles=["neutron"],
    destination="cross_sections",
)

model = openmc.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies,
    plots=plots,
)
model.export_to_model_xml()
model.plot_geometry()
model.run()

sp = openmc.StatePoint("statepoint.100.h5")
t_tally = sp.get_tally(name="tritium tally")
tbr = np.sum(t_tally.get_reshaped_data(value="mean").squeeze())
tbr_err = np.sqrt(
    np.sum(np.square(t_tally.get_reshaped_data(value="std_dev").squeeze()))
)

print("TBR = {:.4f} +/- {:.4f}".format(tbr, tbr_err))
