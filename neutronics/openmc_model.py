import openmc
import pydagmc


def create_mesh_h5m():
    # Assign materials to volumes
    inconel_vol = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    beryllium_vol = [1]
    void_vol = [16]
    flibe_vol = [15]

    # define material map
    group_map = {('mat:inconel625', 2): inconel_vol,
                 ('mat:beryllium', 3): beryllium_vol,
                 ('mat:void', 4): void_vol,
                 ('mat:flibe', 5): flibe_vol}

    # create pydagmc model
    model = pydagmc.Model('libra-one_nosrc.mesh.rtt')
    # assign materials to volumes
    model.add_groups(group_map)
    model.mb.write_file('libra-one_nosrc.mesh.h5m')

    return


def create_mesh_h5m_v1():
    # Assign materials to volumes
    inconel_vol = [1, 2, 3, 4, 5, 6, 7, 13, 14, 15]
    ss304_vol = [10, 11, 12, 16]
    beryllium_vol = [8]
    void_vol = []
    flibe_vol = [9]
    graphite_vol = [17, 18, 19]

    # define material map
    group_map = {('mat:inconel625', 2): inconel_vol,
                 ('mat:beryllium', 3): beryllium_vol,
                 ('mat:void', 4): void_vol,
                 ('mat:flibe', 5): flibe_vol,
                 ('mat:SS304', 6): ss304_vol,
                 ('mat:graphite', 7): graphite_vol}

    # create pydagmc model
    model = pydagmc.Model('libra_one_v1.mesh.rtt')
    # assign materials to volumes
    model.add_groups(group_map)
    model.mb.write_file('libra_one_v1.mesh.h5m')

    return


def create_mesh_h5m_v2():
    # Assign materials to volumes
    inconel_vol = [1, 2, 3, 4, 5, 6, 7, 13, 14, 15]
    ss304_vol = [10, 11, 12, 16]
    beryllium_vol = [8]
    void_vol = [20]
    flibe_vol = [9]
    graphite_vol = [17, 18, 19]

    # define material map
    group_map = {('mat:inconel625', 2): inconel_vol,
                 ('mat:beryllium', 3): beryllium_vol,
                 ('mat:void', 4): void_vol,
                 ('mat:flibe', 5): flibe_vol,
                 ('mat:SS304', 6): ss304_vol,
                 ('mat:graphite', 7): graphite_vol}

    # create pydagmc model
    model = pydagmc.Model('libra_one_v2.mesh.rtt')
    # assign materials to volumes
    model.add_groups(group_map)
    model.mb.write_file('libra_one_v2.mesh.h5m')

    return


def rtt_to_h5m_nomat():
    # create pydagmc model
    model = pydagmc.Model('libra_one_v1.mcnp.mesh.rtt')
    model.mb.write_file('libra_one_v1.mcnp.mesh.h5m')

    return


def openmc_model():

    model = openmc.Model()

    # materials
    inconel625 = openmc.Material(name="inconel625")
    inconel625.add_element("C", 0.000990, "wo")
    inconel625.add_element("Al", 0.003960, "wo")
    inconel625.add_element("Si", 0.004950, "wo")
    inconel625.add_element("P", 0.000148, "wo")
    inconel625.add_element("S", 0.000148, "wo")
    inconel625.add_element("Ti", 0.003960, "wo")
    inconel625.add_element("Cr", 0.215000, "wo")
    inconel625.add_element("Mn", 0.004950, "wo")
    inconel625.add_element("Fe", 0.049495, "wo")
    inconel625.add_element("Co", 0.009899, "wo")
    inconel625.add_element("Ni", 0.580000, "wo")
    inconel625.add_element("Nb", 0.036500, "wo")
    inconel625.add_element("Mo", 0.090000, "wo")
    inconel625.set_density("g/cm3", 8.44)

    # Stainless Steel 304 from PNNL Materials Compendium (PNNL-15870 Rev2)
    SS304 = openmc.Material(name="Stainless Steel 304")
    # SS304.temperature = 700 + 273
    SS304.add_element("C", 0.000800, "wo")
    SS304.add_element("Mn", 0.020000, "wo")
    SS304.add_element("P", 0.000450, "wo")
    SS304.add_element("S", 0.000300, "wo")
    SS304.add_element("Si", 0.010000, "wo")
    SS304.add_element("Cr", 0.190000, "wo")
    SS304.add_element("Ni", 0.095000, "wo")
    SS304.add_element("Fe", 0.683450, "wo")
    SS304.set_density("g/cm3", 8.00)

    # Graphite (reactor-grade) from PNNL Materials Compendium (PNNL-15870 Rev2)
    graphite = openmc.Material(name='Graphite')
    graphite.set_density('g/cm3', 1.7)
    graphite.add_element('B', 0.000001, 'wo')
    graphite.add_element('C', 0.999999, 'wo')

    beryllium = openmc.Material(name="beryllium")
    beryllium.add_element("Be", 1.000000, "wo")
    beryllium.set_density("g/cm3", 1.848)

    flibe = openmc.Material(name="flibe")
    flibe.add_element("Li", 2.0, "ao")
    flibe.add_element("Be", 1.0, "ao")
    flibe.add_element("F", 4.0, "ao")
    flibe.set_density("g/cm3", 1.94)

    model.materials = openmc.Materials([inconel625, beryllium, flibe])

    # geometry
    dag_universe = openmc.DAGMCUniverse(
        "libra-one_nosrc.mesh.h5m").bounded_universe()

    model.geometry = openmc.Geometry(root=dag_universe)

    # settings
    source = openmc.IndependentSource()
    source.space = openmc.stats.Point([0, 0, 0])
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14e6], [1])

    model.settings = openmc.Settings(run_mode='fixed source')
    model.settings.batches = 100
    model.settings.particles = int(1e4)
    model.settings.source = source
    model.settings.output = {'summary': False, 'tallies': False}

    # tallies
    model.tallies = openmc.Tallies()

    unstructured_mesh = openmc.UnstructuredMesh(
        "libra-one_nosrc.mesh.h5m", library='moab')

    unstructured_mesh.output = False
    um_filter = openmc.MeshFilter(unstructured_mesh)

    # tally = openmc.Tally(name='nflux_um')
    # tally.filters = [um_filter]
    # tally.scores = ['flux']
    # model.tallies.append(tally)

    # tally = openmc.Tally(name='tbr_um')
    # tally.filters = [um_filter]
    # tally.scores = ['(n,Xt)']
    # model.tallies.append(tally)

    tally = openmc.Tally(name='tbr_cell')
    tally.filters = [openmc.CellFilter([1, 15])]
    tally.scores = ['(n,Xt)']
    model.tallies.append(tally)

    tally = openmc.Tally(name='tbr_mat')
    tally.filters = [openmc.MaterialFilter([beryllium, flibe])]
    tally.scores = ['(n,Xt)']
    model.tallies.append(tally)

    return model


def main():
    model = openmc_model()
    model.run(threads=8)


if __name__ == "__main__":
    main()
