import numpy as np
import openmc

# Using 2:1 atom ratio of LiF to BeF2, similar to values in
# Seifried, Jeffrey E., et al. ‘A General Approach for Determination of
# Acceptable FLiBe Impurity Concentrations in Fluoride-Salt Cooled High
# Temperature Reactors (FHRs)’. Nuclear Engineering and Design, vol. 343, 2019,
# pp. 85–95, https://doi.org10.1016/j.nucengdes.2018.09.038.
# Also using natural lithium enrichment (~7.5 a% Li6)
flibe_nat = openmc.Material(name="Flibe_nat")
# Flibe_nat.temperature = 700 + 273
flibe_nat.add_element("Be", 0.142857, "ao")
flibe_nat.add_nuclide("Li6", 0.021685, "ao")
flibe_nat.add_nuclide("Li7", 0.264029, "ao")
flibe_nat.add_element("F", 0.571429, "ao")
flibe_nat.set_density("g/cm3", 1.94)

beryllium = openmc.Material(name="Beryllium")
# Estimate Be temperature to be around 100 C
# Be.temperature = 100 + 273
beryllium.add_element("Be", 1.0, "ao")
beryllium.set_density("g/cm3", 1.848)

# Graphite (reactor-grade) from PNNL Materials Compendium (PNNL-15870 Rev2)
graphite = openmc.Material(name="Graphite")
graphite.set_density("g/cm3", 1.7)
graphite.add_element("B", 0.000001, "wo")
graphite.add_element("C", 0.999999, "wo")


def build_libra_xl(
    salt_material=flibe_nat,
    multiplier_material=beryllium,
    multiplier_thickness=2.0 * 2.54,
    reflector_material=graphite,
    salt_thickness=40,
    salt_height=80,
    reflector_thickness=10,
    translation_vector=[0, 0, 0],
):

    ###### Materials ############################

    # Source: PNNL Materials Compendium April 2021
    # PNNL-15870, Rev. 2
    inconel625 = openmc.Material(name="Inconel 625")
    inconel625.set_density("g/cm3", 8.44)
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

    air = openmc.Material(name="Air")
    air.add_element("C", 0.00012399, "wo")
    air.add_element("N", 0.75527, "wo")
    air.add_element("O", 0.23178, "wo")
    air.add_element("Ar", 0.012827, "wo")
    air.set_density("g/cm3", 0.0012)

    beryllium = openmc.Material(name="Beryllium")
    # Estimate Be temperature to be around 100 C
    # Be.temperature = 100 + 273
    beryllium.add_element("Be", 1.0, "ao")
    beryllium.set_density("g/cm3", 1.848)

    lead = openmc.Material(name="Lead")
    lead.add_element("Pb", 1.0, "ao")
    lead.set_density("g/cm3", 11.34)

    # # Stainless Steel 304 from PNNL Materials Compendium (PNNL-15870 Rev2)
    # SS304 = openmc.Material(name="Stainless Steel 304")
    # # SS304.temperature = 700 + 273
    # SS304.add_element('C',  0.000800, "wo")
    # SS304.add_element('Mn', 0.020000, "wo")
    # SS304.add_element('P',  0.000450 , "wo")
    # SS304.add_element('S',  0.000300, "wo")
    # SS304.add_element('Si', 0.010000, "wo")
    # SS304.add_element('Cr', 0.190000, "wo")
    # SS304.add_element('Ni', 0.095000, "wo")
    # SS304.add_element('Fe', 0.683450, "wo")
    # SS304.set_density("g/cm3", 8.00)

    # Stanless Steel 316 from PNNL Materials Compendium
    SS316 = openmc.Material(name="Stainless Steel 316")
    # SS304.temperature = 700 + 273
    SS316.add_element("C", 0.000800, "wo")
    SS316.add_element("Mn", 0.020000, "wo")
    SS316.add_element("P", 0.000450, "wo")
    SS316.add_element("S", 0.000300, "wo")
    SS316.add_element("Si", 0.010000, "wo")
    SS316.add_element("Cr", 0.170000, "wo")
    SS316.add_element("Ni", 0.120000, "wo")
    SS316.add_element("Mo", 0.025000, "wo")
    SS316.add_element("Fe", 0.653450, "wo")
    SS316.set_density("g/cm3", 8.00)

    # Graphite (reactor-grade) from PNNL Materials Compendium (PNNL-15870 Rev2)
    graphite = openmc.Material(name="Graphite")
    graphite.set_density("g/cm3", 1.7)
    graphite.add_element("B", 0.000001, "wo")
    graphite.add_element("C", 0.999999, "wo")

    # Using Microtherm with 1 a% Al2O3, 27 a% ZrO2, and 72 a% SiO2
    # https://www.foundryservice.com/product/microporous-silica-insulating-boards-mintherm-microtherm-1925of-grades/
    firebrick = openmc.Material(name="Firebrick")
    # Estimate average temperature of Firebrick to be around 300 C
    # Firebrick.temperature = 273 + 300
    firebrick.add_element("Al", 0.004, "ao")
    firebrick.add_element("O", 0.666, "ao")
    firebrick.add_element("Si", 0.240, "ao")
    firebrick.add_element("Zr", 0.090, "ao")
    firebrick.set_density("g/cm3", 0.30)

    # High Density Polyethylene
    # Reference:  PNNL Report 15870 (Rev. 1)
    HDPE = openmc.Material(name="HDPE")
    HDPE.set_density("g/cm3", 0.95)
    HDPE.add_element("H", 0.143724, "wo")
    HDPE.add_element("C", 0.856276, "wo")

    # Reference: PNNL Report 15870 (Rev. 1) Low Carbon Steel
    steel_lowC = openmc.Material(name="SteelLowC")
    steel_lowC.add_element("C", 0.0010, "wo")
    steel_lowC.add_element("Mn", 0.0050, "wo")
    steel_lowC.add_element("P", 0.0004, "wo")
    steel_lowC.add_element("S", 0.0005, "wo")
    steel_lowC.add_element("Fe", 0.9931, "wo")
    steel_lowC.set_density("g/cm3", 7.872)

    # Zirconia with the density of Zircar FBD referenced below:
    zirconia = openmc.Material(name="Zirconia")
    zirconia.set_density("g/cm3", 1.4)
    zirconia.add_element("Zr", 1 / 3, "ao")
    zirconia.add_element("O", 2 / 3, "ao")

    # # Zircar FBD zirconia (90% ZrO2, Y2O) insulation with a density of 1.4 g/cm3
    # # Source: https://www.zircarzirconia.com/images/datasheets/ZZ-5000_Rev02_-_ZYFB-3_ZYFB-6___FBD.pdf?type=file
    # # Website: https://www.zircarzirconia.com/products/rigid-materials
    # zircar_fbd = openmc.Material(name='Zircar_FBD')
    # zircar_fbd.set_density('g/cm3', 1.4)
    # zircar_fbd.add_element('')

    # tungsten
    tungsten = openmc.Material(name="Tungsten")
    tungsten.set_density("g/cm3", 19.28)
    tungsten.add_element("W", 1.00, "ao")

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

    # Name: Borated Polyethylene (5% B in via B4C additive)
    # Density: 0.95 g/cm3
    # Reference: PNNL Report 15870 (Rev. 1) but revised to make it 5 wt.% B
    # Describes: General purpose neutron shielding
    BPE = openmc.Material(name="BPE")
    BPE.set_density("g/cm3", 0.95)
    BPE.add_nuclide("H1", 0.1345, "wo")
    BPE.add_element("B", 0.0500, "wo")
    BPE.add_element("C", 0.8155, "wo")

    water = openmc.Material(name="water")
    water.set_density("g/cm3", 0.998)
    water.add_element("H", 2 / 3, "ao")
    water.add_element("O", 1 / 3, "ao")

    materials = openmc.Materials(
        [
            inconel625,
            flibe_nat,
            air,
            beryllium,
            lead,
            graphite,
            firebrick,
            HDPE,
            SS316,
            steel_lowC,
            zirconia,
            tungsten,
            Concrete,
            BPE,
            water,
        ]
    )

    if salt_material:
        materials.append(salt_material)
    if multiplier_material:
        materials.append(multiplier_material)
    if reflector_material:
        materials.append(reflector_material)

    ######## LIBRA Surfaces #################

    libra_wall_th = 0.3175  # 1/8 inch
    gap_thickness = 0.635  # 1/4 inch
    salt_headspace = 8 * 2.54
    multiplier_height = 5.0 * 2 * 2.54
    # shield_thickness = 6*2.54
    support_plate_thickness = 2.54
    # Tank double wall surfaces

    # GEOMETRY

    x0_plane = openmc.XPlane(0.0)
    y0_plane = openmc.YPlane(0.0)

    inner_cyl_1 = openmc.ZCylinder(r=4.5 * 2.54)
    inner_cyl_2 = openmc.ZCylinder(r=4.75 * 2.54)

    z_plane_1 = openmc.ZPlane(0.0)
    z_plane_2 = openmc.ZPlane(libra_wall_th)

    salt_top_plane = openmc.ZPlane(z_plane_2.z0 + salt_height)

    ## Tank top cover surfaces
    tank_top_cover_plane_1 = openmc.ZPlane(salt_top_plane.z0 + salt_headspace)
    tank_top_cover_plane_2 = openmc.ZPlane(
        tank_top_cover_plane_1.z0 + libra_wall_th * 2
    )

    outer_cyl_1 = openmc.ZCylinder(r=inner_cyl_2.r + salt_thickness)
    outer_cyl_2 = openmc.ZCylinder(r=outer_cyl_1.r + libra_wall_th)

    # multiplier surfaces
    source_z_point = np.mean([z_plane_2.z0, salt_top_plane.z0])

    ### Heater reentrant tube surfaces

    # theta measured going counterclockwise from y=0 plane
    # heater_reentrant_1 center: theta = 15 degrees
    # heater_reentrant_2 center: theta = 45 degrees
    # heater_reentrant_3 center: theta = 75 degrees

    heater_R = 1 / 2 * (outer_cyl_1.r + inner_cyl_2.r)

    heater_reentrant_1_in_cyl = openmc.ZCylinder(
        r=2.34,
        x0=heater_R * np.cos(np.deg2rad(15)),
        y0=heater_R * np.sin(np.deg2rad(15)),
    )
    heater_reentrant_1_out_cyl = openmc.ZCylinder(
        r=2.54,
        x0=heater_R * np.cos(np.deg2rad(15)),
        y0=heater_R * np.sin(np.deg2rad(15)),
    )

    heater_reentrant_2_in_cyl = openmc.ZCylinder(
        r=2.34,
        x0=heater_R * np.cos(np.deg2rad(45)),
        y0=heater_R * np.sin(np.deg2rad(45)),
    )
    heater_reentrant_2_out_cyl = openmc.ZCylinder(
        r=2.54,
        x0=heater_R * np.cos(np.deg2rad(45)),
        y0=heater_R * np.sin(np.deg2rad(45)),
    )

    heater_reentrant_3_in_cyl = openmc.ZCylinder(
        r=2.34,
        x0=heater_R * np.cos(np.deg2rad(75)),
        y0=heater_R * np.sin(np.deg2rad(75)),
    )
    heater_reentrant_3_out_cyl = openmc.ZCylinder(
        r=2.54,
        x0=heater_R * np.cos(np.deg2rad(75)),
        y0=heater_R * np.sin(np.deg2rad(75)),
    )
    heater_reentrant_bot_plane_1 = openmc.ZPlane(25.60)
    heater_reentrant_bot_plane_2 = openmc.ZPlane(25.80)

    # print(source_z_point)

    multiplier_top_th = multiplier_thickness

    multiplier_bot_plane = openmc.ZPlane(source_z_point - multiplier_height / 2)
    multiplier_top_plane = openmc.ZPlane(source_z_point + multiplier_height / 2)

    multiplier_inner_cyl = openmc.ZCylinder(r=0.5 * 2.54)
    multiplier_outer_cyl = openmc.ZCylinder(
        r=multiplier_inner_cyl.r + multiplier_thickness
    )

    vacuum_insulation_inner_cyl = openmc.ZCylinder(
        r=multiplier_outer_cyl.r + (1 / 8) * 2.54
    )
    vacuum_insulation_outer_cyl = openmc.ZCylinder(r=inner_cyl_1.r - (1 / 8) * 2.54)

    support_plate_bot_plane = openmc.ZPlane(z_plane_1.z0 - support_plate_thickness)

    reflector_bot_plane = openmc.ZPlane(
        support_plate_bot_plane.z0 - reflector_thickness
    )
    reflector_top_plane = openmc.ZPlane(tank_top_cover_plane_2.z0 + reflector_thickness)
    reflector_outer_cyl = openmc.ZCylinder(r=outer_cyl_2.r + reflector_thickness)

    # print(multiplier_cyl)
    # print(inner_cyl_1)
    ## Void surfaces:

    # void_bot_plane = openmc.ZPlane(-50.0, boundary_type='vacuum')
    # void_top_plane = openmc.ZPlane(salt_gas_tube_top_plane.z0 + 50.0, boundary_type='vacuum')
    # void_cyl = openmc.ZCylinder(r=outer_cyl_4.r+50, boundary_type='vacuum')

    libra_bot_plane = reflector_bot_plane
    libra_top_plane = reflector_top_plane
    libra_out_cyl = reflector_outer_cyl
    ####### Regions and Cells ###############

    inner_side_wall_out_reg = (
        +inner_cyl_1 & -inner_cyl_2 & +z_plane_2 & -tank_top_cover_plane_1
    )
    inner_side_wall_out_cell = openmc.Cell(
        region=inner_side_wall_out_reg, fill=inconel625, name="inner side wall out"
    )

    outer_side_wall_out_reg = (
        +outer_cyl_1 & -outer_cyl_2 & +z_plane_2 & -tank_top_cover_plane_1
    )
    outer_side_wall_out_cell = openmc.Cell(
        region=outer_side_wall_out_reg, fill=inconel625, name="outer side wall out"
    )

    bot_wall_out_reg = -outer_cyl_2 & +inner_cyl_1 & +z_plane_1 & -z_plane_2
    bot_wall_out_cell = openmc.Cell(
        region=bot_wall_out_reg, fill=inconel625, name="bottom outer wall"
    )

    ## heater reentrant tubes
    heater_reentrant_1_reg = (
        +heater_reentrant_1_in_cyl
        & -heater_reentrant_1_out_cyl
        & -tank_top_cover_plane_2
        & +heater_reentrant_bot_plane_2
    ) | (
        -heater_reentrant_1_out_cyl
        & +heater_reentrant_bot_plane_1
        & -heater_reentrant_bot_plane_2
    )
    heater_reentrant_1_cell = openmc.Cell(
        region=heater_reentrant_1_reg, fill=inconel625, name="Heater Reentrant Tube 1"
    )
    heater_fill_1_reg = (
        -heater_reentrant_1_in_cyl
        & -tank_top_cover_plane_2
        & +heater_reentrant_bot_plane_2
    )
    heater_fill_1_cell = openmc.Cell(
        region=heater_fill_1_reg, fill=air, name="Heater 1"
    )

    heater_overall_1_reg = (
        -heater_reentrant_1_out_cyl
        & +heater_reentrant_bot_plane_1
        & -tank_top_cover_plane_2
    )

    heater_reentrant_2_reg = (
        +heater_reentrant_2_in_cyl
        & -heater_reentrant_2_out_cyl
        & -tank_top_cover_plane_2
        & +heater_reentrant_bot_plane_2
    ) | (
        -heater_reentrant_2_out_cyl
        & +heater_reentrant_bot_plane_1
        & -heater_reentrant_bot_plane_2
    )
    heater_reentrant_2_cell = openmc.Cell(
        region=heater_reentrant_2_reg, fill=inconel625, name="Heater Reentrant Tube 2"
    )
    heater_fill_2_reg = (
        -heater_reentrant_2_in_cyl
        & -tank_top_cover_plane_2
        & +heater_reentrant_bot_plane_2
    )
    heater_fill_2_cell = openmc.Cell(
        region=heater_fill_2_reg, fill=air, name="Heater 2"
    )

    heater_overall_2_reg = (
        -heater_reentrant_2_out_cyl
        & +heater_reentrant_bot_plane_1
        & -tank_top_cover_plane_2
    )

    heater_reentrant_3_reg = (
        +heater_reentrant_3_in_cyl
        & -heater_reentrant_3_out_cyl
        & -tank_top_cover_plane_2
        & +heater_reentrant_bot_plane_2
    ) | (
        -heater_reentrant_3_out_cyl
        & +heater_reentrant_bot_plane_1
        & -heater_reentrant_bot_plane_2
    )
    heater_reentrant_3_cell = openmc.Cell(
        region=heater_reentrant_3_reg, fill=inconel625, name="Heater Reentrant Tube 3"
    )
    heater_fill_3_reg = (
        -heater_reentrant_3_in_cyl
        & -tank_top_cover_plane_2
        & +heater_reentrant_bot_plane_2
    )
    heater_fill_3_cell = openmc.Cell(
        region=heater_fill_3_reg, fill=air, name="Heater 3"
    )

    heater_overall_3_reg = (
        -heater_reentrant_3_out_cyl
        & +heater_reentrant_bot_plane_1
        & -tank_top_cover_plane_2
    )

    ## Inner tank top cover
    inner_tank_cover_reg = (
        +tank_top_cover_plane_1
        & -tank_top_cover_plane_2
        & +inner_cyl_1
        & -outer_cyl_2
        & +heater_reentrant_1_out_cyl
        & +heater_reentrant_2_out_cyl
        & +heater_reentrant_3_out_cyl
    )

    inner_tank_cover_cell = openmc.Cell(
        region=inner_tank_cover_reg, fill=inconel625, name="Inner Tank Top Cover"
    )
    # top_air_1_reg = +tank_top_cover_plane_2 & -double_wall_top_plane_2\
    #             & +inner_cyl_2 & -outer_cyl_1 \
    #             & +heater_reentrant_1_out_cyl & +heater_reentrant_2_out_cyl \
    #             & +heater_reentrant_3_out_cyl
    # top_air_1_cell = openmc.Cell(region=top_air_1_reg, fill=air, name='top_air_1')

    ## Outer tank top cover
    # Inner region

    ## Salt region and cell
    salt_reg = (
        +inner_cyl_2
        & -outer_cyl_1
        & +z_plane_2
        & -salt_top_plane
        & ~heater_overall_1_reg
        & ~heater_overall_2_reg
        & ~heater_overall_3_reg
    )
    salt_cell = openmc.Cell(region=salt_reg, fill=salt_material, name="Salt")

    salt_vol = (
        np.pi
        * (salt_top_plane.z0 - z_plane_2.z0)
        * (outer_cyl_1.r**2 - inner_cyl_2.r**2)
    ) - 12 * np.pi * (
        salt_top_plane.z0 - heater_reentrant_bot_plane_1.z0
    ) * heater_reentrant_1_out_cyl.r**2

    ## Inner tank air
    inner_tank_air_reg = (
        +inner_cyl_2
        & -outer_cyl_1
        & +salt_top_plane
        & -tank_top_cover_plane_1
        & ~heater_overall_1_reg
        & ~heater_overall_2_reg
        & ~heater_overall_3_reg
    )
    inner_tank_air_cell = openmc.Cell(
        region=inner_tank_air_reg, fill=air, name="Inner Tank Headspace"
    )

    ## Multiplier
    multiplier_reg = (
        +multiplier_inner_cyl
        & -multiplier_outer_cyl
        & +multiplier_bot_plane
        & -multiplier_top_plane
    )
    multiplier_cell = openmc.Cell(
        region=multiplier_reg, fill=multiplier_material, name="Multiplier"
    )

    vacuum_insulation_inner_wall_reg = (
        +multiplier_outer_cyl
        & -vacuum_insulation_inner_cyl
        & +reflector_bot_plane
        & -tank_top_cover_plane_2
    )
    vacuum_insulation_inner_wall_cell = openmc.Cell(
        region=vacuum_insulation_inner_wall_reg,
        fill=SS316,
        name="Vacuum Insulation Inner Wall",
    )

    vacuum_insulation_reg = (
        +vacuum_insulation_inner_cyl
        & -vacuum_insulation_outer_cyl
        & +reflector_bot_plane
        & -tank_top_cover_plane_2
    )
    vacuum_insulation_cell = openmc.Cell(
        region=vacuum_insulation_reg, fill=None, name="Vacuum Insulation"
    )

    vacuum_insulation_outer_wall_reg = (
        +vacuum_insulation_outer_cyl
        & -inner_cyl_1
        & +reflector_bot_plane
        & -tank_top_cover_plane_2
    )
    vacuum_insulation_outer_wall_cell = openmc.Cell(
        region=vacuum_insulation_outer_wall_reg,
        fill=SS316,
        name="Vacuum Insulation Outer Wall",
    )

    center_air_reg = (
        -multiplier_outer_cyl
        & +reflector_bot_plane
        & -tank_top_cover_plane_2
        & ~multiplier_reg
    )
    center_air_cell = openmc.Cell(
        region=center_air_reg, fill=air, name="center cylinder air"
    )

    # libra_outside_air_reg = -outer_cyl_4 & +tank_top_cover_plane_1 & -salt_gas_tube_top_plane & +x_plane_1 & +y_plane_1 \
    #                 & ~outer_tank_cover_reg \
    #                 & ~inner_tank_cover_reg & ~fill_tube_reg \
    #                 & ~heater_overall_1_reg & ~heater_overall_2_reg \
    #                 & ~heater_overall_3_reg \
    #                 & ~salt_gas_tube_1_reg & ~salt_gas_tube_2_reg \
    #                 & ~thermocouple_tube_1_reg & ~thermocouple_tube_2_reg \
    #                 & ~thermocouple_tube_3_reg \
    #                 & ~center_tank_out_wall_reg
    # libra_outside_air_reg = top_air_reg

    # libra_outside_air_cell = openmc.Cell(region=libra_outside_air_reg, fill=air, name='LIBRA Outside Air')

    support_plate_reg = (
        +inner_cyl_1 & -outer_cyl_2 & +support_plate_bot_plane & -z_plane_1
    )
    support_plate_cell = openmc.Cell(
        region=support_plate_reg, fill=steel_lowC, name="Support Plate"
    )

    reflector_reg = (
        (
            +outer_cyl_2
            & -reflector_outer_cyl
            & -tank_top_cover_plane_2
            & +support_plate_bot_plane
        )
        | (-reflector_outer_cyl & +tank_top_cover_plane_2 & -reflector_top_plane)
        | (
            -reflector_outer_cyl
            & +inner_cyl_1
            & -support_plate_bot_plane
            & +reflector_bot_plane
        )
    )
    reflector_cell = openmc.Cell(
        region=reflector_reg, fill=reflector_material, name="Reflector"
    )

    libra_quarter_1_reg = (
        -libra_out_cyl & +libra_bot_plane & -libra_top_plane & +x0_plane & +y0_plane
    )

    libra_quarter_1_cells = [
        outer_side_wall_out_cell,
        inner_side_wall_out_cell,
        bot_wall_out_cell,
        heater_reentrant_1_cell,
        heater_reentrant_2_cell,
        heater_reentrant_3_cell,
        heater_fill_1_cell,
        heater_fill_2_cell,
        heater_fill_3_cell,
        inner_tank_cover_cell,
        # top_air_1_cell,
        salt_cell,
        inner_tank_air_cell,
        multiplier_cell,
        vacuum_insulation_inner_wall_cell,
        vacuum_insulation_cell,
        vacuum_insulation_outer_wall_cell,
        center_air_cell,
        support_plate_cell,
        reflector_cell,
    ]

    libra_quarter_1_universe = openmc.Universe(cells=libra_quarter_1_cells)

    libra_quarter_1_cell = openmc.Cell(
        region=libra_quarter_1_reg, fill=libra_quarter_1_universe, name="Quad 1"
    )

    # temp_reg = -libra_out_cyl & +libra_bot_plane & -libra_top_plane & ~libra_quarter_1_reg
    # temp_cell = openmc.Cell(region=temp_reg, fill=None, name='Temporary Cell')

    libra_quarter_2_reg = libra_quarter_1_reg.rotate((0, 0, 90))
    libra_quarter_2_cell = openmc.Cell(
        region=libra_quarter_2_reg, fill=libra_quarter_1_universe, name="Quad 2"
    )
    libra_quarter_2_cell.rotation = [0.0, 0.0, 90.0]

    libra_quarter_3_reg = libra_quarter_1_reg.rotate((0, 0, 180))
    libra_quarter_3_cell = openmc.Cell(
        region=libra_quarter_3_reg, fill=libra_quarter_1_universe, name="Quad 3"
    )
    libra_quarter_3_cell.rotation = [0.0, 0.0, 180.0]

    libra_quarter_4_reg = libra_quarter_1_reg.rotate((0, 0, 270))
    libra_quarter_4_cell = openmc.Cell(
        region=libra_quarter_4_reg, fill=libra_quarter_1_universe, name="Quad 4"
    )
    libra_quarter_4_cell.rotation = [0.0, 0.0, 270.0]

    libra_universe = openmc.Universe(
        cells=[
            libra_quarter_1_cell,
            libra_quarter_2_cell,
            libra_quarter_3_cell,
            libra_quarter_4_cell,
        ]
    )
    # libra_universe = openmc.Universe(cells=[libra_quarter_1_cell, temp_cell])

    libra_reg = -libra_out_cyl & +libra_bot_plane & -libra_top_plane
    libra_reg = libra_reg.translate(translation_vector)
    libra_system_cell = openmc.Cell(fill=libra_universe, region=libra_reg, name="LIBRA")
    libra_system_cell.translation = translation_vector

    point = openmc.stats.Point(
        (
            0.01 + translation_vector[0],
            0.01 + translation_vector[1],
            source_z_point + translation_vector[2],
        )
    )
    src = openmc.IndependentSource(space=point)
    src.energy = openmc.stats.Discrete([14.1e6], [1.0])
    src.strength = 1.0

    return (
        libra_reg,
        libra_system_cell,
        materials,
        src,
        salt_cell,
        salt_material,
        salt_vol,
    )
