import numpy as np
import openmc

source_rate = 1.0  # Define the source rate, adjust as needed


def um_tally_to_vtk(filename: str, statepoint: str, tally_name: str):

    sp = openmc.StatePoint(statepoint)
    tally = sp.get_tally(name=tally_name)
    meshgrid = tally.find_filter(openmc.MeshFilter).mesh
    mean = tally.get_reshaped_data(value='mean', expand_dims=True).squeeze()
    rstd = tally.get_reshaped_data(value='rel_err', expand_dims=True).squeeze()
    std_dev = tally.get_reshaped_data(
        value='std_dev', expand_dims=True).squeeze()
    # # writes the mesh tally as a vtk file
    meshgrid.write_data_to_vtk(
        filename=filename,
        datasets={"mean": np.nan_to_num(mean*source_rate, 0.),
                  "std_dev": np.nan_to_num(std_dev*source_rate, 0.),
                  "rstd": np.nan_to_num(rstd*source_rate, nan=0.)})

    return


if __name__ == "__main__":

    tally_names = ['flux_um', 'tbr_um']

    for t in tally_names:

        um_tally_to_vtk(filename=f'{t}_map.vtk',
                        statepoint='statepoint.100.h5',
                        tally_name=t)
