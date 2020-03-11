import numpy as np

class DomainGeometry(object):
    pass

def domain_geometry(num_cells_x, num_cells_y, num_stim_x, num_stim_y,
    stim_start_x, stim_start_y, dx, dy, dz, configuration=None):
    """Generate the geometry of a 3D domain with an number of connected cells
    Input parameters:
        num_cells_x:  Number of cells in the x-direction
        num_cells_y:  Number of cells in the y-direction
        num_stim_x:   Number of cells to stimulate in the x-direction
        num_stim_y:   Number of cells to stimulate in the y-direction
        stim_start_x: First cell to stimulate in the x-direction
        stim_start_y: First cell to stimulate in the y-direction
        dx, dy, dz:   Discretization parameters
    """

    G = DomainGeometry()

    G.dx = dx
    G.dy = dy
    G.dz = dz
    G.num_cells_x = num_cells_x
    G.num_cells_y = num_cells_y
    G.num_cells = G.num_cells_x*G.num_cells_y

    # Set up the cell geometry
    G.cell_length_x = 100e-4
    G.cell_length_y = 18e-4
    G.cell_length_z = 18e-4
    G.j_length_x = 4e-4   # length of Omega_w and Omega_e
    G.j_length_y = 4e-4   # length of Omega_s and Omega_n
    G.j_width_y_x = 8e-4  # width (in x-direction) of gap junction in y-direction
    G.j_width_y_z = 8e-4  # width (in z-direction) of gap junction in y-direction
    G.j_width_x_y = 8e-4  # width (in y-direction) of gap junction in x-direction
    G.j_width_x_z = 8e-4  # width (in z-direction) of gap junction in x-direction

    if configuration == "springer_briefs_splitting_chapter":
        G.cell_length_x = 100e-4
        G.cell_length_y = 18e-4
        G.cell_length_z = 18e-4
        G.j_length_x = 4e-4   # length of Omega_w and Omega_e
        G.j_length_y = 4e-4   # length of Omega_s and Omega_n
        G.j_width_y_x = 10e-4  # width (in x-direction) of gap junction in y-direction
        G.j_width_y_z = 10e-4  # width (in z-direction) of gap junction in y-direction
        G.j_width_x_y = 10e-4  # width (in y-direction) of gap junction in x-direction
        G.j_width_x_z = 10e-4  # width (in z-direction) of gap junction in x-direction


    # Place the cells in space
    distance_to_boundary = 12e-4 # mimimal distance between the intracellular domain
                                # and the boundary of the extracellular space
    distance_to_boundary_z = 4e-4
    G.cell_start_x = distance_to_boundary + G.j_length_x + np.arange(0, num_cells_x)*(G.cell_length_x+2*G.j_length_x)
    G.cell_start_y = distance_to_boundary + G.j_length_y + np.arange(0, num_cells_y)*(G.cell_length_y+2*G.j_length_y)
    G.cell_start_z = distance_to_boundary_z

    # Calculate the domain size
    G.Lx = G.cell_start_x[-1] + G.cell_length_x + G.j_length_x + distance_to_boundary
    G.Ly = G.cell_start_y[-1] + G.cell_length_y + G.j_length_y + distance_to_boundary
    G.Lz = G.cell_start_z + G.cell_length_z + distance_to_boundary_z

    # Compute number of nodes
    G.Nx = int(round(G.Lx/dx))+1
    G.Ny = int(round(G.Ly/dy))+1
    G.Nz = int(round(G.Lz/dz))+1
    G.N = G.Nx*G.Ny*G.Nz

    # Fraction of cell length to stimulate for each cell
    p_stim = np.zeros((G.num_cells_y, G.num_cells_x), dtype=np.float64)

    stim_end_x = stim_start_x + num_stim_x
    stim_end_y = stim_start_y + num_stim_y
    p_stim[stim_start_y:stim_end_y, stim_start_x:stim_end_x] = 1  # Stimulate entire cells
    p_stim = p_stim.reshape(G.num_cells, 1)

    # Fraction of cell length to start stimulating
    s_stim = np.zeros((1, G.num_cells), dtype=np.float64)  # Start stimulating at the beginning of each cell

    # Compute stimulation location
    G.num_stim = np.round(p_stim*(G.cell_length_x+2*G.j_length_x)/G.dx)
    G.stim_start = np.round(s_stim*(G.cell_length_x+2*G.j_length_x)/G.dx)
    G.num_stim = G.num_stim.astype(dtype=np.int64).flatten()
    G.stim_start = G.stim_start.astype(dtype=np.int64).flatten()

    # Set up conductivity
    G.sigma_i = 4   # mS/cm
    G.sigma_e = 20  # mS/cm

    # Define gap junction model (passive)
    G.junc_model = 'passive'
    G.Cg = 0.5       # uF/cm^2
    G.Rg = 0.0045#0.0015  # kOhm cm^2
    G.V_gap = 0    # mV

    # Define ionic model (Grandi)
    G.ion_model = 'Grandi'
    G.Cm = 1.0

    # Define boundary condition
    G.bc = 'Neumann_zero_z'

    return G
