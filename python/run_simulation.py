import argparse
import logging
import math
import operator
import time

import numpy as np

import util
from lib_util import LibraryNotFoundException
from domain_geometry import domain_geometry
from set_up_mesh import set_up_mesh
from enforce_outer_boundary_conditions import enforce_outer_boundary_conditions
from dump_parameters import dump_parameters
from solve_system import solve_system
from ionic_models import Grandi_NumPy
from timer import TimerCollection, format_output
import linear_solvers

try:
    from ionic_models import Grandi_C
    has_grandi_c_solver = True
except ImportError:
    has_grandi_c_solver = False

try:
    import linear_solvers_viennacl_cpu
    has_viennacl_cpu_solvers = True
except LibraryNotFoundException:
    has_viennacl_cpu_solvers = False


valid_boundary_conditions = ["Dirichlet", "Neumann", "Neumann_zero_z", "Neumann_zero_yz"]
valid_ode_schemes = ["FE", "GRL1"]

# Run a simulation of the EMI model using the splitting approach with a
# symmetric extracellular matrix
Tstop = 0.2        # Total simulation time
#Tstop = 1.0        # Total simulation time
#Tstop = 0.05        # Total simulation time

# Set up discretization parameters
dt = 0.02         # Time step (ms)
dx = 2e-4         # cm
dy = 2e-4         # cm
dz = 2e-4         # cm
dt_ode = min(0.001, dt) # Time step for ODE scheme (ms)
ode_scheme = "FE"
bc = "Neumann_zero_yz"
geometry_configuration = "springer_briefs_splitting_chapter"

# Set up stimulation parameters
num_stim_x = 2    # Number of cells to stimulate in the x-direction
num_stim_y = 2    # Number of cells to stimulate in the y-direction
stim_start_x = 0  # Index of the first cell to stimulate in the x-direction
stim_start_y = 0  # Index of the first cell to stimulate in the y-direction
stim_amp = 1      # Factor with which the original stimulus amplitude is multipled

# Set up geometry parameters
num_cells_x = 3
num_cells_y = 2

num_cells_x = 8
num_cells_y = 8

num_ii_it = 2 # M_it
num_ie_it = 2 # N_it


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("--bc", type=str, default=bc, choices=valid_boundary_conditions, help="outer boundary condition")
parser.add_argument("--Tstop", type=float, default=Tstop, help="end time for simulation")
parser.add_argument("--dt", type=float, default=dt, help="main timestep")
parser.add_argument("--dt_ode", type=float, default=dt_ode, help="timestep for the ionic model (ODE system)")
parser.add_argument("--ode_scheme", type=str, default=ode_scheme, choices=valid_ode_schemes, help="Solver scheme for the ionic model (ODE system)")
parser.add_argument("--num_cells_x", type=int, default=num_cells_x, help="number of cells in the x dimension")
parser.add_argument("--num_cells_y", type=int, default=num_cells_y, help="number of cells in the y dimension")
parser.add_argument("--num_ii_it", type=int, default=num_ii_it, help="iterations for refining intracellular coupling via gap junctions")
parser.add_argument("--num_ie_it", type=int, default=num_ie_it, help="iterations for refining intracellular-extracellular coupling")
parser.add_argument("--progress_update_period_timesteps", type=int, default=10, help="period (in number of timesteps) between progress updates")
#parser.add_argument("--tqdm", dest='use_tqdm', action="store_true", help="use tqdm library to display a progress bar (not recommended for batch jobs)")
parser.add_argument("--no_tqdm", dest='use_tqdm', action="store_false", help="do not use tqdm library to display a progress bar (disabling tqdm is recommended for batch jobs)")
parser.add_argument("--plot", dest='plot', action="store_true", help="plot solution")
#parser.add_argument("--no_plot", dest='plot', action="store_false", help="don't plot solution")
parser.add_argument("--no_stability_check", dest='check_ionic_model_stability', action="store_false", help="do not check that the solution for the ionic model is stable (this check has a minor performance penalty)")
parser.set_defaults(use_tqdm=True)
parser.set_defaults(plot=False)
parser.set_defaults(check_ionic_model_stability=True)
args = parser.parse_args()

Tstop = args.Tstop
bc = args.bc
dt = args.dt
dt_ode = args.dt_ode
ode_scheme = args.ode_scheme
num_cells_x = args.num_cells_x
num_cells_y = args.num_cells_y
num_ii_it = args.num_ii_it
num_ie_it = args.num_ie_it
progress_update_period_timesteps = args.progress_update_period_timesteps
use_tqdm = args.use_tqdm
check_ionic_model_stability = args.check_ionic_model_stability
plot = args.plot

if dt_ode > dt:
    print(f"Reducing dt_ode to dt={dt}")
    dt_ode = dt
nt = int(math.ceil(dt/dt_ode)) # number of substeps
dt_ode = dt / nt # adjust dt_ode to divide dt

#fprintf('%d x %d cells\n', num_cells_x, num_cells_y);
print('%d x %d cells' % (num_cells_x, num_cells_y))

logger = logging.getLogger("EMI")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
#ch.setLevel(logging.DEBUG)
logger.addHandler(ch)

timers = TimerCollection()
pre = time.perf_counter()

# Generate the domain geometry from the geometry parameters
timers.start("domain geometry")
G = domain_geometry(num_cells_x, num_cells_y, num_stim_x,
    num_stim_y, stim_start_x, stim_start_y, dx, dy, dz, configuration=geometry_configuration)
timers.stop("domain geometry")

# Save some additional parameters in G
G.dt = dt
G.dt_ode = dt_ode
G.nt = nt
G.stim_amp = stim_amp
G.DT = G.dt
G.DT_V = G.dt
G.Tstop = Tstop
G.Nt = int(round(G.Tstop/G.dt))
G.bc = bc
G.progress_update_period_timesteps = progress_update_period_timesteps

# Select number of iterations in the splitting algorithm
G.num_ii_it = num_ii_it # M_it
G.num_ie_it = num_ie_it # N_it

print(f"dx={G.dx:g}, dt={G.dt:g}, dt_ode={G.dt_ode:g}, nt={G.nt:d}, Nt={G.Nt:d}, num_ii_it={G.num_ii_it:d}, num_ie_it={G.num_ie_it:d}, ode_scheme={ode_scheme}")
print()

# Set up mesh
timers.start("set_up_mesh")
mesh, cells = set_up_mesh(G)
timers.stop("set_up_mesh")
logger.info("Completed mesh setup")

timers.start("enforce_outer_boundary_conditions")
mesh, cells = enforce_outer_boundary_conditions(G, mesh, cells)
timers.stop("enforce_outer_boundary_conditions")
logger.info("Enforced outer boundary conditions")

G.nv = len(mesh.v)
G.nw = len(mesh.w)
G.ni = len(mesh.i_all)
G.ne = len(mesh.e_all)

dump_parameters(G, 'parameters.txt')

# Set up cell model
if G.ion_model == "Grandi":
    if has_grandi_c_solver:
        ionic_model = Grandi_C()
    else:
        logger.warning("Grandi C solver has not been built. Falling back to numpy.")
        ionic_model = Grandi_NumPy()
else:
    raise Exception(f"Unsupported ion model {G.ion_model}")


extracellular_solver_class = linear_solvers.Scipy_Direct_Solver
extracellular_solver_class = linear_solvers.Scipy_CG_Solver
if has_viennacl_cpu_solvers:
    pass
    #extracellular_solver_class = linear_solvers_viennacl_cpu.ViennaCL_CG_Solver_CPU
    extracellular_solver_class = linear_solvers_viennacl_cpu.ViennaCL_CG_with_AMG_Solver_CPU

print("Solver for extracellular system: ", extracellular_solver_class)
print("Solver for ionic model: ", type(ionic_model))

# Run a simulation
U, V, W, loop_times, total_loop_time = solve_system(
    G,
    mesh,
    cells,
    ionic_model,
    extracellular_solver_class,
    timers,
    logger,
    ode_scheme=ode_scheme,
    use_tqdm=use_tqdm,
    check_ionic_model_stability=check_ionic_model_stability,
)
post = time.perf_counter()
time_taken = post - pre
total_times = timers.total_times.copy()

unprofiled_time = time_taken - timers.total_time
total_times['*Outside profiled regions*'] = unprofiled_time
relative_times = {k: t/time_taken for k,t in total_times.items()}

relative_times_items = list(relative_times.items())
relative_times_items.sort(key=operator.itemgetter(1), reverse=True)
print()
print("{0:44s} {1:>10s} {2:>12s}".format("Step", "Rel. time", "Abs. time"))
for name, relative_time in relative_times_items:
    total_t = total_times[name]
    print("{0:44s} {1:10.5f} {2:10.3f} s".format(name, relative_time, total_t))
print("-"*(40 + 1 + 10 + 1 + 12))
print("{0:44s} {1:10.5f} {2:10.3f} s".format("Total", sum(relative_times.values()), time_taken))
print()



print(f"Time taken: {time_taken:g} s")
print("Mean loop time: {0:g} s".format(loop_times.mean()))
if len(loop_times) < 100:
    print("Loop times:")
    print(loop_times)
print()

print("U", U[-1,:].min(), U[-1,:].mean(), U[-1,:].max())
print("V", V[-1,:].min(), V[-1,:].mean(), V[-1,:].max())
Ue_final = U[-1, mesh.e_all]
print("Ue", Ue_final.min(), Ue_final.mean(), Ue_final.max())


# Plot the solution
if plot:
    from plot_solution import plot_solution
    plot_solution(U[-1,:], V[-1, :], G, mesh)
