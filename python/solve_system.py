from typing import Type
import time

import numpy as np
import scipy.sparse.linalg

from set_up_stim_param import set_up_stim_param
from set_up_symmetric_extracellular_matrix import set_up_symmetric_extracellular_matrix
from set_up_membrane_factors import set_up_membrane_factors
from set_up_Im_matrix import set_up_Im_matrix
import linear_solvers

try:
    import tqdm
    has_tqdm = True
except ImportError:
    has_tqdm = False

class UnstableSolutionException(Exception):
    pass

V_stable_min = -150 # mV
V_stable_max = 100 # mV

def solve_system(
    G,
    mesh,
    cells,
    ionic_model,
    extracellular_solver_class: Type[linear_solvers.LinearSolver],
    timers,
    logger,
    ode_scheme = None,
    use_tqdm = True,
    check_ionic_model_stability = True
):
    "Run a simulation of the EMI model"

    # Make sure the discretization parameters are the same in all directions
    if np.var([G.dx, G.dy, G.dz]) > 1e-10:
        raise ValueError(
            f"Please use dx=dy=dz for a symmetric extracellular system. dx: {dx:g}  dy: {dy:g}  dz: {dz:g}"
        )

    timers.start("setup (solve_system)")
    # Set up cell model parameters
    P = set_up_stim_param(ionic_model, G, mesh)

    # Set up initial conditions
    t = 0
    ionic_model.initialise(G.nv, t_start=t, parameters=P)
    V = ionic_model.V  # transmembrane potential
    t = ionic_model.t
    float_t = np.float64
    W = np.zeros(G.nw, dtype=float_t)  # potential over gap-junction
    W_prev = np.zeros_like(W)
    Ui = np.zeros(G.ni, dtype=float_t)  # intracellular potential
    Ue = np.zeros(G.ne, dtype=float_t)  # extracellular potential
    Ig = np.zeros(G.ni, dtype=float_t)  # gap junction currents?

    # Create LU factorisation and rhs-vectors for each cell
    archetype_to_LU = {}
    for n in range(G.num_cells):
        # cells(n).A = set_up_cell_matrix(G, cells(n))

        if not cells[n].archetype in archetype_to_LU:
            archetype_to_LU[cells[n].archetype] = scipy.sparse.linalg.splu(cells[n].A.tocsc())

        cells[n].LU = archetype_to_LU[cells[n].archetype]

        cells[n].b = np.zeros(len(cells[n].c_all), dtype=float_t)
        Ui[cells[n].c_of_i] = V[cells[n].v_idx].mean()

    # Set up matrix and rhs-vector for the extracellular system
    A = set_up_symmetric_extracellular_matrix(G, mesh)
    b = np.zeros(G.ne, dtype=float_t)

    if np.abs(A - A.T).max() > 1e-10:
        raise Exception("Extracellular system not symmetric\n")

    extracellular_solver = extracellular_solver_class(A)

    # Set up factors for the membrane currents in the extracellular system
    membrane_factors = set_up_membrane_factors(G, mesh, float_t=float_t)

    # Set up matrix for computing the membrane currents
    M = set_up_Im_matrix(G, mesh)

    # Specify which time steps to save Ui, Ue and W
    if hasattr(G, "DT"):
        G.num_save = round(G.Tstop / G.DT) + 1
        G.save_step = round(G.DT / G.dt)
    else:
        G.num_save = G.Nt + 1
        G.save_step = 1
        G.DT = G.dt

    # Specify which time steps to save V
    if hasattr(G, "DT_V"):
        G.num_save_v = round(G.Tstop / G.DT_V) + 1
        G.save_step_v = round(G.DT_V / G.dt)
    else:
        G.num_save_v = G.num_save
        G.save_step_v = 1
        G.DT_V = G.DT

    # Matrices for saving the solution
    U_emi = np.zeros((G.num_save, G.N), dtype=float_t)
    V_emi = np.zeros((G.num_save_v, G.nv), dtype=V.dtype)
    W_emi = np.zeros((G.num_save, G.nw), dtype=float_t)
    V_emi[0, :] = V
    W_emi[0, :] = W

    timers.stop("setup (solve_system)")

    # Run simulation
    print("Setup done. Starting simulation...")
    n_print = G.progress_update_period_timesteps

    t1 = time.perf_counter()

    outer_iterator = range(1, G.Nt + 1)
    use_tqdm = has_tqdm and use_tqdm
    if use_tqdm:
        outer_iterator = tqdm.trange(1, G.Nt + 1, leave=False)

    loop_times = []
    timestamp_before_loop = time.perf_counter()

    for n in outer_iterator:
        timestamp_loop_start = time.perf_counter()
        if use_tqdm:
            outer_iterator.set_description(f"t={t:g} ms")

        # Step 1: Solve ode-system for membrane model
        timers.start("solve ionic model")
        ionic_model.forward_n_steps(G.nt, V=V, dt=G.dt_ode)
        t = ionic_model.t
        timers.stop("solve ionic model")
        # print("ionic model t=%g" % ionic_model.t)

        # check that the solution is stable
        if check_ionic_model_stability:
            timers.start("check stability of ionic model")
            ionic_model_solution_is_stable = np.all(np.logical_and(
                V > V_stable_min,
                V < V_stable_max,
            ))
            logger.debug(f"min(V)={V.min()}  mean(V)={V.mean()}  max(V)={V.max()}")
            #print(f"min(V)={V.min()}  mean(V)={V.mean()}  max(V)={V.max()}")
            if not ionic_model_solution_is_stable:
                msg = f"Solution of the ionic model became unstable. Values in V exceed the interval ({V_stable_min}, {V_stable_max}) mV. "
                msg += "\nmin(V)={0:g}  max(V)={1:g}".format(V.min(), V.max())
                raise UnstableSolutionException(msg)
            timers.stop("check stability of ionic model")

        # Store W from previous time step
        timers.start("copy W to W_prev")
        W_prev[:] = W
        timers.stop("copy W to W_prev")

        # Iterations for intracellular-extracellular coupling
        for q in range(G.num_ie_it):

            # Iterations for intracellular-intracellular coupling
            for j in range(G.num_ii_it):

                # Step 2: Intracellular PDEs
                timers.start("update gap junction currents")
                Ig[mesh.w_of_i] = (1 / G.Rg) * W + G.Cg * (
                    W - W_prev
                ) / G.dt  # Approximation of I_{1,2}
                timers.stop("update gap junction currents")

                for k in range(G.num_cells):
                    # Update rhs-vector with current approximation of Ig
                    timers.start("update b for intracellular system")
                    # cells[k].b[:] = 0 # if initialised correctly, this is not needed
                    cells[k].b[cells[k].v_of_c] = (
                        G.Cm / G.dt * (V[cells[k].v_idx] + Ue[cells[k].v_of_e])
                    )
                    cells[k].b[cells[k].gxw_of_c] = Ig[cells[k].gxw_of_i]
                    cells[k].b[cells[k].gys_of_c] = Ig[cells[k].gys_of_i]
                    cells[k].b[cells[k].gxe_of_c] = Ig[cells[k].gxe_of_i]
                    cells[k].b[cells[k].gyn_of_c] = Ig[cells[k].gyn_of_i]
                    timers.stop("update b for intracellular system")

                    # Solve intracellular system
                    timers.start("solve intracellular system")
                    X = cells[k].LU.solve(cells[k].b)
                    timers.stop("solve intracellular system")

                    # Update W
                    timers.start("update W (gap junctions)")
                    W[cells[k].gxw_of_w] = Ui[cells[k].gxw_of_i] - X[cells[k].gxw_of_c]
                    W[cells[k].gys_of_w] = Ui[cells[k].gys_of_i] - X[cells[k].gys_of_c]
                    timers.stop("update W (gap junctions)")

                    # Update Ui
                    timers.start("update Ui from X")
                    Ui[cells[k].c_of_i] = X
                    timers.stop("update Ui from X")

            # Correct Ui for comparison to coupled method
            timers.start("correct Ui for comparison to coupled method")
            Ui[mesh.w_of_i] += W
            timers.stop("correct Ui for comparison to coupled method")

            # Step 3: Extracellular PDEs
            timers.start("update b for extracellular system")
            Im = M @ Ui  # Compute membrane currents from solution of Ui

            # Update rhs-vector
            b[mesh.v_of_e] = membrane_factors * Im[mesh.v_of_i]
            timers.stop("update b for extracellular system")

            # Solve extracellular system
            timers.start("solve extracellular system")
            Ue = extracellular_solver.solve(b)
            timers.stop("solve extracellular system")

        # Step 4: Update V
        timers.start("update V from Ui and Ue")
        V = Ui[mesh.v_of_i] - Ue[mesh.v_of_e]
        timers.stop("update V from Ui and Ue")

        # Update the membrane potential in the cell model
        # This is now done in the call to forward_n_steps
        # ionic_model.V = V

        # Save the solution
        timers.start("save solution")
        if n % G.save_step == 0:
            U_emi[n // G.save_step, mesh.i_all] = Ui
            U_emi[n // G.save_step, mesh.e_all] = Ue
            W_emi[n // G.save_step, :] = W

        if n % G.save_step_v == 0:
            V_emi[n // G.save_step_v, :] = V
        timers.stop("save solution")

        # Estimate remaining simulation time
        if n % n_print == 0:
            # Print current point in time
            s = "t = %g ms. " % (t)

            # Print estimated simulation time
            t2 = time.perf_counter() - t1
            # t2 = toc(t1)                # Time usage for n_print time steps
            t_rem = t2 * (G.Nt - n) / n_print  # Estimated remaining simulation time
            s += "Approximately "
            if t_rem > 86400 * 2:
                s += "%.1f days" % (t_rem / 86400)
            elif t_rem > 3600:
                s += "%.1f h" % (t_rem / 3600)
            elif t_rem > 60:
                s += "%.1f min" % (t_rem / 60)
            else:
                s += "%.1f sec" % t_rem

            s += " remaining..."
            if use_tqdm:
                pass
                #outer_iterator.write(s)
            else:
                print(s)
            t1 = time.perf_counter()
        loop_times.append(time.perf_counter() - timestamp_loop_start)

    total_loop_time = time.perf_counter() - timestamp_before_loop
    loop_times = np.array(loop_times)
    if use_tqdm:
        outer_iterator.close()

    return U_emi, V_emi, W_emi, loop_times, total_loop_time
