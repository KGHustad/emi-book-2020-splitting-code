def dump_parameters(G, filename):
    with open(filename, "w") as f:
        fmt_str_float = '%-20s %g\n'
        fmt_str_int = '%-20s %d\n'
        fmt_str_str = '%-20s %s\n'
        f.write(fmt_str_str % ('bc', G.bc))
        f.write(fmt_str_float % ('Cm', G.Cm))
        f.write(fmt_str_float % ('Cg', G.Cg))
        f.write(fmt_str_int % ('N', G.N))
        f.write(fmt_str_int % ('Nx', G.Nx))
        f.write(fmt_str_int % ('Ny', G.Ny))
        f.write(fmt_str_int % ('Nz', G.Nz))
        f.write(fmt_str_int % ('ne', G.ne))
        f.write(fmt_str_int % ('ni', G.ni))
        f.write(fmt_str_int % ('nv', G.nv))
        f.write(fmt_str_int % ('nw', G.nw))
        f.write(fmt_str_int % ('num_ie_it', G.num_ie_it))
        f.write(fmt_str_int % ('num_ii_it', G.num_ii_it))
        f.write(fmt_str_float % ('dt', G.dt))
        f.write(fmt_str_float % ('dt_ode', G.dt_ode))
        f.write(fmt_str_float % ('dx', G.dx))
        f.write(fmt_str_float % ('dy', G.dy))
        f.write(fmt_str_float % ('dz', G.dz))
        f.write(fmt_str_int % ('num_cells_x', G.num_cells_x))
        f.write(fmt_str_int % ('num_cells_y', G.num_cells_y))
        #f.write(fmt_str_int % ('num_cells_z', G.num_cells_z))

