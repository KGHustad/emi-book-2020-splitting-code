function [mesh, cells] = combine_cells(G, mesh)
%[mesh, cells] = combine_cells(G, mesh) Set up indices for the 
%intracellular points
% 
% Output arguments:
%      mesh: object containing information about the entire mesh
%      cells: collection of objects containing information about the 
%             individual cells
%
% Input arguments:
%      G: object containing information about the domain geometry

% Set up empty arrays for the different types of intracellular grid points
mesh.m_lsw = [];
mesh.m_lse = [];
mesh.m_lnw = [];
mesh.m_lne = [];
mesh.m_hsw = [];
mesh.m_hse = [];
mesh.m_hnw = [];
mesh.m_hne = [];

mesh.m_hw = [];
mesh.m_he = [];
mesh.m_hs = [];
mesh.m_hn = [];
mesh.m_lw = [];
mesh.m_le = [];
mesh.m_ls = [];
mesh.m_ln = [];
mesh.m_ne = [];
mesh.m_sw = [];
mesh.m_se = [];
mesh.m_nw = [];

mesh.m_w = [];
mesh.m_e = [];
mesh.m_s = [];
mesh.m_n = [];
mesh.m_h = [];
mesh.m_l = [];

mesh.i_lsw = [];
mesh.i_lse = [];
mesh.i_lnw = [];
mesh.i_lne = [];
mesh.i_hsw = [];
mesh.i_hse = [];
mesh.i_hnw = [];
mesh.i_hne = [];

mesh.i_hw = [];
mesh.i_he = [];
mesh.i_hs = [];
mesh.i_hn = [];
mesh.i_lw = [];
mesh.i_le = [];
mesh.i_ls = [];
mesh.i_ln = [];
mesh.i_ne = [];
mesh.i_sw = [];
mesh.i_se = [];
mesh.i_nw = [];

mesh.i_w = [];
mesh.i_e = [];
mesh.i_s = [];
mesh.i_n = [];
mesh.i_h = [];
mesh.i_l = [];
mesh.i = [];
mesh.i_all = [];
mesh.to_stim = [];

mesh.gx_e = [];
mesh.gx_he = [];
mesh.gx_le = [];
mesh.gx_se = [];
mesh.gx_ne = [];

mesh.gy_hn = [];
mesh.gy_ln = [];
mesh.gy_nw = [];
mesh.gy_ne = [];
mesh.gy_n = [];

mesh.gx_lse = [];
mesh.gx_lne = [];
mesh.gx_hse = [];
mesh.gx_hne = [];

mesh.gy_lnw = [];
mesh.gy_lne = [];
mesh.gy_hnw = [];
mesh.gy_hne = [];

mesh.gxi_w = [];
mesh.gxi_hw = [];
mesh.gxi_lw = [];
mesh.gxi_sw = [];
mesh.gxi_nw = [];

mesh.gyi_hs = [];
mesh.gyi_ls = [];
mesh.gyi_sw = [];
mesh.gyi_se = [];
mesh.gyi_s = [];

mesh.gxi_lsw = [];
mesh.gxi_lnw = [];
mesh.gxi_hsw = [];
mesh.gxi_hnw = [];

mesh.gyi_lsw = [];
mesh.gyi_lse = [];
mesh.gyi_hsw = [];
mesh.gyi_hse = [];


% Set up cells
for j=1:G.num_cells_y
    for i=1:G.num_cells_x
        
        % Current cell number
        n = (j-1)*G.num_cells_x + i;
        
        % Set up the mesh for a single cell
        cell = set_up_cell(G, n);

        % Save indices for the cell
        cells(n).m_lsw = cell.m_lsw;
        cells(n).m_lse = cell.m_lse;
        cells(n).m_lnw = cell.m_lnw;
        cells(n).m_lne = cell.m_lne;
        cells(n).m_hsw = cell.m_hsw;
        cells(n).m_hse = cell.m_hse;
        cells(n).m_hnw = cell.m_hnw;
        cells(n).m_hne = cell.m_hne;

        cells(n).m_hw = cell.m_hw;
        cells(n).m_he = cell.m_he;
        cells(n).m_hs = cell.m_hs;
        cells(n).m_hn = cell.m_hn;
        cells(n).m_lw = cell.m_lw;
        cells(n).m_le = cell.m_le;
        cells(n).m_ls = cell.m_ls;
        cells(n).m_ln = cell.m_ln;
        cells(n).m_ne = cell.m_ne;
        cells(n).m_sw = cell.m_sw;
        cells(n).m_se = cell.m_se;
        cells(n).m_nw = cell.m_nw;

        cells(n).m_w = cell.m_w;
        cells(n).m_e = cell.m_e;
        cells(n).m_s = cell.m_s;
        cells(n).m_n = cell.m_n;
        cells(n).m_h = cell.m_h;
        cells(n).m_l = cell.m_l;

        cells(n).i_lsw = cell.i_lsw;
        cells(n).i_lse = cell.i_lse;
        cells(n).i_lnw = cell.i_lnw;
        cells(n).i_lne = cell.i_lne;
        cells(n).i_hsw = cell.i_hsw;
        cells(n).i_hse = cell.i_hse;
        cells(n).i_hnw = cell.i_hnw;
        cells(n).i_hne = cell.i_hne;

        cells(n).i_hw = cell.i_hw;
        cells(n).i_he = cell.i_he;
        cells(n).i_hs = cell.i_hs;
        cells(n).i_hn = cell.i_hn;
        cells(n).i_lw = cell.i_lw;
        cells(n).i_le = cell.i_le;
        cells(n).i_ls = cell.i_ls;
        cells(n).i_ln = cell.i_ln;
        cells(n).i_ne = cell.i_ne;
        cells(n).i_sw = cell.i_sw;
        cells(n).i_se = cell.i_se;
        cells(n).i_nw = cell.i_nw;

        cells(n).i_w = cell.i_w;
        cells(n).i_e = cell.i_e;
        cells(n).i_s = cell.i_s;
        cells(n).i_n = cell.i_n;
        cells(n).i_h = cell.i_h;
        cells(n).i_l = cell.i_l;
        cells(n).i = cell.i;
        cells(n).i_all = cell.i_all;

        cells(n).gx_w = cell.gx_w;
        cells(n).gx_e = cell.gx_e;
        cells(n).gx_hw = cell.gx_hw;
        cells(n).gx_lw = cell.gx_lw;
        cells(n).gx_sw = cell.gx_sw;
        cells(n).gx_nw = cell.gx_nw;
        cells(n).gx_he = cell.gx_he;
        cells(n).gx_le = cell.gx_le;
        cells(n).gx_se = cell.gx_se;
        cells(n).gx_ne = cell.gx_ne;

        cells(n).gy_hs = cell.gy_hs;
        cells(n).gy_ls = cell.gy_ls;
        cells(n).gy_sw = cell.gy_sw;
        cells(n).gy_se = cell.gy_se;
        cells(n).gy_hn = cell.gy_hn;
        cells(n).gy_ln = cell.gy_ln;
        cells(n).gy_nw = cell.gy_nw;
        cells(n).gy_ne = cell.gy_ne;
        cells(n).gy_n = cell.gy_n;
        cells(n).gy_s = cell.gy_s;

        cells(n).gx_lsw = cell.gx_lsw;
        cells(n).gx_lse = cell.gx_lse;
        cells(n).gx_lnw = cell.gx_lne;
        cells(n).gx_lne = cell.gx_lne;
        cells(n).gx_hsw = cell.gx_hsw;
        cells(n).gx_hse = cell.gx_hse;
        cells(n).gx_hnw = cell.gx_hnw;
        cells(n).gx_hne = cell.gx_hne;

        cells(n).gy_lsw = cell.gy_lsw;
        cells(n).gy_lse = cell.gy_lse;
        cells(n).gy_lnw = cell.gy_lnw;
        cells(n).gy_lne = cell.gy_lne;
        cells(n).gy_hsw = cell.gy_hsw;
        cells(n).gy_hse = cell.gy_hse;
        cells(n).gy_hnw = cell.gy_hnw;
        cells(n).gy_hne = cell.gy_hne;

        cells(n).gxi_w = cell.gxi_w;
        cells(n).gxi_e = cell.gxi_e;
        cells(n).gxi_hw = cell.gxi_hw;
        cells(n).gxi_lw = cell.gxi_lw;
        cells(n).gxi_sw = cell.gxi_sw;
        cells(n).gxi_nw = cell.gxi_nw;
        cells(n).gxi_he = cell.gxi_he;
        cells(n).gxi_le = cell.gxi_le;
        cells(n).gxi_se = cell.gxi_se;
        cells(n).gxi_ne = cell.gxi_ne;

        cells(n).gyi_hs = cell.gyi_hs;
        cells(n).gyi_ls = cell.gyi_ls;
        cells(n).gyi_sw = cell.gyi_sw;
        cells(n).gyi_se = cell.gyi_se;
        cells(n).gyi_hn = cell.gyi_hn;
        cells(n).gyi_ln = cell.gyi_ln;
        cells(n).gyi_nw = cell.gyi_nw;
        cells(n).gyi_ne = cell.gyi_ne;
        cells(n).gyi_n = cell.gyi_n;
        cells(n).gyi_s = cell.gyi_s;

        cells(n).gxi_lsw = cell.gxi_lsw;
        cells(n).gxi_lse = cell.gxi_lse;
        cells(n).gxi_lnw = cell.gxi_lnw;
        cells(n).gxi_lne = cell.gxi_lne;
        cells(n).gxi_hsw = cell.gxi_hsw;
        cells(n).gxi_hse = cell.gxi_hse;
        cells(n).gxi_hnw = cell.gxi_hnw;
        cells(n).gxi_hne = cell.gxi_hne;

        cells(n).gyi_lsw = cell.gyi_lsw;
        cells(n).gyi_lse = cell.gyi_lse;
        cells(n).gyi_lnw = cell.gyi_lnw;
        cells(n).gyi_lne = cell.gyi_lne;
        cells(n).gyi_hsw = cell.gyi_hsw;
        cells(n).gyi_hse = cell.gyi_hse;
        cells(n).gyi_hnw = cell.gyi_hnw;
        cells(n).gyi_hne = cell.gyi_hne;
            


        % Make adjustments for cells on the boundary
        if j == 1 
            cells(n).gy_s = [];
            cells(n).gy_sw = [];
            cells(n).gy_se = [];
            cells(n).gy_ls = [];
            cells(n).gy_hs = [];
            cells(n).gy_hsw = [];
            cells(n).gy_lse = [];
            cells(n).gy_lsw = [];
            cells(n).gy_hse = [];

            cells(n).m_s = [cells(n).m_s, cell.gy_s];
            cells(n).m_sw = [cells(n).m_sw, cell.gy_sw];
            cells(n).m_se = [cells(n).m_se, cell.gy_se];
            cells(n).m_ls = [cells(n).m_ls, cell.gy_ls];
            cells(n).m_hs = [cells(n).m_hs, cell.gy_hs];
            cells(n).m_hsw = [cells(n).m_hsw, cell.gy_hsw];
            cells(n).m_lse = [cells(n).m_lse, cell.gy_lse];
            cells(n).m_lsw = [cells(n).m_lsw, cell.gy_lsw];
            cells(n).m_hse = [cells(n).m_hse, cell.gy_hse];  

            cells(n).gyi_s = [];
            cells(n).gyi_sw = [];
            cells(n).gyi_se = [];
            cells(n).gyi_ls = [];
            cells(n).gyi_hs = [];
            cells(n).gyi_hsw = [];
            cells(n).gyi_lse = [];
            cells(n).gyi_lsw = [];
            cells(n).gyi_hse = [];

            cells(n).i_s = [cells(n).i_s, cell.gyi_s];
            cells(n).i_sw = [cells(n).i_sw, cell.gyi_sw];
            cells(n).i_se = [cells(n).i_se, cell.gyi_se];
            cells(n).i_ls = [cells(n).i_ls, cell.gyi_ls];
            cells(n).i_hs = [cells(n).i_hs, cell.gyi_hs];
            cells(n).i_hsw = [cells(n).i_hsw, cell.gyi_hsw];
            cells(n).i_lse = [cells(n).i_lse, cell.gyi_lse];
            cells(n).i_lsw = [cells(n).i_lsw, cell.gyi_lsw];
            cells(n).i_hse = [cells(n).i_hse, cell.gyi_hse];
        end
        if j == G.num_cells_y 
            cells(n).gy_n = [];
            cells(n).gy_nw = [];
            cells(n).gy_ne = [];
            cells(n).gy_ln = [];
            cells(n).gy_hn = [];
            cells(n).gy_hnw = [];
            cells(n).gy_lne = [];
            cells(n).gy_lnw = [];
            cells(n).gy_hne = [];

            cells(n).m_n = [cells(n).m_n, cell.gy_n];
            cells(n).m_nw = [cells(n).m_nw, cell.gy_nw];
            cells(n).m_ne = [cells(n).m_ne, cell.gy_ne];
            cells(n).m_ln = [cells(n).m_ln, cell.gy_ln];
            cells(n).m_hn = [cells(n).m_hn, cell.gy_hn];
            cells(n).m_hnw = [cells(n).m_hnw, cell.gy_hnw];
            cells(n).m_lne = [cells(n).m_lne, cell.gy_lne];
            cells(n).m_lnw = [cells(n).m_lnw, cell.gy_lnw];
            cells(n).m_hne = [cells(n).m_hne, cell.gy_hne];  

            cells(n).gyi_n = [];
            cells(n).gyi_nw = [];
            cells(n).gyi_ne = [];
            cells(n).gyi_ln = [];
            cells(n).gyi_hn = [];
            cells(n).gyi_hnw = [];
            cells(n).gyi_lne = [];
            cells(n).gyi_lnw = [];
            cells(n).gyi_hne = [];

            cells(n).i_n = [cells(n).i_n, cell.gyi_n];
            cells(n).i_nw = [cells(n).i_nw, cell.gyi_nw];
            cells(n).i_ne = [cells(n).i_ne, cell.gyi_ne];
            cells(n).i_ln = [cells(n).i_ln, cell.gyi_ln];
            cells(n).i_hn = [cells(n).i_hn, cell.gyi_hn];
            cells(n).i_hnw = [cells(n).i_hnw, cell.gyi_hnw];
            cells(n).i_lne = [cells(n).i_lne, cell.gyi_lne];
            cells(n).i_lnw = [cells(n).i_lnw, cell.gyi_lnw];
            cells(n).i_hne = [cells(n).i_hne, cell.gyi_hne];
        end
        if i == 1
            cells(n).gx_w = [];
            cells(n).gx_sw = [];
            cells(n).gx_nw = [];
            cells(n).gx_lw = [];
            cells(n).gx_hw = [];
            cells(n).gx_hsw = [];
            cells(n).gx_lnw = [];
            cells(n).gx_lsw = [];
            cells(n).gx_hnw = [];

            cells(n).m_w = [cells(n).m_w, cell.gx_w];
            cells(n).m_sw = [cells(n).m_sw, cell.gx_sw];
            cells(n).m_nw = [cells(n).m_nw, cell.gx_nw];
            cells(n).m_lw = [cells(n).m_lw, cell.gx_lw];
            cells(n).m_hw = [cells(n).m_hw, cell.gx_hw];
            cells(n).m_hsw = [cells(n).m_hsw, cell.gx_hsw];
            cells(n).m_lnw = [cells(n).m_lnw, cell.gx_lnw];
            cells(n).m_lsw = [cells(n).m_lsw, cell.gx_lsw];
            cells(n).m_hnw = [cells(n).m_hnw, cell.gx_hnw];           

            cells(n).gxi_w = [];
            cells(n).gxi_sw = [];
            cells(n).gxi_nw = [];
            cells(n).gxi_lw = [];
            cells(n).gxi_hw = [];
            cells(n).gxi_hsw = [];
            cells(n).gxi_lnw = [];
            cells(n).gxi_lsw = [];
            cells(n).gxi_hnw = [];

            cells(n).i_w = [cells(n).i_w, cell.gxi_w];
            cells(n).i_sw = [cells(n).i_sw, cell.gxi_sw];
            cells(n).i_nw = [cells(n).i_nw, cell.gxi_nw];
            cells(n).i_lw = [cells(n).i_lw, cell.gxi_lw];
            cells(n).i_hw = [cells(n).i_hw, cell.gxi_hw];
            cells(n).i_hsw = [cells(n).i_hsw, cell.gxi_hsw];
            cells(n).i_lnw = [cells(n).i_lnw, cell.gxi_lnw];
            cells(n).i_lsw = [cells(n).i_lsw, cell.gxi_lsw];
            cells(n).i_hnw = [cells(n).i_hnw, cell.gxi_hnw];

        end
        if i == G.num_cells_x
            cells(n).gx_e = [];
            cells(n).gx_se = [];
            cells(n).gx_ne = [];
            cells(n).gx_le = [];
            cells(n).gx_he = [];
            cells(n).gx_hse = [];
            cells(n).gx_lne = [];
            cells(n).gx_lse = [];
            cells(n).gx_hne = [];

            cells(n).m_e = [cells(n).m_e, cell.gx_e];
            cells(n).m_se = [cells(n).m_se, cell.gx_se];
            cells(n).m_ne = [cells(n).m_ne, cell.gx_ne];
            cells(n).m_le = [cells(n).m_le, cell.gx_le];
            cells(n).m_he = [cells(n).m_he, cell.gx_he];
            cells(n).m_hse = [cells(n).m_hse, cell.gx_hse];
            cells(n).m_lne = [cells(n).m_lne, cell.gx_lne];
            cells(n).m_lse = [cells(n).m_lse, cell.gx_lse];
            cells(n).m_hne = [cells(n).m_hne, cell.gx_hne];           

            cells(n).gxi_e = [];
            cells(n).gxi_se = [];
            cells(n).gxi_ne = [];
            cells(n).gxi_le = [];
            cells(n).gxi_he = [];
            cells(n).gxi_hse = [];
            cells(n).gxi_lne = [];
            cells(n).gxi_lse = [];
            cells(n).gxi_hne = [];

            cells(n).i_e = [cells(n).i_e, cell.gxi_e];
            cells(n).i_se = [cells(n).i_se, cell.gxi_se];
            cells(n).i_ne = [cells(n).i_ne, cell.gxi_ne];
            cells(n).i_le = [cells(n).i_le, cell.gxi_le];
            cells(n).i_he = [cells(n).i_he, cell.gxi_he];
            cells(n).i_hse = [cells(n).i_hse, cell.gxi_hse];
            cells(n).i_lne = [cells(n).i_lne, cell.gxi_lne];
            cells(n).i_lse = [cells(n).i_lse, cell.gxi_lse];
            cells(n).i_hne = [cells(n).i_hne, cell.gxi_hne];

        end

        % Set up membrane 
        cells(n).v = sort([cells(n).m_lsw, cells(n).m_lse, cells(n).m_lnw, cells(n).m_lne, cells(n).m_hsw, cells(n).m_hse, ...
            cells(n).m_hnw, cells(n).m_hne, cells(n).m_hw, cells(n).m_he, cells(n).m_hs, cells(n).m_hn, ...
            cells(n).m_lw, cells(n).m_le, cells(n).m_ls, cells(n).m_ln, cells(n).m_ne, cells(n).m_sw, ...
            cells(n).m_se, cells(n).m_nw, cells(n).m_w, cells(n).m_e, cells(n).m_s, cells(n).m_n, cells(n).m_h, ...
            cells(n).m_l]);

        % Set up i_all
        cells(n).c_all = sort(unique([cells(n).v, cells(n).i, cells(n).i_lsw, ...
            cells(n).i_lse, cells(n).i_lnw, cells(n).i_lne, cells(n).i_hsw, ...
            cells(n).i_hse, cells(n).i_hnw, cells(n).i_hne, cells(n).i_hw, ...
            cells(n).i_he, cells(n).i_hs, cells(n).i_hn, cells(n).i_lw, ...
            cells(n).i_le, cells(n).i_ls, cells(n).i_ln, cells(n).i_ne, ...
            cells(n).i_sw, cells(n).i_se, cells(n).i_nw, cells(n).i_w, ...
            cells(n).i_e, cells(n).i_s, cells(n).i_n, cells(n).i_h, cells(n).i_l, ...
            cells(n).gx_w, cells(n).gx_e, cells(n).gxi_w, cells(n).gxi_sw, ...
            cells(n).gxi_nw, cells(n).gxi_lw, cells(n).gxi_hw, cells(n).gxi_hsw, ...
            cells(n).gxi_lnw, cells(n).gxi_lsw, cells(n).gxi_hnw, cells(n).gxi_e, ...
            cells(n).gxi_se, cells(n).gxi_ne, cells(n).gxi_le, cells(n).gxi_he, ...
            cells(n).gxi_hse, cells(n).gxi_lne, cells(n).gxi_lse, cells(n).gxi_hne, ...
            cells(n).gyi_n, cells(n).gyi_nw, cells(n).gyi_ne, cells(n).gyi_ln, ...
            cells(n).gyi_hn, cells(n).gyi_hnw, cells(n).gyi_lne, cells(n).gyi_lnw, ...
            cells(n).gyi_hne, cells(n).gyi_s, cells(n).gyi_sw, cells(n).gyi_se, ...
            cells(n).gyi_ls, cells(n).gyi_hs, cells(n).gyi_hsw, cells(n).gyi_lse, ...
            cells(n).gyi_lsw, cells(n).gyi_hse, cells(n).gy_n, cells(n).gy_s]));

        cells(n).i = sort(unique([cells(n).i, cells(n).i_lsw, ...
            cells(n).i_lse, cells(n).i_lnw, cells(n).i_lne, cells(n).i_hsw, ...
            cells(n).i_hse, cells(n).i_hnw, cells(n).i_hne, cells(n).i_hw, ...
            cells(n).i_he, cells(n).i_hs, cells(n).i_hn, cells(n).i_lw, ...
            cells(n).i_le, cells(n).i_ls, cells(n).i_ln, cells(n).i_ne, ...
            cells(n).i_sw, cells(n).i_se, cells(n).i_nw, cells(n).i_w, ...
            cells(n).i_e, cells(n).i_s, cells(n).i_n, cells(n).i_h, cells(n).i_l, ...
            cells(n).gxi_w, cells(n).gxi_sw, cells(n).gxi_nw, cells(n).gxi_lw, ...
            cells(n).gxi_hw, cells(n).gxi_hsw, cells(n).gxi_lnw, ...
            cells(n).gxi_lsw, cells(n).gxi_hnw, cells(n).gxi_e, ...
            cells(n).gxi_se, cells(n).gxi_ne, cells(n).gxi_le, cells(n).gxi_he, ...
            cells(n).gxi_hse, cells(n).gxi_lne, cells(n).gxi_lse, cells(n).gxi_hne, ...
            cells(n).gyi_n, cells(n).gyi_nw, cells(n).gyi_ne, cells(n).gyi_ln, ...
            cells(n).gyi_hn, cells(n).gyi_hnw, cells(n).gyi_lne, cells(n).gyi_lnw, ...
            cells(n).gyi_hne, cells(n).gyi_s, cells(n).gyi_sw, cells(n).gyi_se, ...
            cells(n).gyi_ls, cells(n).gyi_hs, cells(n).gyi_hsw, cells(n).gyi_lse, ...
            cells(n).gyi_lsw, cells(n).gyi_hse]));

        % Set up local indices
        vec = zeros(G.N, 1);
        vec(cells(n).v) = 1;
        vec = vec(cells(n).c_all);
        cells(n).v_of_c = find(vec);

        % gr
        vec = zeros(G.N, 1);
        vec(cells(n).gx_e) = 1;
        vec = vec(cells(n).c_all);
        cells(n).gxe_of_c = find(vec);

        vec = zeros(G.N, 1);
        vec(cells(n).gy_n) = 1;
        vec = vec(cells(n).c_all);
        cells(n).gyn_of_c = find(vec);

        % gn
        vec = zeros(G.N, 1);
        vec(cells(n).gx_w) = 1;
        vec = vec(cells(n).c_all);
        cells(n).gxw_of_c = find(vec);

        vec = zeros(G.N, 1);
        vec(cells(n).gy_s) = 1;
        vec = vec(cells(n).c_all);
        cells(n).gys_of_c = find(vec);


        % Set up indices to stimulate
        if G.num_stim(n) == 0
            cells(n).to_stim = [];
        else
            [x, y, z] = meshgrid(cell.j_start_ix+G.stim_start(n):(cell.j_start_ix+G.stim_start(n)+G.num_stim(n)), 1:G.Ny, 1:G.Nz);
            to_stim = sort(reshape(sub2ind([G.Nx,G.Ny,G.Nz], x, y, z), size(x,1)*size(x,2)*size(x, 3), 1));
            indices = zeros(1,G. N);
            indices(to_stim) = indices(to_stim) + 1;
            indices(cells(n).v) = indices(cells(n).v) + 1;   
            cells(n).to_stim = round(find((indices==2)));
        end

        % Save indices for all cells
        mesh.m_lsw = [mesh.m_lsw, cells(n).m_lsw];
        mesh.m_lse = [mesh.m_lse, cells(n).m_lse];
        mesh.m_lnw = [mesh.m_lnw, cells(n).m_lnw];
        mesh.m_lne = [mesh.m_lne, cells(n).m_lne];
        mesh.m_hsw = [mesh.m_hsw, cells(n).m_hsw];
        mesh.m_hse = [mesh.m_hse, cells(n).m_hse];
        mesh.m_hnw = [mesh.m_hnw, cells(n).m_hnw];
        mesh.m_hne = [mesh.m_hne, cells(n).m_hne];

        mesh.m_hw = [mesh.m_hw, cells(n).m_hw];
        mesh.m_he = [mesh.m_he, cells(n).m_he];
        mesh.m_hs = [mesh.m_hs, cells(n).m_hs];
        mesh.m_hn = [mesh.m_hn, cells(n).m_hn];
        mesh.m_lw = [mesh.m_lw, cells(n).m_lw];
        mesh.m_le = [mesh.m_le, cells(n).m_le];
        mesh.m_ls = [mesh.m_ls, cells(n).m_ls];
        mesh.m_ln = [mesh.m_ln, cells(n).m_ln];
        mesh.m_ne = [mesh.m_ne, cells(n).m_ne];
        mesh.m_sw = [mesh.m_sw, cells(n).m_sw];
        mesh.m_se = [mesh.m_se, cells(n).m_se];
        mesh.m_nw = [mesh.m_nw, cells(n).m_nw];

        mesh.m_w = [mesh.m_w, cells(n).m_w];
        mesh.m_e = [mesh.m_e, cells(n).m_e];
        mesh.m_s = [mesh.m_s, cells(n).m_s];
        mesh.m_n = [mesh.m_n, cells(n).m_n];
        mesh.m_h = [mesh.m_h, cells(n).m_h];
        mesh.m_l = [mesh.m_l, cells(n).m_l];

        mesh.i_lsw = [mesh.i_lsw, cells(n).i_lsw];
        mesh.i_lse = [mesh.i_lse, cells(n).i_lse];
        mesh.i_lnw = [mesh.i_lnw, cells(n).i_lnw];
        mesh.i_lne = [mesh.i_lne, cells(n).i_lne];
        mesh.i_hsw = [mesh.i_hsw, cells(n).i_hsw];
        mesh.i_hse = [mesh.i_hse, cells(n).i_hse];
        mesh.i_hnw = [mesh.i_hnw, cells(n).i_hnw];
        mesh.i_hne = [mesh.i_hne, cells(n).i_hne];

        mesh.i_hw = [mesh.i_hw, cells(n).i_hw];
        mesh.i_he = [mesh.i_he, cells(n).i_he];
        mesh.i_hs = [mesh.i_hs, cells(n).i_hs];
        mesh.i_hn = [mesh.i_hn, cells(n).i_hn];
        mesh.i_lw = [mesh.i_lw, cells(n).i_lw];
        mesh.i_le = [mesh.i_le, cells(n).i_le];
        mesh.i_ls = [mesh.i_ls, cells(n).i_ls];
        mesh.i_ln = [mesh.i_ln, cells(n).i_ln];
        mesh.i_ne = [mesh.i_ne, cells(n).i_ne];
        mesh.i_sw = [mesh.i_sw, cells(n).i_sw];
        mesh.i_se = [mesh.i_se, cells(n).i_se];
        mesh.i_nw = [mesh.i_nw, cells(n).i_nw];

        mesh.i_w = [mesh.i_w, cells(n).i_w];
        mesh.i_e = [mesh.i_e, cells(n).i_e];
        mesh.i_s = [mesh.i_s, cells(n).i_s];
        mesh.i_n = [mesh.i_n, cells(n).i_n];
        mesh.i_h = [mesh.i_h, cells(n).i_h];
        mesh.i_l = [mesh.i_l, cells(n).i_l];
        mesh.i = [mesh.i, cells(n).i];
        mesh.i_all = [mesh.i_all, cells(n).i_all];
        mesh.to_stim = [mesh.to_stim, cells(n).to_stim];

        mesh.gx_lse = [mesh.gx_lse, cells(n).gx_lse];
        mesh.gx_lne = [mesh.gx_lne, cells(n).gx_lne];
        mesh.gx_hse = [mesh.gx_hse, cells(n).gx_hse];
        mesh.gx_hne = [mesh.gx_hne, cells(n).gx_hne];

        mesh.gx_he = [mesh.gx_he, cells(n).gx_he];
        mesh.gx_le = [mesh.gx_le, cells(n).gx_le];
        mesh.gx_ne = [mesh.gx_ne, cells(n).gx_ne];
        mesh.gx_se = [mesh.gx_se, cells(n).gx_se];

        mesh.gx_e = [mesh.gx_e, cells(n).gx_e];

        mesh.gy_lnw = [mesh.gy_lnw, cells(n).gy_lnw];
        mesh.gy_lne = [mesh.gy_lne, cells(n).gy_lne];
        mesh.gy_hnw = [mesh.gy_hnw, cells(n).gy_hnw];
        mesh.gy_hne = [mesh.gy_hne, cells(n).gy_hne];

        mesh.gy_hn = [mesh.gy_hn, cells(n).gy_hn];
        mesh.gy_ln = [mesh.gy_ln, cells(n).gy_ln];
        mesh.gy_ne = [mesh.gy_ne, cells(n).gy_ne];
        mesh.gy_nw = [mesh.gy_nw, cells(n).gy_nw];

        mesh.gy_n = [mesh.gy_n, cells(n).gy_n];


        mesh.gxi_lsw = [mesh.gxi_lsw, cells(n).gxi_lsw];
        mesh.i_ls = [mesh.i_ls, cells(n).gxi_lse];
        mesh.gxi_lnw = [mesh.gxi_lnw, cells(n).gxi_lnw];
        mesh.i_ln = [mesh.i_ln, cells(n).gxi_lne];
        mesh.gxi_hsw = [mesh.gxi_hsw, cells(n).gxi_hsw];
        mesh.i_hs = [mesh.i_hs, cells(n).gxi_hse];
        mesh.gxi_hnw = [mesh.gxi_hnw, cells(n).gxi_hnw];
        mesh.i_hn = [mesh.i_hn, cells(n).gxi_hne];

        mesh.gxi_hw = [mesh.gxi_hw, cells(n).gxi_hw];
        mesh.i_h = [mesh.i_h, cells(n).gxi_he];
        mesh.gxi_lw = [mesh.gxi_lw, cells(n).gxi_lw];
        mesh.i_l = [mesh.i_l, cells(n).gxi_le];
        mesh.i_n = [mesh.i_n, cells(n).gxi_ne];
        mesh.gxi_sw = [mesh.gxi_sw, cells(n).gxi_sw];
        mesh.i_s = [mesh.i_s, cells(n).gxi_se];
        mesh.gxi_nw = [mesh.gxi_nw, cells(n).gxi_nw];

        mesh.gxi_w = [mesh.gxi_w, cells(n).gxi_w];
        mesh.i = [mesh.i, cells(n).gxi_e];

        mesh.gyi_lsw = [mesh.gyi_lsw, cells(n).gyi_lsw];
        mesh.gyi_lse = [mesh.gyi_lse, cells(n).gyi_lse];
        mesh.i_lw = [mesh.i_lw, cells(n).gyi_lnw];
        mesh.i_le = [mesh.i_le, cells(n).gyi_lne];
        mesh.gyi_hsw = [mesh.gyi_hsw, cells(n).gyi_hsw];
        mesh.gyi_hse = [mesh.gyi_hse, cells(n).gyi_hse];
        mesh.i_hw = [mesh.i_hw, cells(n).gyi_hnw];
        mesh.i_he = [mesh.i_he, cells(n).gyi_hne];

        mesh.gyi_hs = [mesh.gyi_hs, cells(n).gyi_hs];
        mesh.i_h = [mesh.i_h, cells(n).gyi_hn];
        mesh.gyi_ls = [mesh.gyi_ls, cells(n).gyi_ls];
        mesh.i_l = [mesh.i_l, cells(n).gyi_ln];
        mesh.i_e = [mesh.i_e, cells(n).gyi_ne];
        mesh.gyi_sw = [mesh.gyi_sw, cells(n).gyi_sw];
        mesh.gyi_se = [mesh.gyi_se, cells(n).gyi_se];
        mesh.i_w = [mesh.i_w, cells(n).gyi_nw];

        mesh.gyi_s = [mesh.gyi_s, cells(n).gyi_s];
        mesh.i = [mesh.i, cells(n).gyi_n];
    end
end

