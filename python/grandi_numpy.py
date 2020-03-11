import numpy as np

# Gotran generated code for the  "grandi" model
def init_state_values_single(**values):
    """
    Initialize state values
    """
    # Init values
    # m=0.003793087414436, h=0.626221949492493, j=0.624553572490432,
    # x_kr=0.0210022533039071, x_ks=0.00428016666258923,
    # x_to_s=0.000440445885642567, y_to_s=0.785115828275182,
    # x_to_f=0.000440438103758954, y_to_f=0.999995844038706,
    # d=2.92407183949469e-06, f=0.995135796703515,
    # f_Ca_Bj=0.0246760872105795, f_Ca_Bsl=0.0152723084239416,
    # Ry_Rr=0.890806040818203, Ry_Ro=7.40481128853622e-07,
    # Ry_Ri=9.07666168960848e-08, Na_Bj=3.4543773303328,
    # Na_Bsl=0.753740951477775, Tn_CL=0.00893455096919132,
    # Tn_CHc=0.117412025936615, Tn_CHm=0.0106160166692932,
    # CaM=0.000295573424135051, Myo_c=0.00192322252438022,
    # Myo_m=0.137560495022823, SRB=0.00217360235649355,
    # SLL_j=0.00740524521680039, SLL_sl=0.00990339304377132,
    # SLH_j=0.0735890020284214, SLH_sl=0.114583623436917,
    # Csqn_b=1.19723145924432, Ca_sr=0.554760499828172,
    # Na_j=8.40537012592918, Na_sl=8.40491910001025, Na_i=8.40513364344858,
    # K_i=120, Ca_j=0.000175882395147342, Ca_sl=0.000106779509977354,
    # Ca_i=8.72509677797499e-05, V_m=-81.4552030512661
    init_values = np.array(
        [
            0.003793087414436,
            0.626221949492493,
            0.624553572490432,
            0.0210022533039071,
            0.00428016666258923,
            0.000440445885642567,
            0.785115828275182,
            0.000440438103758954,
            0.999995844038706,
            2.92407183949469e-06,
            0.995135796703515,
            0.0246760872105795,
            0.0152723084239416,
            0.890806040818203,
            7.40481128853622e-07,
            9.07666168960848e-08,
            3.4543773303328,
            0.753740951477775,
            0.00893455096919132,
            0.117412025936615,
            0.0106160166692932,
            0.000295573424135051,
            0.00192322252438022,
            0.137560495022823,
            0.00217360235649355,
            0.00740524521680039,
            0.00990339304377132,
            0.0735890020284214,
            0.114583623436917,
            1.19723145924432,
            0.554760499828172,
            8.40537012592918,
            8.40491910001025,
            8.40513364344858,
            120,
            0.000175882395147342,
            0.000106779509977354,
            8.72509677797499e-05,
            -81.4552030512661,
        ],
        dtype=np.float_,
    )

    # State indices and limit checker
    state_ind = dict(
        [
            ("m", 0),
            ("h", 1),
            ("j", 2),
            ("x_kr", 3),
            ("x_ks", 4),
            ("x_to_s", 5),
            ("y_to_s", 6),
            ("x_to_f", 7),
            ("y_to_f", 8),
            ("d", 9),
            ("f", 10),
            ("f_Ca_Bj", 11),
            ("f_Ca_Bsl", 12),
            ("Ry_Rr", 13),
            ("Ry_Ro", 14),
            ("Ry_Ri", 15),
            ("Na_Bj", 16),
            ("Na_Bsl", 17),
            ("Tn_CL", 18),
            ("Tn_CHc", 19),
            ("Tn_CHm", 20),
            ("CaM", 21),
            ("Myo_c", 22),
            ("Myo_m", 23),
            ("SRB", 24),
            ("SLL_j", 25),
            ("SLL_sl", 26),
            ("SLH_j", 27),
            ("SLH_sl", 28),
            ("Csqn_b", 29),
            ("Ca_sr", 30),
            ("Na_j", 31),
            ("Na_sl", 32),
            ("Na_i", 33),
            ("K_i", 34),
            ("Ca_j", 35),
            ("Ca_sl", 36),
            ("Ca_i", 37),
            ("V_m", 38),
        ]
    )

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{0} is not a state.".format(state_name))
        ind = state_ind[state_name]

        # Assign value
        init_values[ind] = value

    return init_values


def init_parameter_values_single(**values):
    """
    Initialize parameter values
    """
    # Param values
    # Fjunc=0.11, Fjunc_CaL=0.9, cellLength=100, cellRadius=10.25,
    # distJuncSL=0.5, distSLcyto=0.45, junctionLength=0.16,
    # junctionRadius=0.015, GNa=23, GNaB=0.000597, IbarNaK=1.8,
    # KmKo=1.5, KmNaip=11, Q10KmNai=1.39, Q10NaK=1.63, GKr=0.035,
    # GKp=0.002, GKs=0.0035, pNaK=0.01833, GK1=0.35, Gto=0.13, epi=1,
    # GClB=0.009, GClCa=0.0548125, KdClCa=0.1, GCaL=0.5, Q10CaL=1.8,
    # pCa=0.00054, pK=2.7e-07, pNa=1.5e-08, IbarNCX=4.5, Kdact=0.00015,
    # KmCai=0.00359, KmCao=1.3, KmNai=12.29, KmNao=87.5, Q10NCX=1.57,
    # ksat=0.32, nu=0.27, IbarSLCaP=0.0673, KmPCa=0.0005,
    # Q10SLCaP=2.35, GCaB=0.0005513, Kmf=0.000246, Kmr=1.7, MaxSR=15,
    # MinSR=1, Q10SRCaP=2.6, Vmax_SRCaP=0.0053114, ec50SR=0.45,
    # hillSRCaP=1.787, kiCa=0.5, kim=0.005, koCa=10, kom=0.06, ks=25,
    # Bmax_Naj=7.561, Bmax_Nasl=1.65, koff_na=0.001, kon_na=0.0001,
    # Bmax_CaM=0.024, Bmax_SR=0.0171, Bmax_TnChigh=0.14,
    # Bmax_TnClow=0.07, Bmax_myosin=0.14, koff_cam=0.238,
    # koff_myoca=0.00046, koff_myomg=5.7e-05, koff_sr=0.06,
    # koff_tnchca=3.2e-05, koff_tnchmg=0.00333, koff_tncl=0.0196,
    # kon_cam=34, kon_myoca=13.8, kon_myomg=0.0157, kon_sr=100,
    # kon_tnchca=2.37, kon_tnchmg=0.003, kon_tncl=32.7,
    # Bmax_SLhighj0=0.000165, Bmax_SLhighsl0=0.0134, Bmax_SLlowj0=0.00046,
    # Bmax_SLlowsl0=0.0374, koff_slh=0.03, koff_sll=1.3, kon_slh=100,
    # kon_sll=100, Bmax_Csqn0=0.14, DcaJuncSL=1.64e-06,
    # DcaSLcyto=1.22e-06, J_ca_juncsl=8.2413e-13, J_ca_slmyo=3.7243e-12,
    # koff_csqn=65, kon_csqn=100, DnaJuncSL=1.09e-05, DnaSLcyto=1.79e-05,
    # J_na_juncsl=1.8313e-14, J_na_slmyo=1.6386e-12, Nao=140, Ko=5.4,
    # Cao=1.8, Cli=15, Clo=150, Mgi=1, Cmem=1.381e-10, Frdy=96485,
    # R=8314, Temp=310, stim_amplitude=40.0, stim_duration=1.0,
    # stim_period=1000.0, stim_start=0.0
    init_values = np.array(
        [
            0.11,
            0.9,
            100,
            10.25,
            0.5,
            0.45,
            0.16,
            0.015,
            23,
            0.000597,
            1.8,
            1.5,
            11,
            1.39,
            1.63,
            0.035,
            0.002,
            0.0035,
            0.01833,
            0.35,
            0.13,
            1,
            0.009,
            0.0548125,
            0.1,
            0.5,
            1.8,
            0.00054,
            2.7e-07,
            1.5e-08,
            4.5,
            0.00015,
            0.00359,
            1.3,
            12.29,
            87.5,
            1.57,
            0.32,
            0.27,
            0.0673,
            0.0005,
            2.35,
            0.0005513,
            0.000246,
            1.7,
            15,
            1,
            2.6,
            0.0053114,
            0.45,
            1.787,
            0.5,
            0.005,
            10,
            0.06,
            25,
            7.561,
            1.65,
            0.001,
            0.0001,
            0.024,
            0.0171,
            0.14,
            0.07,
            0.14,
            0.238,
            0.00046,
            5.7e-05,
            0.06,
            3.2e-05,
            0.00333,
            0.0196,
            34,
            13.8,
            0.0157,
            100,
            2.37,
            0.003,
            32.7,
            0.000165,
            0.0134,
            0.00046,
            0.0374,
            0.03,
            1.3,
            100,
            100,
            0.14,
            1.64e-06,
            1.22e-06,
            8.2413e-13,
            3.7243e-12,
            65,
            100,
            1.09e-05,
            1.79e-05,
            1.8313e-14,
            1.6386e-12,
            140,
            5.4,
            1.8,
            15,
            150,
            1,
            1.381e-10,
            96485,
            8314,
            310,
            40.0,
            1.0,
            1000.0,
            0.0,
        ],
        dtype=np.float_,
    )

    # Parameter indices and limit checker
    param_ind = dict(
        [
            ("Fjunc", 0),
            ("Fjunc_CaL", 1),
            ("cellLength", 2),
            ("cellRadius", 3),
            ("distJuncSL", 4),
            ("distSLcyto", 5),
            ("junctionLength", 6),
            ("junctionRadius", 7),
            ("GNa", 8),
            ("GNaB", 9),
            ("IbarNaK", 10),
            ("KmKo", 11),
            ("KmNaip", 12),
            ("Q10KmNai", 13),
            ("Q10NaK", 14),
            ("GKr", 15),
            ("GKp", 16),
            ("GKs", 17),
            ("pNaK", 18),
            ("GK1", 19),
            ("Gto", 20),
            ("epi", 21),
            ("GClB", 22),
            ("GClCa", 23),
            ("KdClCa", 24),
            ("GCaL", 25),
            ("Q10CaL", 26),
            ("pCa", 27),
            ("pK", 28),
            ("pNa", 29),
            ("IbarNCX", 30),
            ("Kdact", 31),
            ("KmCai", 32),
            ("KmCao", 33),
            ("KmNai", 34),
            ("KmNao", 35),
            ("Q10NCX", 36),
            ("ksat", 37),
            ("nu", 38),
            ("IbarSLCaP", 39),
            ("KmPCa", 40),
            ("Q10SLCaP", 41),
            ("GCaB", 42),
            ("Kmf", 43),
            ("Kmr", 44),
            ("MaxSR", 45),
            ("MinSR", 46),
            ("Q10SRCaP", 47),
            ("Vmax_SRCaP", 48),
            ("ec50SR", 49),
            ("hillSRCaP", 50),
            ("kiCa", 51),
            ("kim", 52),
            ("koCa", 53),
            ("kom", 54),
            ("ks", 55),
            ("Bmax_Naj", 56),
            ("Bmax_Nasl", 57),
            ("koff_na", 58),
            ("kon_na", 59),
            ("Bmax_CaM", 60),
            ("Bmax_SR", 61),
            ("Bmax_TnChigh", 62),
            ("Bmax_TnClow", 63),
            ("Bmax_myosin", 64),
            ("koff_cam", 65),
            ("koff_myoca", 66),
            ("koff_myomg", 67),
            ("koff_sr", 68),
            ("koff_tnchca", 69),
            ("koff_tnchmg", 70),
            ("koff_tncl", 71),
            ("kon_cam", 72),
            ("kon_myoca", 73),
            ("kon_myomg", 74),
            ("kon_sr", 75),
            ("kon_tnchca", 76),
            ("kon_tnchmg", 77),
            ("kon_tncl", 78),
            ("Bmax_SLhighj0", 79),
            ("Bmax_SLhighsl0", 80),
            ("Bmax_SLlowj0", 81),
            ("Bmax_SLlowsl0", 82),
            ("koff_slh", 83),
            ("koff_sll", 84),
            ("kon_slh", 85),
            ("kon_sll", 86),
            ("Bmax_Csqn0", 87),
            ("DcaJuncSL", 88),
            ("DcaSLcyto", 89),
            ("J_ca_juncsl", 90),
            ("J_ca_slmyo", 91),
            ("koff_csqn", 92),
            ("kon_csqn", 93),
            ("DnaJuncSL", 94),
            ("DnaSLcyto", 95),
            ("J_na_juncsl", 96),
            ("J_na_slmyo", 97),
            ("Nao", 98),
            ("Ko", 99),
            ("Cao", 100),
            ("Cli", 101),
            ("Clo", 102),
            ("Mgi", 103),
            ("Cmem", 104),
            ("Frdy", 105),
            ("R", 106),
            ("Temp", 107),
            ("stim_amplitude", 108),
            ("stim_duration", 109),
            ("stim_period", 110),
            ("stim_start", 111),
        ]
    )

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{0} is not a parameter.".format(param_name))
        ind = param_ind[state_name]

        # Assign value
        init_values[ind] = value

    return init_values


def state_indices(*states):
    """
    State indices
    """
    state_inds = dict(
        [
            ("m", 0),
            ("h", 1),
            ("j", 2),
            ("x_kr", 3),
            ("x_ks", 4),
            ("x_to_s", 5),
            ("y_to_s", 6),
            ("x_to_f", 7),
            ("y_to_f", 8),
            ("d", 9),
            ("f", 10),
            ("f_Ca_Bj", 11),
            ("f_Ca_Bsl", 12),
            ("Ry_Rr", 13),
            ("Ry_Ro", 14),
            ("Ry_Ri", 15),
            ("Na_Bj", 16),
            ("Na_Bsl", 17),
            ("Tn_CL", 18),
            ("Tn_CHc", 19),
            ("Tn_CHm", 20),
            ("CaM", 21),
            ("Myo_c", 22),
            ("Myo_m", 23),
            ("SRB", 24),
            ("SLL_j", 25),
            ("SLL_sl", 26),
            ("SLH_j", 27),
            ("SLH_sl", 28),
            ("Csqn_b", 29),
            ("Ca_sr", 30),
            ("Na_j", 31),
            ("Na_sl", 32),
            ("Na_i", 33),
            ("K_i", 34),
            ("Ca_j", 35),
            ("Ca_sl", 36),
            ("Ca_i", 37),
            ("V_m", 38),
        ]
    )

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices) > 1:
        return indices
    else:
        return indices[0]


def parameter_indices(*params):
    """
    Parameter indices
    """
    param_inds = dict(
        [
            ("Fjunc", 0),
            ("Fjunc_CaL", 1),
            ("cellLength", 2),
            ("cellRadius", 3),
            ("distJuncSL", 4),
            ("distSLcyto", 5),
            ("junctionLength", 6),
            ("junctionRadius", 7),
            ("GNa", 8),
            ("GNaB", 9),
            ("IbarNaK", 10),
            ("KmKo", 11),
            ("KmNaip", 12),
            ("Q10KmNai", 13),
            ("Q10NaK", 14),
            ("GKr", 15),
            ("GKp", 16),
            ("GKs", 17),
            ("pNaK", 18),
            ("GK1", 19),
            ("Gto", 20),
            ("epi", 21),
            ("GClB", 22),
            ("GClCa", 23),
            ("KdClCa", 24),
            ("GCaL", 25),
            ("Q10CaL", 26),
            ("pCa", 27),
            ("pK", 28),
            ("pNa", 29),
            ("IbarNCX", 30),
            ("Kdact", 31),
            ("KmCai", 32),
            ("KmCao", 33),
            ("KmNai", 34),
            ("KmNao", 35),
            ("Q10NCX", 36),
            ("ksat", 37),
            ("nu", 38),
            ("IbarSLCaP", 39),
            ("KmPCa", 40),
            ("Q10SLCaP", 41),
            ("GCaB", 42),
            ("Kmf", 43),
            ("Kmr", 44),
            ("MaxSR", 45),
            ("MinSR", 46),
            ("Q10SRCaP", 47),
            ("Vmax_SRCaP", 48),
            ("ec50SR", 49),
            ("hillSRCaP", 50),
            ("kiCa", 51),
            ("kim", 52),
            ("koCa", 53),
            ("kom", 54),
            ("ks", 55),
            ("Bmax_Naj", 56),
            ("Bmax_Nasl", 57),
            ("koff_na", 58),
            ("kon_na", 59),
            ("Bmax_CaM", 60),
            ("Bmax_SR", 61),
            ("Bmax_TnChigh", 62),
            ("Bmax_TnClow", 63),
            ("Bmax_myosin", 64),
            ("koff_cam", 65),
            ("koff_myoca", 66),
            ("koff_myomg", 67),
            ("koff_sr", 68),
            ("koff_tnchca", 69),
            ("koff_tnchmg", 70),
            ("koff_tncl", 71),
            ("kon_cam", 72),
            ("kon_myoca", 73),
            ("kon_myomg", 74),
            ("kon_sr", 75),
            ("kon_tnchca", 76),
            ("kon_tnchmg", 77),
            ("kon_tncl", 78),
            ("Bmax_SLhighj0", 79),
            ("Bmax_SLhighsl0", 80),
            ("Bmax_SLlowj0", 81),
            ("Bmax_SLlowsl0", 82),
            ("koff_slh", 83),
            ("koff_sll", 84),
            ("kon_slh", 85),
            ("kon_sll", 86),
            ("Bmax_Csqn0", 87),
            ("DcaJuncSL", 88),
            ("DcaSLcyto", 89),
            ("J_ca_juncsl", 90),
            ("J_ca_slmyo", 91),
            ("koff_csqn", 92),
            ("kon_csqn", 93),
            ("DnaJuncSL", 94),
            ("DnaSLcyto", 95),
            ("J_na_juncsl", 96),
            ("J_na_slmyo", 97),
            ("Nao", 98),
            ("Ko", 99),
            ("Cao", 100),
            ("Cli", 101),
            ("Clo", 102),
            ("Mgi", 103),
            ("Cmem", 104),
            ("Frdy", 105),
            ("R", 106),
            ("Temp", 107),
            ("stim_amplitude", 108),
            ("stim_duration", 109),
            ("stim_period", 110),
            ("stim_start", 111),
        ]
    )

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    if len(indices) > 1:
        return indices
    else:
        return indices[0]


def monitor_indices(*monitored):
    """
    Monitor indices
    """
    monitor_inds = dict(
        [
            ("Vcell", 0),
            ("Vmyo", 1),
            ("Vsr", 2),
            ("Vsl", 3),
            ("Vjunc", 4),
            ("SAjunc", 5),
            ("SAsl", 6),
            ("Fsl", 7),
            ("Fsl_CaL", 8),
            ("mss", 9),
            ("taum", 10),
            ("ah", 11),
            ("bh", 12),
            ("tauh", 13),
            ("hss", 14),
            ("aj", 15),
            ("bj", 16),
            ("tauj", 17),
            ("jss", 18),
            ("I_Na_junc", 19),
            ("I_Na_sl", 20),
            ("I_Na", 21),
            ("I_nabk_junc", 22),
            ("I_nabk_sl", 23),
            ("I_nabk", 24),
            ("sigma", 25),
            ("fnak", 26),
            ("I_nak_junc", 27),
            ("I_nak_sl", 28),
            ("I_nak", 29),
            ("gkr", 30),
            ("xrss", 31),
            ("tauxr", 32),
            ("rkr", 33),
            ("I_kr", 34),
            ("kp_kp", 35),
            ("I_kp_junc", 36),
            ("I_kp_sl", 37),
            ("I_kp", 38),
            ("eks", 39),
            ("gks_junc", 40),
            ("gks_sl", 41),
            ("xsss", 42),
            ("tauxs", 43),
            ("I_ks_junc", 44),
            ("I_ks_sl", 45),
            ("I_ks", 46),
            ("GtoSlow", 47),
            ("GtoFast", 48),
            ("xtoss", 49),
            ("ytoss", 50),
            ("tauxtos", 51),
            ("tauytos", 52),
            ("I_tos", 53),
            ("tauxtof", 54),
            ("tauytof", 55),
            ("I_tof", 56),
            ("I_to", 57),
            ("I_ClCa_junc", 58),
            ("I_ClCa_sl", 59),
            ("I_ClCa", 60),
            ("I_Clbk", 61),
            ("fss", 62),
            ("dss", 63),
            ("taud", 64),
            ("tauf", 65),
            ("fcaCaMSL", 66),
            ("fcaCaj", 67),
            ("ibarca_j", 68),
            ("ibarca_sl", 69),
            ("ibark", 70),
            ("ibarna_j", 71),
            ("ibarna_sl", 72),
            ("I_Ca_junc", 73),
            ("I_Ca_sl", 74),
            ("I_Ca", 75),
            ("I_CaK", 76),
            ("I_CaNa_junc", 77),
            ("I_CaNa_sl", 78),
            ("I_CaNa", 79),
            ("I_Catot", 80),
            ("Ka_junc", 81),
            ("Ka_sl", 82),
            ("s1_junc", 83),
            ("s1_sl", 84),
            ("s2_junc", 85),
            ("s3_junc", 86),
            ("s2_sl", 87),
            ("s3_sl", 88),
            ("I_ncx_junc", 89),
            ("I_ncx_sl", 90),
            ("I_ncx", 91),
            ("I_pca_junc", 92),
            ("I_pca_sl", 93),
            ("I_pca", 94),
            ("I_cabk_junc", 95),
            ("I_cabk_sl", 96),
            ("I_cabk", 97),
            ("kCaSR", 98),
            ("koSRCa", 99),
            ("kiSRCa", 100),
            ("RI", 101),
            ("J_SRCarel", 102),
            ("J_serca", 103),
            ("J_SRleak", 104),
            ("J_CaB_cytosol", 105),
            ("Bmax_SLlowsl", 106),
            ("Bmax_SLlowj", 107),
            ("Bmax_SLhighsl", 108),
            ("Bmax_SLhighj", 109),
            ("J_CaB_junction", 110),
            ("J_CaB_sl", 111),
            ("Bmax_Csqn", 112),
            ("I_Na_tot_junc", 113),
            ("I_Na_tot_sl", 114),
            ("I_Na_tot_sl2", 115),
            ("I_Na_tot_junc2", 116),
            ("I_K_tot", 117),
            ("I_Ca_tot_junc", 118),
            ("I_Ca_tot_sl", 119),
            ("i_Stim", 120),
            ("I_Na_tot", 121),
            ("I_Cl_tot", 122),
            ("I_Ca_tot", 123),
            ("I_tot", 124),
            ("FoRT", 125),
            ("ena_junc", 126),
            ("ena_sl", 127),
            ("ek", 128),
            ("eca_junc", 129),
            ("eca_sl", 130),
            ("ecl", 131),
            ("Qpow", 132),
            ("aki", 133),
            ("bki", 134),
            ("kiss", 135),
            ("I_K1", 136),
            ("dm_dt", 137),
            ("dh_dt", 138),
            ("dj_dt", 139),
            ("dx_kr_dt", 140),
            ("dx_ks_dt", 141),
            ("dx_to_s_dt", 142),
            ("dy_to_s_dt", 143),
            ("dx_to_f_dt", 144),
            ("dy_to_f_dt", 145),
            ("dd_dt", 146),
            ("df_dt", 147),
            ("df_Ca_Bj_dt", 148),
            ("df_Ca_Bsl_dt", 149),
            ("dRy_Rr_dt", 150),
            ("dRy_Ro_dt", 151),
            ("dRy_Ri_dt", 152),
            ("dNa_Bj_dt", 153),
            ("dNa_Bsl_dt", 154),
            ("dTn_CL_dt", 155),
            ("dTn_CHc_dt", 156),
            ("dTn_CHm_dt", 157),
            ("dCaM_dt", 158),
            ("dMyo_c_dt", 159),
            ("dMyo_m_dt", 160),
            ("dSRB_dt", 161),
            ("dSLL_j_dt", 162),
            ("dSLL_sl_dt", 163),
            ("dSLH_j_dt", 164),
            ("dSLH_sl_dt", 165),
            ("dCsqn_b_dt", 166),
            ("dCa_sr_dt", 167),
            ("dNa_j_dt", 168),
            ("dNa_sl_dt", 169),
            ("dNa_i_dt", 170),
            ("dK_i_dt", 171),
            ("dCa_j_dt", 172),
            ("dCa_sl_dt", 173),
            ("dCa_i_dt", 174),
            ("dV_m_dt", 175),
        ]
    )

    indices = []
    for monitor in monitored:
        if monitor not in monitor_inds:
            raise ValueError("Unknown monitored: '{0}'".format(monitor))
        indices.append(monitor_inds[monitor])
    if len(indices) > 1:
        return indices
    else:
        return indices[0]


def rhs(states, t, parameters, values=None):
    """
    Compute the right hand side of the grandi ODE
    """
    # Assign states
    assert len(states) == 39
    (
        m,
        h,
        j,
        x_kr,
        x_ks,
        x_to_s,
        y_to_s,
        x_to_f,
        y_to_f,
        d,
        f,
        f_Ca_Bj,
        f_Ca_Bsl,
        Ry_Rr,
        Ry_Ro,
        Ry_Ri,
        Na_Bj,
        Na_Bsl,
        Tn_CL,
        Tn_CHc,
        Tn_CHm,
        CaM,
        Myo_c,
        Myo_m,
        SRB,
        SLL_j,
        SLL_sl,
        SLH_j,
        SLH_sl,
        Csqn_b,
        Ca_sr,
        Na_j,
        Na_sl,
        Na_i,
        K_i,
        Ca_j,
        Ca_sl,
        Ca_i,
        V_m,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    Fjunc = parameters[0]
    Fjunc_CaL = parameters[1]
    cellLength = parameters[2]
    cellRadius = parameters[3]
    GNa = parameters[8]
    GNaB = parameters[9]
    IbarNaK = parameters[10]
    KmKo = parameters[11]
    KmNaip = parameters[12]
    GKr = parameters[15]
    GKp = parameters[16]
    GKs = parameters[17]
    pNaK = parameters[18]
    GK1 = parameters[19]
    Gto = parameters[20]
    epi = parameters[21]
    GClB = parameters[22]
    GClCa = parameters[23]
    KdClCa = parameters[24]
    GCaL = parameters[25]
    Q10CaL = parameters[26]
    pCa = parameters[27]
    pK = parameters[28]
    pNa = parameters[29]
    IbarNCX = parameters[30]
    Kdact = parameters[31]
    KmCai = parameters[32]
    KmCao = parameters[33]
    KmNai = parameters[34]
    KmNao = parameters[35]
    Q10NCX = parameters[36]
    ksat = parameters[37]
    nu = parameters[38]
    IbarSLCaP = parameters[39]
    KmPCa = parameters[40]
    Q10SLCaP = parameters[41]
    GCaB = parameters[42]
    Kmf = parameters[43]
    Kmr = parameters[44]
    MaxSR = parameters[45]
    MinSR = parameters[46]
    Q10SRCaP = parameters[47]
    Vmax_SRCaP = parameters[48]
    ec50SR = parameters[49]
    hillSRCaP = parameters[50]
    kiCa = parameters[51]
    kim = parameters[52]
    koCa = parameters[53]
    kom = parameters[54]
    ks = parameters[55]
    Bmax_Naj = parameters[56]
    Bmax_Nasl = parameters[57]
    koff_na = parameters[58]
    kon_na = parameters[59]
    Bmax_CaM = parameters[60]
    Bmax_SR = parameters[61]
    Bmax_TnChigh = parameters[62]
    Bmax_TnClow = parameters[63]
    Bmax_myosin = parameters[64]
    koff_cam = parameters[65]
    koff_myoca = parameters[66]
    koff_myomg = parameters[67]
    koff_sr = parameters[68]
    koff_tnchca = parameters[69]
    koff_tnchmg = parameters[70]
    koff_tncl = parameters[71]
    kon_cam = parameters[72]
    kon_myoca = parameters[73]
    kon_myomg = parameters[74]
    kon_sr = parameters[75]
    kon_tnchca = parameters[76]
    kon_tnchmg = parameters[77]
    kon_tncl = parameters[78]
    Bmax_SLhighj0 = parameters[79]
    Bmax_SLhighsl0 = parameters[80]
    Bmax_SLlowj0 = parameters[81]
    Bmax_SLlowsl0 = parameters[82]
    koff_slh = parameters[83]
    koff_sll = parameters[84]
    kon_slh = parameters[85]
    kon_sll = parameters[86]
    Bmax_Csqn0 = parameters[87]
    J_ca_juncsl = parameters[90]
    J_ca_slmyo = parameters[91]
    koff_csqn = parameters[92]
    kon_csqn = parameters[93]
    J_na_juncsl = parameters[96]
    J_na_slmyo = parameters[97]
    Nao = parameters[98]
    Ko = parameters[99]
    Cao = parameters[100]
    Cli = parameters[101]
    Clo = parameters[102]
    Mgi = parameters[103]
    Cmem = parameters[104]
    Frdy = parameters[105]
    R = parameters[106]
    Temp = parameters[107]
    stim_amplitude = parameters[108]
    stim_duration = parameters[109]
    stim_period = parameters[110]
    stim_start = parameters[111]

    # Init return args
    if values is None:
        values = np.zeros_like(states)
    else:
        assert isinstance(values, np.ndarray) and values.shape == states.shape

    # Expressions for the Geometry component
    Vcell = 1e-15 * np.pi * cellLength * (cellRadius * cellRadius)
    Vmyo = 0.65 * Vcell
    Vsr = 0.035 * Vcell
    Vsl = 0.02 * Vcell
    Vjunc = 0.0005390000000000001 * Vcell
    Fsl = 1 - Fjunc
    Fsl_CaL = 1 - Fjunc_CaL

    # Expressions for the Reversal potentials component
    FoRT = Frdy / (R * Temp)
    ena_junc = np.log(Nao / Na_j) / FoRT
    ena_sl = np.log(Nao / Na_sl) / FoRT
    ek = np.log(Ko / K_i) / FoRT
    eca_junc = np.log(Cao / Ca_j) / (2 * FoRT)
    eca_sl = np.log(Cao / Ca_sl) / (2 * FoRT)
    ecl = np.log(Cli / Clo) / FoRT
    Qpow = -31 + Temp / 10

    # Expressions for the I_Na component
    mss = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
    )
    taum = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * V_m)
            * (2.9465894465894467 + 0.06435006435006435 * V_m)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * V_m)
            * (-0.09434663536776214 + 0.019561815336463225 * V_m)
        )
    )
    ah = np.where(
        V_m >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V_m)
    )
    bh = np.where(
        V_m >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V_m)),
        310000.0 * np.exp(0.3485 * V_m) + 2.7 * np.exp(0.079 * V_m),
    )
    tauh = 1.0 / (ah + bh)
    hss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    aj = np.where(
        V_m >= -40,
        0,
        (37.78 + V_m)
        * (-25428.0 * np.exp(0.2444 * V_m) - 6.948e-06 * np.exp(-0.04391 * V_m))
        / (1 + 50262745825.95399 * np.exp(0.311 * V_m)),
    )
    bj = np.where(
        V_m >= -40,
        0.6 * np.exp(0.057 * V_m) / (1 + 0.040762203978366204 * np.exp(-0.1 * V_m)),
        0.02424
        * np.exp(-0.01052 * V_m)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * V_m)),
    )
    tauj = 1.0 / (aj + bj)
    jss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    values[0] = (-m + mss) / taum
    values[1] = (-h + hss) / tauh
    values[2] = (-j + jss) / tauj
    I_Na_junc = Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j
    I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j

    # Expressions for the I_NaBK component
    I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m)
    I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl

    # Expressions for the I_NaK component
    sigma = -1 / 7 + np.exp(0.01485884101040119 * Nao) / 7
    fnak = 1.0 / (
        1 + 0.1245 * np.exp(-0.1 * FoRT * V_m) + 0.0365 * np.exp(-FoRT * V_m) * sigma
    )
    I_nak_junc = (
        Fjunc
        * IbarNaK
        * Ko
        * fnak
        / ((1 + np.power(KmNaip, 4) / np.power(Na_j, 4)) * (KmKo + Ko))
    )
    I_nak_sl = (
        IbarNaK
        * Ko
        * Fsl
        * fnak
        / ((1 + np.power(KmNaip, 4) / np.power(Na_sl, 4)) * (KmKo + Ko))
    )
    I_nak = I_nak_junc + I_nak_sl

    # Expressions for the I_Kr component
    gkr = 0.4303314829119352 * GKr * np.sqrt(Ko)
    xrss = 1.0 / (1 + np.exp(-2 - V_m / 5))
    tauxr = 230 / (1 + np.exp(2 + V_m / 20)) + 3300 / (
        (1 + np.exp(-22 / 9 - V_m / 9)) * (1 + np.exp(11 / 9 + V_m / 9))
    )
    values[3] = (-x_kr + xrss) / tauxr
    rkr = 1.0 / (1 + np.exp(37 / 12 + V_m / 24))
    I_kr = (-ek + V_m) * gkr * rkr * x_kr

    # Expressions for the I_Kp component
    kp_kp = 1.0 / (1 + 1786.4755653786237 * np.exp(-0.16722408026755853 * V_m))
    I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp
    I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp
    I_kp = I_kp_junc + I_kp_sl

    # Expressions for the I_Ks component
    eks = np.log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT
    gks_junc = GKs
    gks_sl = GKs
    xsss = 1.0 / (1 + 0.7659283383646487 * np.exp(-0.07017543859649122 * V_m))
    tauxs = 990.1 / (1 + 0.8415404088681017 * np.exp(-0.0708215297450425 * V_m))
    values[4] = (-x_ks + xsss) / tauxs
    I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc
    I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl
    I_ks = I_ks_junc + I_ks_sl

    # Expressions for the I_to component
    GtoSlow = np.where(epi == 1, 0.12 * Gto, 0.28919999999999996 * Gto)
    GtoFast = np.where(epi == 1, 0.88 * Gto, 0.010799999999999999 * Gto)
    xtoss = 1.0 / (1 + np.exp(19 / 13 - V_m / 13))
    ytoss = 1.0 / (1 + 49.40244910553019 * np.exp(V_m / 5))
    tauxtos = 0.5 + 9 / (1 + np.exp(1 / 5 + V_m / 15))
    tauytos = 30 + 800 / (1 + np.exp(6 + V_m / 10))
    values[5] = (-x_to_s + xtoss) / tauxtos
    values[6] = (-y_to_s + ytoss) / tauytos
    I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s
    tauxtof = 0.5 + 8.5 * np.exp(-((9 / 10 + V_m / 50) * (9 / 10 + V_m / 50)))
    tauytof = 7 + 85 * np.exp(-((40 + V_m) * (40 + V_m)) / 220)
    values[7] = (-x_to_f + xtoss) / tauxtof
    values[8] = (-y_to_f + ytoss) / tauytof
    I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f
    I_to = I_tof + I_tos

    # Expressions for the I_K1 component
    aki = 1.02 / (1 + 7.35454251046446e-07 * np.exp(0.2385 * V_m - 0.2385 * ek))
    bki = (
        0.7626240065063081 * np.exp(0.08032 * V_m - 0.08032 * ek)
        + 1.1534056351865558e-16 * np.exp(0.06175 * V_m - 0.06175 * ek)
    ) / (1 + 0.08677229415769332 * np.exp(0.5143 * ek - 0.5143 * V_m))
    kiss = aki / (aki + bki)
    I_K1 = 0.4303314829119352 * GK1 * np.sqrt(Ko) * (-ek + V_m) * kiss

    # Expressions for the I_ClCa component
    I_ClCa_junc = Fjunc * GClCa * (-ecl + V_m) / (1 + KdClCa / Ca_j)
    I_ClCa_sl = GClCa * (-ecl + V_m) * Fsl / (1 + KdClCa / Ca_sl)
    I_ClCa = I_ClCa_junc + I_ClCa_sl
    I_Clbk = GClB * (-ecl + V_m)

    # Expressions for the I_Ca component
    fss = 1.0 / (1 + np.exp(35 / 9 + V_m / 9)) + 0.6 / (1 + np.exp(5 / 2 - V_m / 20))
    dss = 1.0 / (1 + np.exp(-5 / 6 - V_m / 6))
    taud = (1 - np.exp(-5 / 6 - V_m / 6)) * dss / (0.17500000000000002 + 0.035 * V_m)
    tauf = 1.0 / (
        0.02
        + 0.0197
        * np.exp(
            -(
                (0.48865000000000003 + 0.0337 * V_m)
                * (0.48865000000000003 + 0.0337 * V_m)
            )
        )
    )
    values[9] = (-d + dss) / taud
    values[10] = (-f + fss) / tauf
    values[11] = -0.0119 * f_Ca_Bj + 1.7 * (1 - f_Ca_Bj) * Ca_j
    values[12] = -0.0119 * f_Ca_Bsl + 1.7 * (1 - f_Ca_Bsl) * Ca_sl
    fcaCaMSL = 0
    fcaCaj = 0
    ibarca_j = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_j * np.exp(2 * FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    ibarca_sl = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_sl * np.exp(2 * FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    ibark = (
        Frdy
        * GCaL
        * pK
        * (-0.75 * Ko + 0.75 * K_i * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    ibarna_j = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_j * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    ibarna_sl = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_sl * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    I_Ca_junc = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaj - f_Ca_Bj)
        * d
        * f
        * ibarca_j
    )
    I_Ca_sl = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaMSL - f_Ca_Bsl)
        * Fsl_CaL
        * d
        * f
        * ibarca_sl
    )
    I_CaK = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (Fjunc_CaL * (1 + fcaCaj - f_Ca_Bj) + (1 + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
        * d
        * f
        * ibark
    )
    I_CaNa_junc = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaj - f_Ca_Bj)
        * d
        * f
        * ibarna_j
    )
    I_CaNa_sl = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaMSL - f_Ca_Bsl)
        * Fsl_CaL
        * d
        * f
        * ibarna_sl
    )

    # Expressions for the I_NCX component
    Ka_junc = 1.0 / (1 + (Kdact * Kdact) / (Ca_j * Ca_j))
    Ka_sl = 1.0 / (1 + (Kdact * Kdact) / (Ca_sl * Ca_sl))
    s1_junc = Cao * (Na_j * Na_j * Na_j) * np.exp(nu * FoRT * V_m)
    s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * np.exp(nu * FoRT * V_m)
    s2_junc = (Nao * Nao * Nao) * Ca_j * np.exp((-1 + nu) * FoRT * V_m)
    s3_junc = (
        Cao * (Na_j * Na_j * Na_j)
        + KmCao * (Na_j * Na_j * Na_j)
        + (Nao * Nao * Nao) * Ca_j
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_j / KmCai) * Ca_j
    )
    s2_sl = (Nao * Nao * Nao) * Ca_sl * np.exp((-1 + nu) * FoRT * V_m)
    s3_sl = (
        Cao * (Na_sl * Na_sl * Na_sl)
        + KmCao * (Na_sl * Na_sl * Na_sl)
        + (Nao * Nao * Nao) * Ca_sl
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_sl / KmCai) * Ca_sl
    )
    I_ncx_junc = (
        Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_junc + s1_junc)
        * Ka_junc
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc)
    )
    I_ncx_sl = (
        IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_sl + s1_sl)
        * Fsl
        * Ka_sl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl)
    )

    # Expressions for the I_PCa component
    I_pca_junc = (
        Fjunc
        * IbarSLCaP
        * np.power(Q10SLCaP, Qpow)
        * np.power(Ca_j, 1.6)
        / (np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6))
    )
    I_pca_sl = (
        IbarSLCaP
        * np.power(Q10SLCaP, Qpow)
        * np.power(Ca_sl, 1.6)
        * Fsl
        / (np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6))
    )

    # Expressions for the I_CaBK component
    I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m)
    I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl

    # Expressions for the SR Fluxes component
    kCaSR = MaxSR - (MaxSR - MinSR) / (1 + np.power(ec50SR / Ca_sr, 2.5))
    koSRCa = koCa / kCaSR
    kiSRCa = kiCa * kCaSR
    RI = 1 - Ry_Ri - Ry_Ro - Ry_Rr
    values[13] = (
        kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa
    )
    values[14] = (
        kim * Ry_Ri
        - kom * Ry_Ro
        + (Ca_j * Ca_j) * Ry_Rr * koSRCa
        - Ca_j * Ry_Ro * kiSRCa
    )
    values[15] = (
        -kim * Ry_Ri - kom * Ry_Ri + (Ca_j * Ca_j) * RI * koSRCa + Ca_j * Ry_Ro * kiSRCa
    )
    J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro
    J_serca = (
        Vmax_SRCaP
        * np.power(Q10SRCaP, Qpow)
        * (np.power(Ca_i / Kmf, hillSRCaP) - np.power(Ca_sr / Kmr, hillSRCaP))
        / (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
    )
    J_SRleak = 5.348e-06 * Ca_sr - 5.348e-06 * Ca_j

    # Expressions for the Na Buffers component
    values[16] = -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j
    values[17] = -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl

    # Expressions for the Cytosolic Ca Buffers component
    values[18] = -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    values[19] = (
        -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
    )
    values[20] = -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (
        Bmax_TnChigh - Tn_CHc - Tn_CHm
    )
    values[21] = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i
    values[22] = -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
    values[23] = -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
    values[24] = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i
    J_CaB_cytosol = (
        -koff_cam * CaM
        - koff_myoca * Myo_c
        - koff_myomg * Myo_m
        - koff_sr * SRB
        - koff_tnchca * Tn_CHc
        - koff_tnchmg * Tn_CHm
        - koff_tncl * Tn_CL
        + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
        + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
        + kon_cam * (Bmax_CaM - CaM) * Ca_i
        + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
        + kon_sr * (Bmax_SR - SRB) * Ca_i
        + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
        + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    )

    # Expressions for the Junctional and SL Ca Buffers component
    Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl
    Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc
    Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl
    Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc
    values[25] = -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j
    values[26] = -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl
    values[27] = -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
    values[28] = -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
    J_CaB_junction = (
        -koff_slh * SLH_j
        - koff_sll * SLL_j
        + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
        + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j
    )
    J_CaB_sl = (
        -koff_slh * SLH_sl
        - koff_sll * SLL_sl
        + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
        + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl
    )

    # Expressions for the SR Ca Concentrations component
    Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr
    values[29] = -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
    values[30] = (
        -J_SRCarel
        + koff_csqn * Csqn_b
        - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
        - J_SRleak * Vmyo / Vsr
        + J_serca
    )

    # Expressions for the Na Concentrations component
    I_Na_tot_junc = (
        3 * I_nak_junc + 3 * I_ncx_junc + I_CaNa_junc + I_Na_junc + I_nabk_junc
    )
    I_Na_tot_sl = 3 * I_nak_sl + 3 * I_ncx_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl
    values[31] = (
        -values[16]
        + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
        - Cmem * I_Na_tot_junc / (Frdy * Vjunc)
    )
    values[32] = (
        -values[17]
        + J_na_juncsl * (-Na_sl + Na_j) / Vsl
        + J_na_slmyo * (-Na_sl + Na_i) / Vsl
        - Cmem * I_Na_tot_sl / (Frdy * Vsl)
    )
    values[33] = J_na_slmyo * (-Na_i + Na_sl) / Vmyo

    # Expressions for the K Concentration component
    I_K_tot = -2 * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to
    values[34] = 0

    # Expressions for the Ca Concentrations component
    I_Ca_tot_junc = -2 * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc
    I_Ca_tot_sl = -2 * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl
    values[35] = (
        -J_CaB_junction
        + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
        + J_SRCarel * Vsr / Vjunc
        + J_SRleak * Vmyo / Vjunc
        - Cmem * I_Ca_tot_junc / (2 * Frdy * Vjunc)
    )
    values[36] = (
        -J_CaB_sl
        + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
        + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
        - Cmem * I_Ca_tot_sl / (2 * Frdy * Vsl)
    )
    values[37] = (
        -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo
    )

    # Expressions for the Membrane potential component
    i_Stim = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
            t - stim_period * np.floor(t / stim_period) >= stim_start,
        ),
        -stim_amplitude,
        0,
    )
    I_Na_tot = I_Na_tot_junc + I_Na_tot_sl
    I_Cl_tot = I_ClCa + I_Clbk
    I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl
    I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot
    values[38] = -I_tot - i_Stim

    # Return results
    return values


def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the grandi ODE
    """
    # Assign states
    assert len(states) == 39
    (
        m,
        h,
        j,
        x_kr,
        x_ks,
        x_to_s,
        y_to_s,
        x_to_f,
        y_to_f,
        d,
        f,
        f_Ca_Bj,
        f_Ca_Bsl,
        Ry_Rr,
        Ry_Ro,
        Ry_Ri,
        Na_Bj,
        Na_Bsl,
        Tn_CL,
        Tn_CHc,
        Tn_CHm,
        CaM,
        Myo_c,
        Myo_m,
        SRB,
        SLL_j,
        SLL_sl,
        SLH_j,
        SLH_sl,
        Csqn_b,
        Ca_sr,
        Na_j,
        Na_sl,
        Na_i,
        K_i,
        Ca_j,
        Ca_sl,
        Ca_i,
        V_m,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    Fjunc = parameters[0]
    Fjunc_CaL = parameters[1]
    cellLength = parameters[2]
    cellRadius = parameters[3]
    junctionLength = parameters[6]
    junctionRadius = parameters[7]
    GNa = parameters[8]
    GNaB = parameters[9]
    IbarNaK = parameters[10]
    KmKo = parameters[11]
    KmNaip = parameters[12]
    GKr = parameters[15]
    GKp = parameters[16]
    GKs = parameters[17]
    pNaK = parameters[18]
    GK1 = parameters[19]
    Gto = parameters[20]
    epi = parameters[21]
    GClB = parameters[22]
    GClCa = parameters[23]
    KdClCa = parameters[24]
    GCaL = parameters[25]
    Q10CaL = parameters[26]
    pCa = parameters[27]
    pK = parameters[28]
    pNa = parameters[29]
    IbarNCX = parameters[30]
    Kdact = parameters[31]
    KmCai = parameters[32]
    KmCao = parameters[33]
    KmNai = parameters[34]
    KmNao = parameters[35]
    Q10NCX = parameters[36]
    ksat = parameters[37]
    nu = parameters[38]
    IbarSLCaP = parameters[39]
    KmPCa = parameters[40]
    Q10SLCaP = parameters[41]
    GCaB = parameters[42]
    Kmf = parameters[43]
    Kmr = parameters[44]
    MaxSR = parameters[45]
    MinSR = parameters[46]
    Q10SRCaP = parameters[47]
    Vmax_SRCaP = parameters[48]
    ec50SR = parameters[49]
    hillSRCaP = parameters[50]
    kiCa = parameters[51]
    kim = parameters[52]
    koCa = parameters[53]
    kom = parameters[54]
    ks = parameters[55]
    Bmax_Naj = parameters[56]
    Bmax_Nasl = parameters[57]
    koff_na = parameters[58]
    kon_na = parameters[59]
    Bmax_CaM = parameters[60]
    Bmax_SR = parameters[61]
    Bmax_TnChigh = parameters[62]
    Bmax_TnClow = parameters[63]
    Bmax_myosin = parameters[64]
    koff_cam = parameters[65]
    koff_myoca = parameters[66]
    koff_myomg = parameters[67]
    koff_sr = parameters[68]
    koff_tnchca = parameters[69]
    koff_tnchmg = parameters[70]
    koff_tncl = parameters[71]
    kon_cam = parameters[72]
    kon_myoca = parameters[73]
    kon_myomg = parameters[74]
    kon_sr = parameters[75]
    kon_tnchca = parameters[76]
    kon_tnchmg = parameters[77]
    kon_tncl = parameters[78]
    Bmax_SLhighj0 = parameters[79]
    Bmax_SLhighsl0 = parameters[80]
    Bmax_SLlowj0 = parameters[81]
    Bmax_SLlowsl0 = parameters[82]
    koff_slh = parameters[83]
    koff_sll = parameters[84]
    kon_slh = parameters[85]
    kon_sll = parameters[86]
    Bmax_Csqn0 = parameters[87]
    J_ca_juncsl = parameters[90]
    J_ca_slmyo = parameters[91]
    koff_csqn = parameters[92]
    kon_csqn = parameters[93]
    J_na_juncsl = parameters[96]
    J_na_slmyo = parameters[97]
    Nao = parameters[98]
    Ko = parameters[99]
    Cao = parameters[100]
    Cli = parameters[101]
    Clo = parameters[102]
    Mgi = parameters[103]
    Cmem = parameters[104]
    Frdy = parameters[105]
    R = parameters[106]
    Temp = parameters[107]
    stim_amplitude = parameters[108]
    stim_duration = parameters[109]
    stim_period = parameters[110]
    stim_start = parameters[111]

    # Init return args
    if monitored is None:
        monitored = np.zeros((176,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (176,)

    # Expressions for the Geometry component
    monitored[0] = 1e-15 * np.pi * cellLength * (cellRadius * cellRadius)
    monitored[1] = 0.65 * monitored[0]
    monitored[2] = 0.035 * monitored[0]
    monitored[3] = 0.02 * monitored[0]
    monitored[4] = 0.0005390000000000001 * monitored[0]
    monitored[5] = 40300 * np.pi * junctionLength * junctionRadius
    monitored[6] = 2 * np.pi * cellLength * cellRadius
    monitored[7] = 1 - Fjunc
    monitored[8] = 1 - Fjunc_CaL

    # Expressions for the Reversal potentials component
    monitored[125] = Frdy / (R * Temp)
    monitored[126] = np.log(Nao / Na_j) / monitored[125]
    monitored[127] = np.log(Nao / Na_sl) / monitored[125]
    monitored[128] = np.log(Ko / K_i) / monitored[125]
    monitored[129] = np.log(Cao / Ca_j) / (2 * monitored[125])
    monitored[130] = np.log(Cao / Ca_sl) / (2 * monitored[125])
    monitored[131] = np.log(Cli / Clo) / monitored[125]
    monitored[132] = -31 + Temp / 10

    # Expressions for the I_Na component
    monitored[9] = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
    )
    monitored[10] = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * V_m)
            * (2.9465894465894467 + 0.06435006435006435 * V_m)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * V_m)
            * (-0.09434663536776214 + 0.019561815336463225 * V_m)
        )
    )
    monitored[11] = np.where(
        V_m >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V_m)
    )
    monitored[12] = np.where(
        V_m >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V_m)),
        310000.0 * np.exp(0.3485 * V_m) + 2.7 * np.exp(0.079 * V_m),
    )
    monitored[13] = 1.0 / (monitored[11] + monitored[12])
    monitored[14] = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    monitored[15] = np.where(
        V_m >= -40,
        0,
        (37.78 + V_m)
        * (-25428.0 * np.exp(0.2444 * V_m) - 6.948e-06 * np.exp(-0.04391 * V_m))
        / (1 + 50262745825.95399 * np.exp(0.311 * V_m)),
    )
    monitored[16] = np.where(
        V_m >= -40,
        0.6 * np.exp(0.057 * V_m) / (1 + 0.040762203978366204 * np.exp(-0.1 * V_m)),
        0.02424
        * np.exp(-0.01052 * V_m)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * V_m)),
    )
    monitored[17] = 1.0 / (monitored[15] + monitored[16])
    monitored[18] = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    monitored[137] = (-m + monitored[9]) / monitored[10]
    monitored[138] = (-h + monitored[14]) / monitored[13]
    monitored[139] = (-j + monitored[18]) / monitored[17]
    monitored[19] = Fjunc * GNa * (m * m * m) * (-monitored[126] + V_m) * h * j
    monitored[20] = GNa * (m * m * m) * (-monitored[127] + V_m) * h * j * monitored[7]
    monitored[21] = monitored[19] + monitored[20]

    # Expressions for the I_NaBK component
    monitored[22] = Fjunc * GNaB * (-monitored[126] + V_m)
    monitored[23] = GNaB * (-monitored[127] + V_m) * monitored[7]
    monitored[24] = monitored[22] + monitored[23]

    # Expressions for the I_NaK component
    monitored[25] = -1 / 7 + np.exp(0.01485884101040119 * Nao) / 7
    monitored[26] = 1.0 / (
        1
        + 0.1245 * np.exp(-0.1 * V_m * monitored[125])
        + 0.0365 * np.exp(-V_m * monitored[125]) * monitored[25]
    )
    monitored[27] = (
        Fjunc
        * IbarNaK
        * Ko
        * monitored[26]
        / ((1 + np.power(KmNaip, 4) / np.power(Na_j, 4)) * (KmKo + Ko))
    )
    monitored[28] = (
        IbarNaK
        * Ko
        * monitored[26]
        * monitored[7]
        / ((1 + np.power(KmNaip, 4) / np.power(Na_sl, 4)) * (KmKo + Ko))
    )
    monitored[29] = monitored[27] + monitored[28]

    # Expressions for the I_Kr component
    monitored[30] = 0.4303314829119352 * GKr * np.sqrt(Ko)
    monitored[31] = 1.0 / (1 + np.exp(-2 - V_m / 5))
    monitored[32] = 230 / (1 + np.exp(2 + V_m / 20)) + 3300 / (
        (1 + np.exp(-22 / 9 - V_m / 9)) * (1 + np.exp(11 / 9 + V_m / 9))
    )
    monitored[140] = (-x_kr + monitored[31]) / monitored[32]
    monitored[33] = 1.0 / (1 + np.exp(37 / 12 + V_m / 24))
    monitored[34] = (-monitored[128] + V_m) * monitored[30] * monitored[33] * x_kr

    # Expressions for the I_Kp component
    monitored[35] = 1.0 / (1 + 1786.4755653786237 * np.exp(-0.16722408026755853 * V_m))
    monitored[36] = Fjunc * GKp * (-monitored[128] + V_m) * monitored[35]
    monitored[37] = GKp * (-monitored[128] + V_m) * monitored[35] * monitored[7]
    monitored[38] = monitored[36] + monitored[37]

    # Expressions for the I_Ks component
    monitored[39] = np.log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / monitored[125]
    monitored[40] = GKs
    monitored[41] = GKs
    monitored[42] = 1.0 / (1 + 0.7659283383646487 * np.exp(-0.07017543859649122 * V_m))
    monitored[43] = 990.1 / (1 + 0.8415404088681017 * np.exp(-0.0708215297450425 * V_m))
    monitored[141] = (-x_ks + monitored[42]) / monitored[43]
    monitored[44] = Fjunc * (x_ks * x_ks) * (-monitored[39] + V_m) * monitored[40]
    monitored[45] = (
        (x_ks * x_ks) * (-monitored[39] + V_m) * monitored[41] * monitored[7]
    )
    monitored[46] = monitored[44] + monitored[45]

    # Expressions for the I_to component
    monitored[47] = np.where(epi == 1, 0.12 * Gto, 0.28919999999999996 * Gto)
    monitored[48] = np.where(epi == 1, 0.88 * Gto, 0.010799999999999999 * Gto)
    monitored[49] = 1.0 / (1 + np.exp(19 / 13 - V_m / 13))
    monitored[50] = 1.0 / (1 + 49.40244910553019 * np.exp(V_m / 5))
    monitored[51] = 0.5 + 9 / (1 + np.exp(1 / 5 + V_m / 15))
    monitored[52] = 30 + 800 / (1 + np.exp(6 + V_m / 10))
    monitored[142] = (-x_to_s + monitored[49]) / monitored[51]
    monitored[143] = (-y_to_s + monitored[50]) / monitored[52]
    monitored[53] = (-monitored[128] + V_m) * monitored[47] * x_to_s * y_to_s
    monitored[54] = 0.5 + 8.5 * np.exp(-((9 / 10 + V_m / 50) * (9 / 10 + V_m / 50)))
    monitored[55] = 7 + 85 * np.exp(-((40 + V_m) * (40 + V_m)) / 220)
    monitored[144] = (-x_to_f + monitored[49]) / monitored[54]
    monitored[145] = (-y_to_f + monitored[50]) / monitored[55]
    monitored[56] = (-monitored[128] + V_m) * monitored[48] * x_to_f * y_to_f
    monitored[57] = monitored[53] + monitored[56]

    # Expressions for the I_K1 component
    monitored[133] = 1.02 / (
        1 + 7.35454251046446e-07 * np.exp(0.2385 * V_m - 0.2385 * monitored[128])
    )
    monitored[134] = (
        0.7626240065063081 * np.exp(0.08032 * V_m - 0.08032 * monitored[128])
        + 1.1534056351865558e-16 * np.exp(0.06175 * V_m - 0.06175 * monitored[128])
    ) / (1 + 0.08677229415769332 * np.exp(0.5143 * monitored[128] - 0.5143 * V_m))
    monitored[135] = monitored[133] / (monitored[133] + monitored[134])
    monitored[136] = (
        0.4303314829119352
        * GK1
        * np.sqrt(Ko)
        * (-monitored[128] + V_m)
        * monitored[135]
    )

    # Expressions for the I_ClCa component
    monitored[58] = Fjunc * GClCa * (-monitored[131] + V_m) / (1 + KdClCa / Ca_j)
    monitored[59] = (
        GClCa * (-monitored[131] + V_m) * monitored[7] / (1 + KdClCa / Ca_sl)
    )
    monitored[60] = monitored[58] + monitored[59]
    monitored[61] = GClB * (-monitored[131] + V_m)

    # Expressions for the I_Ca component
    monitored[62] = 1.0 / (1 + np.exp(35 / 9 + V_m / 9)) + 0.6 / (
        1 + np.exp(5 / 2 - V_m / 20)
    )
    monitored[63] = 1.0 / (1 + np.exp(-5 / 6 - V_m / 6))
    monitored[64] = (
        (1 - np.exp(-5 / 6 - V_m / 6))
        * monitored[63]
        / (0.17500000000000002 + 0.035 * V_m)
    )
    monitored[65] = 1.0 / (
        0.02
        + 0.0197
        * np.exp(
            -(
                (0.48865000000000003 + 0.0337 * V_m)
                * (0.48865000000000003 + 0.0337 * V_m)
            )
        )
    )
    monitored[146] = (-d + monitored[63]) / monitored[64]
    monitored[147] = (-f + monitored[62]) / monitored[65]
    monitored[148] = -0.0119 * f_Ca_Bj + 1.7 * (1 - f_Ca_Bj) * Ca_j
    monitored[149] = -0.0119 * f_Ca_Bsl + 1.7 * (1 - f_Ca_Bsl) * Ca_sl
    monitored[66] = 0
    monitored[67] = 0
    monitored[68] = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_j * np.exp(2 * V_m * monitored[125]))
        * V_m
        * monitored[125]
        / (-1 + np.exp(2 * V_m * monitored[125]))
    )
    monitored[69] = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_sl * np.exp(2 * V_m * monitored[125]))
        * V_m
        * monitored[125]
        / (-1 + np.exp(2 * V_m * monitored[125]))
    )
    monitored[70] = (
        Frdy
        * GCaL
        * pK
        * (-0.75 * Ko + 0.75 * K_i * np.exp(V_m * monitored[125]))
        * V_m
        * monitored[125]
        / (-1 + np.exp(V_m * monitored[125]))
    )
    monitored[71] = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_j * np.exp(V_m * monitored[125]))
        * V_m
        * monitored[125]
        / (-1 + np.exp(V_m * monitored[125]))
    )
    monitored[72] = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_sl * np.exp(V_m * monitored[125]))
        * V_m
        * monitored[125]
        / (-1 + np.exp(V_m * monitored[125]))
    )
    monitored[73] = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, monitored[132])
        * (1 + monitored[67] - f_Ca_Bj)
        * d
        * f
        * monitored[68]
    )
    monitored[74] = (
        0.45
        * np.power(Q10CaL, monitored[132])
        * (1 + monitored[66] - f_Ca_Bsl)
        * d
        * f
        * monitored[69]
        * monitored[8]
    )
    monitored[75] = monitored[73] + monitored[74]
    monitored[76] = (
        0.45
        * np.power(Q10CaL, monitored[132])
        * (
            Fjunc_CaL * (1 + monitored[67] - f_Ca_Bj)
            + (1 + monitored[66] - f_Ca_Bsl) * monitored[8]
        )
        * d
        * f
        * monitored[70]
    )
    monitored[77] = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, monitored[132])
        * (1 + monitored[67] - f_Ca_Bj)
        * d
        * f
        * monitored[71]
    )
    monitored[78] = (
        0.45
        * np.power(Q10CaL, monitored[132])
        * (1 + monitored[66] - f_Ca_Bsl)
        * d
        * f
        * monitored[72]
        * monitored[8]
    )
    monitored[79] = monitored[77] + monitored[78]
    monitored[80] = monitored[75] + monitored[76] + monitored[79]

    # Expressions for the I_NCX component
    monitored[81] = 1.0 / (1 + (Kdact * Kdact) / (Ca_j * Ca_j))
    monitored[82] = 1.0 / (1 + (Kdact * Kdact) / (Ca_sl * Ca_sl))
    monitored[83] = Cao * (Na_j * Na_j * Na_j) * np.exp(nu * V_m * monitored[125])
    monitored[84] = Cao * (Na_sl * Na_sl * Na_sl) * np.exp(nu * V_m * monitored[125])
    monitored[85] = (Nao * Nao * Nao) * Ca_j * np.exp((-1 + nu) * V_m * monitored[125])
    monitored[86] = (
        Cao * (Na_j * Na_j * Na_j)
        + KmCao * (Na_j * Na_j * Na_j)
        + (Nao * Nao * Nao) * Ca_j
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_j / KmCai) * Ca_j
    )
    monitored[87] = (Nao * Nao * Nao) * Ca_sl * np.exp((-1 + nu) * V_m * monitored[125])
    monitored[88] = (
        Cao * (Na_sl * Na_sl * Na_sl)
        + KmCao * (Na_sl * Na_sl * Na_sl)
        + (Nao * Nao * Nao) * Ca_sl
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_sl / KmCai) * Ca_sl
    )
    monitored[89] = (
        Fjunc
        * IbarNCX
        * np.power(Q10NCX, monitored[132])
        * (-monitored[85] + monitored[83])
        * monitored[81]
        / ((1 + ksat * np.exp((-1 + nu) * V_m * monitored[125])) * monitored[86])
    )
    monitored[90] = (
        IbarNCX
        * np.power(Q10NCX, monitored[132])
        * (-monitored[87] + monitored[84])
        * monitored[7]
        * monitored[82]
        / ((1 + ksat * np.exp((-1 + nu) * V_m * monitored[125])) * monitored[88])
    )
    monitored[91] = monitored[89] + monitored[90]

    # Expressions for the I_PCa component
    monitored[92] = (
        Fjunc
        * IbarSLCaP
        * np.power(Q10SLCaP, monitored[132])
        * np.power(Ca_j, 1.6)
        / (np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6))
    )
    monitored[93] = (
        IbarSLCaP
        * np.power(Q10SLCaP, monitored[132])
        * np.power(Ca_sl, 1.6)
        * monitored[7]
        / (np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6))
    )
    monitored[94] = monitored[92] + monitored[93]

    # Expressions for the I_CaBK component
    monitored[95] = Fjunc * GCaB * (-monitored[129] + V_m)
    monitored[96] = GCaB * (-monitored[130] + V_m) * monitored[7]
    monitored[97] = monitored[95] + monitored[96]

    # Expressions for the SR Fluxes component
    monitored[98] = MaxSR - (MaxSR - MinSR) / (1 + np.power(ec50SR / Ca_sr, 2.5))
    monitored[99] = koCa / monitored[98]
    monitored[100] = kiCa * monitored[98]
    monitored[101] = 1 - Ry_Ri - Ry_Ro - Ry_Rr
    monitored[150] = (
        kim * monitored[101]
        + kom * Ry_Ro
        - (Ca_j * Ca_j) * Ry_Rr * monitored[99]
        - Ca_j * Ry_Rr * monitored[100]
    )
    monitored[151] = (
        kim * Ry_Ri
        - kom * Ry_Ro
        + (Ca_j * Ca_j) * Ry_Rr * monitored[99]
        - Ca_j * Ry_Ro * monitored[100]
    )
    monitored[152] = (
        -kim * Ry_Ri
        - kom * Ry_Ri
        + (Ca_j * Ca_j) * monitored[101] * monitored[99]
        + Ca_j * Ry_Ro * monitored[100]
    )
    monitored[102] = ks * (-Ca_j + Ca_sr) * Ry_Ro
    monitored[103] = (
        Vmax_SRCaP
        * np.power(Q10SRCaP, monitored[132])
        * (np.power(Ca_i / Kmf, hillSRCaP) - np.power(Ca_sr / Kmr, hillSRCaP))
        / (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
    )
    monitored[104] = 5.348e-06 * Ca_sr - 5.348e-06 * Ca_j

    # Expressions for the Na Buffers component
    monitored[153] = -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j
    monitored[154] = -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl

    # Expressions for the Cytosolic Ca Buffers component
    monitored[155] = -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    monitored[156] = (
        -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
    )
    monitored[157] = -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (
        Bmax_TnChigh - Tn_CHc - Tn_CHm
    )
    monitored[158] = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i
    monitored[159] = (
        -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
    )
    monitored[160] = -koff_myomg * Myo_m + Mgi * kon_myomg * (
        Bmax_myosin - Myo_c - Myo_m
    )
    monitored[161] = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i
    monitored[105] = (
        -koff_cam * CaM
        - koff_myoca * Myo_c
        - koff_myomg * Myo_m
        - koff_sr * SRB
        - koff_tnchca * Tn_CHc
        - koff_tnchmg * Tn_CHm
        - koff_tncl * Tn_CL
        + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
        + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
        + kon_cam * (Bmax_CaM - CaM) * Ca_i
        + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
        + kon_sr * (Bmax_SR - SRB) * Ca_i
        + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
        + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    )

    # Expressions for the Junctional and SL Ca Buffers component
    monitored[106] = Bmax_SLlowsl0 * monitored[1] / monitored[3]
    monitored[107] = Bmax_SLlowj0 * monitored[1] / monitored[4]
    monitored[108] = Bmax_SLhighsl0 * monitored[1] / monitored[3]
    monitored[109] = Bmax_SLhighj0 * monitored[1] / monitored[4]
    monitored[162] = -koff_sll * SLL_j + kon_sll * (-SLL_j + monitored[107]) * Ca_j
    monitored[163] = -koff_sll * SLL_sl + kon_sll * (-SLL_sl + monitored[106]) * Ca_sl
    monitored[164] = -koff_slh * SLH_j + kon_slh * (-SLH_j + monitored[109]) * Ca_j
    monitored[165] = -koff_slh * SLH_sl + kon_slh * (-SLH_sl + monitored[108]) * Ca_sl
    monitored[110] = (
        -koff_slh * SLH_j
        - koff_sll * SLL_j
        + kon_slh * (-SLH_j + monitored[109]) * Ca_j
        + kon_sll * (-SLL_j + monitored[107]) * Ca_j
    )
    monitored[111] = (
        -koff_slh * SLH_sl
        - koff_sll * SLL_sl
        + kon_slh * (-SLH_sl + monitored[108]) * Ca_sl
        + kon_sll * (-SLL_sl + monitored[106]) * Ca_sl
    )

    # Expressions for the SR Ca Concentrations component
    monitored[112] = Bmax_Csqn0 * monitored[1] / monitored[2]
    monitored[166] = -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + monitored[112]) * Ca_sr
    monitored[167] = (
        -monitored[102]
        + koff_csqn * Csqn_b
        - kon_csqn * (-Csqn_b + monitored[112]) * Ca_sr
        - monitored[104] * monitored[1] / monitored[2]
        + monitored[103]
    )

    # Expressions for the Na Concentrations component
    monitored[113] = (
        3 * monitored[27]
        + 3 * monitored[89]
        + monitored[19]
        + monitored[22]
        + monitored[77]
    )
    monitored[114] = (
        3 * monitored[28]
        + 3 * monitored[90]
        + monitored[20]
        + monitored[23]
        + monitored[78]
    )
    monitored[115] = 3 * monitored[28] + 3 * monitored[90] + monitored[78]
    monitored[116] = 3 * monitored[27] + 3 * monitored[89] + monitored[77]
    monitored[168] = (
        -monitored[153]
        + J_na_juncsl * (-Na_j + Na_sl) / monitored[4]
        - Cmem * monitored[113] / (Frdy * monitored[4])
    )
    monitored[169] = (
        -monitored[154]
        + J_na_juncsl * (-Na_sl + Na_j) / monitored[3]
        + J_na_slmyo * (-Na_sl + Na_i) / monitored[3]
        - Cmem * monitored[114] / (Frdy * monitored[3])
    )
    monitored[170] = J_na_slmyo * (-Na_i + Na_sl) / monitored[1]

    # Expressions for the K Concentration component
    monitored[117] = (
        -2 * monitored[29]
        + monitored[136]
        + monitored[34]
        + monitored[38]
        + monitored[46]
        + monitored[57]
        + monitored[76]
    )
    monitored[171] = 0

    # Expressions for the Ca Concentrations component
    monitored[118] = -2 * monitored[89] + monitored[73] + monitored[92] + monitored[95]
    monitored[119] = -2 * monitored[90] + monitored[74] + monitored[93] + monitored[96]
    monitored[172] = (
        -monitored[110]
        + J_ca_juncsl * (-Ca_j + Ca_sl) / monitored[4]
        + monitored[102] * monitored[2] / monitored[4]
        + monitored[104] * monitored[1] / monitored[4]
        - Cmem * monitored[118] / (2 * Frdy * monitored[4])
    )
    monitored[173] = (
        -monitored[111]
        + J_ca_juncsl * (-Ca_sl + Ca_j) / monitored[3]
        + J_ca_slmyo * (-Ca_sl + Ca_i) / monitored[3]
        - Cmem * monitored[119] / (2 * Frdy * monitored[3])
    )
    monitored[174] = (
        -monitored[105]
        + J_ca_slmyo * (-Ca_i + Ca_sl) / monitored[1]
        - monitored[103] * monitored[2] / monitored[1]
    )

    # Expressions for the Membrane potential component
    monitored[120] = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
            t - stim_period * np.floor(t / stim_period) >= stim_start,
        ),
        -stim_amplitude,
        0,
    )
    monitored[121] = monitored[113] + monitored[114]
    monitored[122] = monitored[60] + monitored[61]
    monitored[123] = monitored[118] + monitored[119]
    monitored[124] = monitored[117] + monitored[121] + monitored[122] + monitored[123]
    monitored[175] = -monitored[120] - monitored[124]

    # Return results
    return monitored


def FE(states, t, dt, parameters):
    """
    Compute a forward step using the explicit Euler algorithm to the grandi ODE
    """
    # Assign states
    assert len(states) == 39
    (
        m,
        h,
        j,
        x_kr,
        x_ks,
        x_to_s,
        y_to_s,
        x_to_f,
        y_to_f,
        d,
        f,
        f_Ca_Bj,
        f_Ca_Bsl,
        Ry_Rr,
        Ry_Ro,
        Ry_Ri,
        Na_Bj,
        Na_Bsl,
        Tn_CL,
        Tn_CHc,
        Tn_CHm,
        CaM,
        Myo_c,
        Myo_m,
        SRB,
        SLL_j,
        SLL_sl,
        SLH_j,
        SLH_sl,
        Csqn_b,
        Ca_sr,
        Na_j,
        Na_sl,
        Na_i,
        K_i,
        Ca_j,
        Ca_sl,
        Ca_i,
        V_m,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    Fjunc = parameters[0]
    Fjunc_CaL = parameters[1]
    cellLength = parameters[2]
    cellRadius = parameters[3]
    GNa = parameters[8]
    GNaB = parameters[9]
    IbarNaK = parameters[10]
    KmKo = parameters[11]
    KmNaip = parameters[12]
    GKr = parameters[15]
    GKp = parameters[16]
    GKs = parameters[17]
    pNaK = parameters[18]
    GK1 = parameters[19]
    Gto = parameters[20]
    epi = parameters[21]
    GClB = parameters[22]
    GClCa = parameters[23]
    KdClCa = parameters[24]
    GCaL = parameters[25]
    Q10CaL = parameters[26]
    pCa = parameters[27]
    pK = parameters[28]
    pNa = parameters[29]
    IbarNCX = parameters[30]
    Kdact = parameters[31]
    KmCai = parameters[32]
    KmCao = parameters[33]
    KmNai = parameters[34]
    KmNao = parameters[35]
    Q10NCX = parameters[36]
    ksat = parameters[37]
    nu = parameters[38]
    IbarSLCaP = parameters[39]
    KmPCa = parameters[40]
    Q10SLCaP = parameters[41]
    GCaB = parameters[42]
    Kmf = parameters[43]
    Kmr = parameters[44]
    MaxSR = parameters[45]
    MinSR = parameters[46]
    Q10SRCaP = parameters[47]
    Vmax_SRCaP = parameters[48]
    ec50SR = parameters[49]
    hillSRCaP = parameters[50]
    kiCa = parameters[51]
    kim = parameters[52]
    koCa = parameters[53]
    kom = parameters[54]
    ks = parameters[55]
    Bmax_Naj = parameters[56]
    Bmax_Nasl = parameters[57]
    koff_na = parameters[58]
    kon_na = parameters[59]
    Bmax_CaM = parameters[60]
    Bmax_SR = parameters[61]
    Bmax_TnChigh = parameters[62]
    Bmax_TnClow = parameters[63]
    Bmax_myosin = parameters[64]
    koff_cam = parameters[65]
    koff_myoca = parameters[66]
    koff_myomg = parameters[67]
    koff_sr = parameters[68]
    koff_tnchca = parameters[69]
    koff_tnchmg = parameters[70]
    koff_tncl = parameters[71]
    kon_cam = parameters[72]
    kon_myoca = parameters[73]
    kon_myomg = parameters[74]
    kon_sr = parameters[75]
    kon_tnchca = parameters[76]
    kon_tnchmg = parameters[77]
    kon_tncl = parameters[78]
    Bmax_SLhighj0 = parameters[79]
    Bmax_SLhighsl0 = parameters[80]
    Bmax_SLlowj0 = parameters[81]
    Bmax_SLlowsl0 = parameters[82]
    koff_slh = parameters[83]
    koff_sll = parameters[84]
    kon_slh = parameters[85]
    kon_sll = parameters[86]
    Bmax_Csqn0 = parameters[87]
    J_ca_juncsl = parameters[90]
    J_ca_slmyo = parameters[91]
    koff_csqn = parameters[92]
    kon_csqn = parameters[93]
    J_na_juncsl = parameters[96]
    J_na_slmyo = parameters[97]
    Nao = parameters[98]
    Ko = parameters[99]
    Cao = parameters[100]
    Cli = parameters[101]
    Clo = parameters[102]
    Mgi = parameters[103]
    Cmem = parameters[104]
    Frdy = parameters[105]
    R = parameters[106]
    Temp = parameters[107]
    stim_amplitude = parameters[108]
    stim_duration = parameters[109]
    stim_period = parameters[110]
    stim_start = parameters[111]

    # Expressions for the Geometry component
    Vcell = 1e-15 * np.pi * cellLength * (cellRadius * cellRadius)
    Vmyo = 0.65 * Vcell
    Vsr = 0.035 * Vcell
    Vsl = 0.02 * Vcell
    Vjunc = 0.0005390000000000001 * Vcell
    Fsl = 1 - Fjunc
    Fsl_CaL = 1 - Fjunc_CaL

    # Expressions for the Reversal potentials component
    FoRT = Frdy / (R * Temp)
    ena_junc = np.log(Nao / Na_j) / FoRT
    ena_sl = np.log(Nao / Na_sl) / FoRT
    ek = np.log(Ko / K_i) / FoRT
    eca_junc = np.log(Cao / Ca_j) / (2 * FoRT)
    eca_sl = np.log(Cao / Ca_sl) / (2 * FoRT)
    ecl = np.log(Cli / Clo) / FoRT
    Qpow = -31 + Temp / 10

    # Expressions for the I_Na component
    mss = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
    )
    taum = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * V_m)
            * (2.9465894465894467 + 0.06435006435006435 * V_m)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * V_m)
            * (-0.09434663536776214 + 0.019561815336463225 * V_m)
        )
    )
    ah = np.where(
        V_m >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V_m)
    )
    bh = np.where(
        V_m >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V_m)),
        310000.0 * np.exp(0.3485 * V_m) + 2.7 * np.exp(0.079 * V_m),
    )
    tauh = 1.0 / (ah + bh)
    hss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    aj = np.where(
        V_m >= -40,
        0,
        (37.78 + V_m)
        * (-25428.0 * np.exp(0.2444 * V_m) - 6.948e-06 * np.exp(-0.04391 * V_m))
        / (1 + 50262745825.95399 * np.exp(0.311 * V_m)),
    )
    bj = np.where(
        V_m >= -40,
        0.6 * np.exp(0.057 * V_m) / (1 + 0.040762203978366204 * np.exp(-0.1 * V_m)),
        0.02424
        * np.exp(-0.01052 * V_m)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * V_m)),
    )
    tauj = 1.0 / (aj + bj)
    jss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    dm_dt = (-m + mss) / taum
    states[0] = dt * dm_dt + m
    dh_dt = (-h + hss) / tauh
    states[1] = dt * dh_dt + h
    dj_dt = (-j + jss) / tauj
    states[2] = dt * dj_dt + j
    I_Na_junc = Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j
    I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j

    # Expressions for the I_NaBK component
    I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m)
    I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl

    # Expressions for the I_NaK component
    sigma = -1 / 7 + np.exp(0.01485884101040119 * Nao) / 7
    fnak = 1.0 / (
        1 + 0.1245 * np.exp(-0.1 * FoRT * V_m) + 0.0365 * np.exp(-FoRT * V_m) * sigma
    )
    I_nak_junc = (
        Fjunc
        * IbarNaK
        * Ko
        * fnak
        / ((1 + np.power(KmNaip, 4) / np.power(Na_j, 4)) * (KmKo + Ko))
    )
    I_nak_sl = (
        IbarNaK
        * Ko
        * Fsl
        * fnak
        / ((1 + np.power(KmNaip, 4) / np.power(Na_sl, 4)) * (KmKo + Ko))
    )
    I_nak = I_nak_junc + I_nak_sl

    # Expressions for the I_Kr component
    gkr = 0.4303314829119352 * GKr * np.sqrt(Ko)
    xrss = 1.0 / (1 + np.exp(-2 - V_m / 5))
    tauxr = 230 / (1 + np.exp(2 + V_m / 20)) + 3300 / (
        (1 + np.exp(-22 / 9 - V_m / 9)) * (1 + np.exp(11 / 9 + V_m / 9))
    )
    dx_kr_dt = (-x_kr + xrss) / tauxr
    states[3] = dt * dx_kr_dt + x_kr
    rkr = 1.0 / (1 + np.exp(37 / 12 + V_m / 24))
    I_kr = (-ek + V_m) * gkr * rkr * x_kr

    # Expressions for the I_Kp component
    kp_kp = 1.0 / (1 + 1786.4755653786237 * np.exp(-0.16722408026755853 * V_m))
    I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp
    I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp
    I_kp = I_kp_junc + I_kp_sl

    # Expressions for the I_Ks component
    eks = np.log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT
    gks_junc = GKs
    gks_sl = GKs
    xsss = 1.0 / (1 + 0.7659283383646487 * np.exp(-0.07017543859649122 * V_m))
    tauxs = 990.1 / (1 + 0.8415404088681017 * np.exp(-0.0708215297450425 * V_m))
    dx_ks_dt = (-x_ks + xsss) / tauxs
    states[4] = dt * dx_ks_dt + x_ks
    I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc
    I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl
    I_ks = I_ks_junc + I_ks_sl

    # Expressions for the I_to component
    GtoSlow = np.where(epi == 1, 0.12 * Gto, 0.28919999999999996 * Gto)
    GtoFast = np.where(epi == 1, 0.88 * Gto, 0.010799999999999999 * Gto)
    xtoss = 1.0 / (1 + np.exp(19 / 13 - V_m / 13))
    ytoss = 1.0 / (1 + 49.40244910553019 * np.exp(V_m / 5))
    tauxtos = 0.5 + 9 / (1 + np.exp(1 / 5 + V_m / 15))
    tauytos = 30 + 800 / (1 + np.exp(6 + V_m / 10))
    dx_to_s_dt = (-x_to_s + xtoss) / tauxtos
    states[5] = dt * dx_to_s_dt + x_to_s
    dy_to_s_dt = (-y_to_s + ytoss) / tauytos
    states[6] = dt * dy_to_s_dt + y_to_s
    I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s
    tauxtof = 0.5 + 8.5 * np.exp(-((9 / 10 + V_m / 50) * (9 / 10 + V_m / 50)))
    tauytof = 7 + 85 * np.exp(-((40 + V_m) * (40 + V_m)) / 220)
    dx_to_f_dt = (-x_to_f + xtoss) / tauxtof
    states[7] = dt * dx_to_f_dt + x_to_f
    dy_to_f_dt = (-y_to_f + ytoss) / tauytof
    states[8] = dt * dy_to_f_dt + y_to_f
    I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f
    I_to = I_tof + I_tos

    # Expressions for the I_K1 component
    aki = 1.02 / (1 + 7.35454251046446e-07 * np.exp(0.2385 * V_m - 0.2385 * ek))
    bki = (
        0.7626240065063081 * np.exp(0.08032 * V_m - 0.08032 * ek)
        + 1.1534056351865558e-16 * np.exp(0.06175 * V_m - 0.06175 * ek)
    ) / (1 + 0.08677229415769332 * np.exp(0.5143 * ek - 0.5143 * V_m))
    kiss = aki / (aki + bki)
    I_K1 = 0.4303314829119352 * GK1 * np.sqrt(Ko) * (-ek + V_m) * kiss

    # Expressions for the I_ClCa component
    I_ClCa_junc = Fjunc * GClCa * (-ecl + V_m) / (1 + KdClCa / Ca_j)
    I_ClCa_sl = GClCa * (-ecl + V_m) * Fsl / (1 + KdClCa / Ca_sl)
    I_ClCa = I_ClCa_junc + I_ClCa_sl
    I_Clbk = GClB * (-ecl + V_m)

    # Expressions for the I_Ca component
    fss = 1.0 / (1 + np.exp(35 / 9 + V_m / 9)) + 0.6 / (1 + np.exp(5 / 2 - V_m / 20))
    dss = 1.0 / (1 + np.exp(-5 / 6 - V_m / 6))
    taud = (1 - np.exp(-5 / 6 - V_m / 6)) * dss / (0.17500000000000002 + 0.035 * V_m)
    tauf = 1.0 / (
        0.02
        + 0.0197
        * np.exp(
            -(
                (0.48865000000000003 + 0.0337 * V_m)
                * (0.48865000000000003 + 0.0337 * V_m)
            )
        )
    )
    dd_dt = (-d + dss) / taud
    states[9] = dt * dd_dt + d
    df_dt = (-f + fss) / tauf
    states[10] = dt * df_dt + f
    df_Ca_Bj_dt = -0.0119 * f_Ca_Bj + 1.7 * (1 - f_Ca_Bj) * Ca_j
    states[11] = dt * df_Ca_Bj_dt + f_Ca_Bj
    df_Ca_Bsl_dt = -0.0119 * f_Ca_Bsl + 1.7 * (1 - f_Ca_Bsl) * Ca_sl
    states[12] = dt * df_Ca_Bsl_dt + f_Ca_Bsl
    fcaCaMSL = 0
    fcaCaj = 0
    ibarca_j = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_j * np.exp(2 * FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    ibarca_sl = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_sl * np.exp(2 * FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    ibark = (
        Frdy
        * GCaL
        * pK
        * (-0.75 * Ko + 0.75 * K_i * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    ibarna_j = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_j * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    ibarna_sl = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_sl * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    I_Ca_junc = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaj - f_Ca_Bj)
        * d
        * f
        * ibarca_j
    )
    I_Ca_sl = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaMSL - f_Ca_Bsl)
        * Fsl_CaL
        * d
        * f
        * ibarca_sl
    )
    I_CaK = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (Fjunc_CaL * (1 + fcaCaj - f_Ca_Bj) + (1 + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
        * d
        * f
        * ibark
    )
    I_CaNa_junc = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaj - f_Ca_Bj)
        * d
        * f
        * ibarna_j
    )
    I_CaNa_sl = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaMSL - f_Ca_Bsl)
        * Fsl_CaL
        * d
        * f
        * ibarna_sl
    )

    # Expressions for the I_NCX component
    Ka_junc = 1.0 / (1 + (Kdact * Kdact) / (Ca_j * Ca_j))
    Ka_sl = 1.0 / (1 + (Kdact * Kdact) / (Ca_sl * Ca_sl))
    s1_junc = Cao * (Na_j * Na_j * Na_j) * np.exp(nu * FoRT * V_m)
    s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * np.exp(nu * FoRT * V_m)
    s2_junc = (Nao * Nao * Nao) * Ca_j * np.exp((-1 + nu) * FoRT * V_m)
    s3_junc = (
        Cao * (Na_j * Na_j * Na_j)
        + KmCao * (Na_j * Na_j * Na_j)
        + (Nao * Nao * Nao) * Ca_j
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_j / KmCai) * Ca_j
    )
    s2_sl = (Nao * Nao * Nao) * Ca_sl * np.exp((-1 + nu) * FoRT * V_m)
    s3_sl = (
        Cao * (Na_sl * Na_sl * Na_sl)
        + KmCao * (Na_sl * Na_sl * Na_sl)
        + (Nao * Nao * Nao) * Ca_sl
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_sl / KmCai) * Ca_sl
    )
    I_ncx_junc = (
        Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_junc + s1_junc)
        * Ka_junc
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc)
    )
    I_ncx_sl = (
        IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_sl + s1_sl)
        * Fsl
        * Ka_sl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl)
    )

    # Expressions for the I_PCa component
    I_pca_junc = (
        Fjunc
        * IbarSLCaP
        * np.power(Q10SLCaP, Qpow)
        * np.power(Ca_j, 1.6)
        / (np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6))
    )
    I_pca_sl = (
        IbarSLCaP
        * np.power(Q10SLCaP, Qpow)
        * np.power(Ca_sl, 1.6)
        * Fsl
        / (np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6))
    )

    # Expressions for the I_CaBK component
    I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m)
    I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl

    # Expressions for the SR Fluxes component
    kCaSR = MaxSR - (MaxSR - MinSR) / (1 + np.power(ec50SR / Ca_sr, 2.5))
    koSRCa = koCa / kCaSR
    kiSRCa = kiCa * kCaSR
    RI = 1 - Ry_Ri - Ry_Ro - Ry_Rr
    dRy_Rr_dt = (
        kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa
    )
    states[13] = dt * dRy_Rr_dt + Ry_Rr
    dRy_Ro_dt = (
        kim * Ry_Ri
        - kom * Ry_Ro
        + (Ca_j * Ca_j) * Ry_Rr * koSRCa
        - Ca_j * Ry_Ro * kiSRCa
    )
    states[14] = dt * dRy_Ro_dt + Ry_Ro
    dRy_Ri_dt = (
        -kim * Ry_Ri - kom * Ry_Ri + (Ca_j * Ca_j) * RI * koSRCa + Ca_j * Ry_Ro * kiSRCa
    )
    states[15] = dt * dRy_Ri_dt + Ry_Ri
    J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro
    J_serca = (
        Vmax_SRCaP
        * np.power(Q10SRCaP, Qpow)
        * (np.power(Ca_i / Kmf, hillSRCaP) - np.power(Ca_sr / Kmr, hillSRCaP))
        / (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
    )
    J_SRleak = 5.348e-06 * Ca_sr - 5.348e-06 * Ca_j

    # Expressions for the Na Buffers component
    dNa_Bj_dt = -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j
    states[16] = dt * dNa_Bj_dt + Na_Bj
    dNa_Bsl_dt = -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl
    states[17] = dt * dNa_Bsl_dt + Na_Bsl

    # Expressions for the Cytosolic Ca Buffers component
    dTn_CL_dt = -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    states[18] = dt * dTn_CL_dt + Tn_CL
    dTn_CHc_dt = (
        -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
    )
    states[19] = dt * dTn_CHc_dt + Tn_CHc
    dTn_CHm_dt = -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (
        Bmax_TnChigh - Tn_CHc - Tn_CHm
    )
    states[20] = dt * dTn_CHm_dt + Tn_CHm
    dCaM_dt = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i
    states[21] = dt * dCaM_dt + CaM
    dMyo_c_dt = -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
    states[22] = dt * dMyo_c_dt + Myo_c
    dMyo_m_dt = -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
    states[23] = dt * dMyo_m_dt + Myo_m
    dSRB_dt = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i
    states[24] = dt * dSRB_dt + SRB
    J_CaB_cytosol = (
        -koff_cam * CaM
        - koff_myoca * Myo_c
        - koff_myomg * Myo_m
        - koff_sr * SRB
        - koff_tnchca * Tn_CHc
        - koff_tnchmg * Tn_CHm
        - koff_tncl * Tn_CL
        + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
        + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
        + kon_cam * (Bmax_CaM - CaM) * Ca_i
        + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
        + kon_sr * (Bmax_SR - SRB) * Ca_i
        + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
        + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    )

    # Expressions for the Junctional and SL Ca Buffers component
    Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl
    Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc
    Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl
    Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc
    dSLL_j_dt = -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j
    states[25] = dt * dSLL_j_dt + SLL_j
    dSLL_sl_dt = -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl
    states[26] = dt * dSLL_sl_dt + SLL_sl
    dSLH_j_dt = -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
    states[27] = dt * dSLH_j_dt + SLH_j
    dSLH_sl_dt = -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
    states[28] = dt * dSLH_sl_dt + SLH_sl
    J_CaB_junction = (
        -koff_slh * SLH_j
        - koff_sll * SLL_j
        + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
        + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j
    )
    J_CaB_sl = (
        -koff_slh * SLH_sl
        - koff_sll * SLL_sl
        + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
        + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl
    )

    # Expressions for the SR Ca Concentrations component
    Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr
    dCsqn_b_dt = -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
    states[29] = dt * dCsqn_b_dt + Csqn_b
    dCa_sr_dt = (
        -J_SRCarel
        + koff_csqn * Csqn_b
        - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
        - J_SRleak * Vmyo / Vsr
        + J_serca
    )
    states[30] = dt * dCa_sr_dt + Ca_sr

    # Expressions for the Na Concentrations component
    I_Na_tot_junc = (
        3 * I_nak_junc + 3 * I_ncx_junc + I_CaNa_junc + I_Na_junc + I_nabk_junc
    )
    I_Na_tot_sl = 3 * I_nak_sl + 3 * I_ncx_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl
    dNa_j_dt = (
        -dNa_Bj_dt
        + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
        - Cmem * I_Na_tot_junc / (Frdy * Vjunc)
    )
    states[31] = dt * dNa_j_dt + Na_j
    dNa_sl_dt = (
        -dNa_Bsl_dt
        + J_na_juncsl * (-Na_sl + Na_j) / Vsl
        + J_na_slmyo * (-Na_sl + Na_i) / Vsl
        - Cmem * I_Na_tot_sl / (Frdy * Vsl)
    )
    states[32] = dt * dNa_sl_dt + Na_sl
    dNa_i_dt = J_na_slmyo * (-Na_i + Na_sl) / Vmyo
    states[33] = dt * dNa_i_dt + Na_i

    # Expressions for the K Concentration component
    I_K_tot = -2 * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to
    dK_i_dt = 0
    states[34] = dt * dK_i_dt + K_i

    # Expressions for the Ca Concentrations component
    I_Ca_tot_junc = -2 * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc
    I_Ca_tot_sl = -2 * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl
    dCa_j_dt = (
        -J_CaB_junction
        + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
        + J_SRCarel * Vsr / Vjunc
        + J_SRleak * Vmyo / Vjunc
        - Cmem * I_Ca_tot_junc / (2 * Frdy * Vjunc)
    )
    states[35] = dt * dCa_j_dt + Ca_j
    dCa_sl_dt = (
        -J_CaB_sl
        + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
        + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
        - Cmem * I_Ca_tot_sl / (2 * Frdy * Vsl)
    )
    states[36] = dt * dCa_sl_dt + Ca_sl
    dCa_i_dt = (
        -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo
    )
    states[37] = dt * dCa_i_dt + Ca_i

    # Expressions for the Membrane potential component
    i_Stim = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
            t - stim_period * np.floor(t / stim_period) >= stim_start,
        ),
        -stim_amplitude,
        0,
    )
    I_Na_tot = I_Na_tot_junc + I_Na_tot_sl
    I_Cl_tot = I_ClCa + I_Clbk
    I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl
    I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot
    dV_m_dt = -I_tot - i_Stim
    states[38] = dt * dV_m_dt + V_m

    # Return results
    return states


def GRL1(states, t, dt, parameters):
    """
    Compute a forward step using the rush larsen algorithm to the grandi ODE
    """
    # Assign states
    assert len(states) == 39
    (
        m,
        h,
        j,
        x_kr,
        x_ks,
        x_to_s,
        y_to_s,
        x_to_f,
        y_to_f,
        d,
        f,
        f_Ca_Bj,
        f_Ca_Bsl,
        Ry_Rr,
        Ry_Ro,
        Ry_Ri,
        Na_Bj,
        Na_Bsl,
        Tn_CL,
        Tn_CHc,
        Tn_CHm,
        CaM,
        Myo_c,
        Myo_m,
        SRB,
        SLL_j,
        SLL_sl,
        SLH_j,
        SLH_sl,
        Csqn_b,
        Ca_sr,
        Na_j,
        Na_sl,
        Na_i,
        K_i,
        Ca_j,
        Ca_sl,
        Ca_i,
        V_m,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    Fjunc = parameters[0]
    Fjunc_CaL = parameters[1]
    cellLength = parameters[2]
    cellRadius = parameters[3]
    GNa = parameters[8]
    GNaB = parameters[9]
    IbarNaK = parameters[10]
    KmKo = parameters[11]
    KmNaip = parameters[12]
    GKr = parameters[15]
    GKp = parameters[16]
    GKs = parameters[17]
    pNaK = parameters[18]
    GK1 = parameters[19]
    Gto = parameters[20]
    epi = parameters[21]
    GClB = parameters[22]
    GClCa = parameters[23]
    KdClCa = parameters[24]
    GCaL = parameters[25]
    Q10CaL = parameters[26]
    pCa = parameters[27]
    pK = parameters[28]
    pNa = parameters[29]
    IbarNCX = parameters[30]
    Kdact = parameters[31]
    KmCai = parameters[32]
    KmCao = parameters[33]
    KmNai = parameters[34]
    KmNao = parameters[35]
    Q10NCX = parameters[36]
    ksat = parameters[37]
    nu = parameters[38]
    IbarSLCaP = parameters[39]
    KmPCa = parameters[40]
    Q10SLCaP = parameters[41]
    GCaB = parameters[42]
    Kmf = parameters[43]
    Kmr = parameters[44]
    MaxSR = parameters[45]
    MinSR = parameters[46]
    Q10SRCaP = parameters[47]
    Vmax_SRCaP = parameters[48]
    ec50SR = parameters[49]
    hillSRCaP = parameters[50]
    kiCa = parameters[51]
    kim = parameters[52]
    koCa = parameters[53]
    kom = parameters[54]
    ks = parameters[55]
    Bmax_Naj = parameters[56]
    Bmax_Nasl = parameters[57]
    koff_na = parameters[58]
    kon_na = parameters[59]
    Bmax_CaM = parameters[60]
    Bmax_SR = parameters[61]
    Bmax_TnChigh = parameters[62]
    Bmax_TnClow = parameters[63]
    Bmax_myosin = parameters[64]
    koff_cam = parameters[65]
    koff_myoca = parameters[66]
    koff_myomg = parameters[67]
    koff_sr = parameters[68]
    koff_tnchca = parameters[69]
    koff_tnchmg = parameters[70]
    koff_tncl = parameters[71]
    kon_cam = parameters[72]
    kon_myoca = parameters[73]
    kon_myomg = parameters[74]
    kon_sr = parameters[75]
    kon_tnchca = parameters[76]
    kon_tnchmg = parameters[77]
    kon_tncl = parameters[78]
    Bmax_SLhighj0 = parameters[79]
    Bmax_SLhighsl0 = parameters[80]
    Bmax_SLlowj0 = parameters[81]
    Bmax_SLlowsl0 = parameters[82]
    koff_slh = parameters[83]
    koff_sll = parameters[84]
    kon_slh = parameters[85]
    kon_sll = parameters[86]
    Bmax_Csqn0 = parameters[87]
    J_ca_juncsl = parameters[90]
    J_ca_slmyo = parameters[91]
    koff_csqn = parameters[92]
    kon_csqn = parameters[93]
    J_na_juncsl = parameters[96]
    J_na_slmyo = parameters[97]
    Nao = parameters[98]
    Ko = parameters[99]
    Cao = parameters[100]
    Cli = parameters[101]
    Clo = parameters[102]
    Mgi = parameters[103]
    Cmem = parameters[104]
    Frdy = parameters[105]
    R = parameters[106]
    Temp = parameters[107]
    stim_amplitude = parameters[108]
    stim_duration = parameters[109]
    stim_period = parameters[110]
    stim_start = parameters[111]

    # Expressions for the Geometry component
    Vcell = 1e-15 * np.pi * cellLength * (cellRadius * cellRadius)
    Vmyo = 0.65 * Vcell
    Vsr = 0.035 * Vcell
    Vsl = 0.02 * Vcell
    Vjunc = 0.0005390000000000001 * Vcell
    Fsl = 1 - Fjunc
    Fsl_CaL = 1 - Fjunc_CaL

    # Expressions for the Reversal potentials component
    FoRT = Frdy / (R * Temp)
    ena_junc = np.log(Nao / Na_j) / FoRT
    ena_sl = np.log(Nao / Na_sl) / FoRT
    ek = np.log(Ko / K_i) / FoRT
    eca_junc = np.log(Cao / Ca_j) / (2 * FoRT)
    eca_sl = np.log(Cao / Ca_sl) / (2 * FoRT)
    ecl = np.log(Cli / Clo) / FoRT
    Qpow = -31 + Temp / 10

    # Expressions for the I_Na component
    mss = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V_m))
    )
    taum = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * V_m)
            * (2.9465894465894467 + 0.06435006435006435 * V_m)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * V_m)
            * (-0.09434663536776214 + 0.019561815336463225 * V_m)
        )
    )
    ah = np.where(
        V_m >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V_m)
    )
    bh = np.where(
        V_m >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V_m)),
        310000.0 * np.exp(0.3485 * V_m) + 2.7 * np.exp(0.079 * V_m),
    )
    tauh = 1.0 / (ah + bh)
    hss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    aj = np.where(
        V_m >= -40,
        0,
        (37.78 + V_m)
        * (-25428.0 * np.exp(0.2444 * V_m) - 6.948e-06 * np.exp(-0.04391 * V_m))
        / (1 + 50262745825.95399 * np.exp(0.311 * V_m)),
    )
    bj = np.where(
        V_m >= -40,
        0.6 * np.exp(0.057 * V_m) / (1 + 0.040762203978366204 * np.exp(-0.1 * V_m)),
        0.02424
        * np.exp(-0.01052 * V_m)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * V_m)),
    )
    tauj = 1.0 / (aj + bj)
    jss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V_m))
    )
    dm_dt = (-m + mss) / taum
    dm_dt_linearized = -1 / taum
    states[0] = (
        np.where(
            np.abs(dm_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized,
            dt * dm_dt,
        )
        + m
    )
    dh_dt = (-h + hss) / tauh
    dh_dt_linearized = -1 / tauh
    states[1] = (
        np.where(
            np.abs(dh_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dh_dt_linearized)) * dh_dt / dh_dt_linearized,
            dt * dh_dt,
        )
        + h
    )
    dj_dt = (-j + jss) / tauj
    dj_dt_linearized = -1 / tauj
    states[2] = (
        np.where(
            np.abs(dj_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized,
            dt * dj_dt,
        )
        + j
    )
    I_Na_junc = Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j
    I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j

    # Expressions for the I_NaBK component
    I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m)
    I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl

    # Expressions for the I_NaK component
    sigma = -1 / 7 + np.exp(0.01485884101040119 * Nao) / 7
    fnak = 1.0 / (
        1 + 0.1245 * np.exp(-0.1 * FoRT * V_m) + 0.0365 * np.exp(-FoRT * V_m) * sigma
    )
    I_nak_junc = (
        Fjunc
        * IbarNaK
        * Ko
        * fnak
        / ((1 + np.power(KmNaip, 4) / np.power(Na_j, 4)) * (KmKo + Ko))
    )
    I_nak_sl = (
        IbarNaK
        * Ko
        * Fsl
        * fnak
        / ((1 + np.power(KmNaip, 4) / np.power(Na_sl, 4)) * (KmKo + Ko))
    )
    I_nak = I_nak_junc + I_nak_sl

    # Expressions for the I_Kr component
    gkr = 0.4303314829119352 * GKr * np.sqrt(Ko)
    xrss = 1.0 / (1 + np.exp(-2 - V_m / 5))
    tauxr = 230 / (1 + np.exp(2 + V_m / 20)) + 3300 / (
        (1 + np.exp(-22 / 9 - V_m / 9)) * (1 + np.exp(11 / 9 + V_m / 9))
    )
    dx_kr_dt = (-x_kr + xrss) / tauxr
    dx_kr_dt_linearized = -1 / tauxr
    states[3] = (
        np.where(
            np.abs(dx_kr_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dx_kr_dt_linearized)) * dx_kr_dt / dx_kr_dt_linearized,
            dt * dx_kr_dt,
        )
        + x_kr
    )
    rkr = 1.0 / (1 + np.exp(37 / 12 + V_m / 24))
    I_kr = (-ek + V_m) * gkr * rkr * x_kr

    # Expressions for the I_Kp component
    kp_kp = 1.0 / (1 + 1786.4755653786237 * np.exp(-0.16722408026755853 * V_m))
    I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp
    I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp
    I_kp = I_kp_junc + I_kp_sl

    # Expressions for the I_Ks component
    eks = np.log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT
    gks_junc = GKs
    gks_sl = GKs
    xsss = 1.0 / (1 + 0.7659283383646487 * np.exp(-0.07017543859649122 * V_m))
    tauxs = 990.1 / (1 + 0.8415404088681017 * np.exp(-0.0708215297450425 * V_m))
    dx_ks_dt = (-x_ks + xsss) / tauxs
    dx_ks_dt_linearized = -1 / tauxs
    states[4] = (
        np.where(
            np.abs(dx_ks_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dx_ks_dt_linearized)) * dx_ks_dt / dx_ks_dt_linearized,
            dt * dx_ks_dt,
        )
        + x_ks
    )
    I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc
    I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl
    I_ks = I_ks_junc + I_ks_sl

    # Expressions for the I_to component
    GtoSlow = np.where(epi == 1, 0.12 * Gto, 0.28919999999999996 * Gto)
    GtoFast = np.where(epi == 1, 0.88 * Gto, 0.010799999999999999 * Gto)
    xtoss = 1.0 / (1 + np.exp(19 / 13 - V_m / 13))
    ytoss = 1.0 / (1 + 49.40244910553019 * np.exp(V_m / 5))
    tauxtos = 0.5 + 9 / (1 + np.exp(1 / 5 + V_m / 15))
    tauytos = 30 + 800 / (1 + np.exp(6 + V_m / 10))
    dx_to_s_dt = (-x_to_s + xtoss) / tauxtos
    dx_to_s_dt_linearized = -1 / tauxtos
    states[5] = (
        np.where(
            np.abs(dx_to_s_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dx_to_s_dt_linearized))
            * dx_to_s_dt
            / dx_to_s_dt_linearized,
            dt * dx_to_s_dt,
        )
        + x_to_s
    )
    dy_to_s_dt = (-y_to_s + ytoss) / tauytos
    dy_to_s_dt_linearized = -1 / tauytos
    states[6] = (
        np.where(
            np.abs(dy_to_s_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dy_to_s_dt_linearized))
            * dy_to_s_dt
            / dy_to_s_dt_linearized,
            dt * dy_to_s_dt,
        )
        + y_to_s
    )
    I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s
    tauxtof = 0.5 + 8.5 * np.exp(-((9 / 10 + V_m / 50) * (9 / 10 + V_m / 50)))
    tauytof = 7 + 85 * np.exp(-((40 + V_m) * (40 + V_m)) / 220)
    dx_to_f_dt = (-x_to_f + xtoss) / tauxtof
    dx_to_f_dt_linearized = -1 / tauxtof
    states[7] = (
        np.where(
            np.abs(dx_to_f_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dx_to_f_dt_linearized))
            * dx_to_f_dt
            / dx_to_f_dt_linearized,
            dt * dx_to_f_dt,
        )
        + x_to_f
    )
    dy_to_f_dt = (-y_to_f + ytoss) / tauytof
    dy_to_f_dt_linearized = -1 / tauytof
    states[8] = (
        np.where(
            np.abs(dy_to_f_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dy_to_f_dt_linearized))
            * dy_to_f_dt
            / dy_to_f_dt_linearized,
            dt * dy_to_f_dt,
        )
        + y_to_f
    )
    I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f
    I_to = I_tof + I_tos

    # Expressions for the I_K1 component
    aki = 1.02 / (1 + 7.35454251046446e-07 * np.exp(0.2385 * V_m - 0.2385 * ek))
    bki = (
        0.7626240065063081 * np.exp(0.08032 * V_m - 0.08032 * ek)
        + 1.1534056351865558e-16 * np.exp(0.06175 * V_m - 0.06175 * ek)
    ) / (1 + 0.08677229415769332 * np.exp(0.5143 * ek - 0.5143 * V_m))
    kiss = aki / (aki + bki)
    I_K1 = 0.4303314829119352 * GK1 * np.sqrt(Ko) * (-ek + V_m) * kiss

    # Expressions for the I_ClCa component
    I_ClCa_junc = Fjunc * GClCa * (-ecl + V_m) / (1 + KdClCa / Ca_j)
    I_ClCa_sl = GClCa * (-ecl + V_m) * Fsl / (1 + KdClCa / Ca_sl)
    I_ClCa = I_ClCa_junc + I_ClCa_sl
    I_Clbk = GClB * (-ecl + V_m)

    # Expressions for the I_Ca component
    fss = 1.0 / (1 + np.exp(35 / 9 + V_m / 9)) + 0.6 / (1 + np.exp(5 / 2 - V_m / 20))
    dss = 1.0 / (1 + np.exp(-5 / 6 - V_m / 6))
    taud = (1 - np.exp(-5 / 6 - V_m / 6)) * dss / (0.17500000000000002 + 0.035 * V_m)
    tauf = 1.0 / (
        0.02
        + 0.0197
        * np.exp(
            -(
                (0.48865000000000003 + 0.0337 * V_m)
                * (0.48865000000000003 + 0.0337 * V_m)
            )
        )
    )
    dd_dt = (-d + dss) / taud
    dd_dt_linearized = -1 / taud
    states[9] = (
        np.where(
            np.abs(dd_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dd_dt_linearized)) * dd_dt / dd_dt_linearized,
            dt * dd_dt,
        )
        + d
    )
    df_dt = (-f + fss) / tauf
    df_dt_linearized = -1 / tauf
    states[10] = (
        np.where(
            np.abs(df_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * df_dt_linearized)) * df_dt / df_dt_linearized,
            dt * df_dt,
        )
        + f
    )
    df_Ca_Bj_dt = -0.0119 * f_Ca_Bj + 1.7 * (1 - f_Ca_Bj) * Ca_j
    df_Ca_Bj_dt_linearized = -0.0119 - 1.7 * Ca_j
    states[11] = (
        np.where(
            np.abs(df_Ca_Bj_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * df_Ca_Bj_dt_linearized))
            * df_Ca_Bj_dt
            / df_Ca_Bj_dt_linearized,
            dt * df_Ca_Bj_dt,
        )
        + f_Ca_Bj
    )
    df_Ca_Bsl_dt = -0.0119 * f_Ca_Bsl + 1.7 * (1 - f_Ca_Bsl) * Ca_sl
    df_Ca_Bsl_dt_linearized = -0.0119 - 1.7 * Ca_sl
    states[12] = (
        np.where(
            np.abs(df_Ca_Bsl_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * df_Ca_Bsl_dt_linearized))
            * df_Ca_Bsl_dt
            / df_Ca_Bsl_dt_linearized,
            dt * df_Ca_Bsl_dt,
        )
        + f_Ca_Bsl
    )
    fcaCaMSL = 0
    fcaCaj = 0
    ibarca_j = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_j * np.exp(2 * FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    ibarca_sl = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_sl * np.exp(2 * FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    ibark = (
        Frdy
        * GCaL
        * pK
        * (-0.75 * Ko + 0.75 * K_i * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    ibarna_j = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_j * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    ibarna_sl = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_sl * np.exp(FoRT * V_m))
        * FoRT
        * V_m
        / (-1 + np.exp(FoRT * V_m))
    )
    I_Ca_junc = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaj - f_Ca_Bj)
        * d
        * f
        * ibarca_j
    )
    I_Ca_sl = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaMSL - f_Ca_Bsl)
        * Fsl_CaL
        * d
        * f
        * ibarca_sl
    )
    I_CaK = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (Fjunc_CaL * (1 + fcaCaj - f_Ca_Bj) + (1 + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
        * d
        * f
        * ibark
    )
    I_CaNa_junc = (
        0.45
        * Fjunc_CaL
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaj - f_Ca_Bj)
        * d
        * f
        * ibarna_j
    )
    I_CaNa_sl = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (1 + fcaCaMSL - f_Ca_Bsl)
        * Fsl_CaL
        * d
        * f
        * ibarna_sl
    )

    # Expressions for the I_NCX component
    Ka_junc = 1.0 / (1 + (Kdact * Kdact) / (Ca_j * Ca_j))
    Ka_sl = 1.0 / (1 + (Kdact * Kdact) / (Ca_sl * Ca_sl))
    s1_junc = Cao * (Na_j * Na_j * Na_j) * np.exp(nu * FoRT * V_m)
    s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * np.exp(nu * FoRT * V_m)
    s2_junc = (Nao * Nao * Nao) * Ca_j * np.exp((-1 + nu) * FoRT * V_m)
    s3_junc = (
        Cao * (Na_j * Na_j * Na_j)
        + KmCao * (Na_j * Na_j * Na_j)
        + (Nao * Nao * Nao) * Ca_j
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_j / KmCai) * Ca_j
    )
    s2_sl = (Nao * Nao * Nao) * Ca_sl * np.exp((-1 + nu) * FoRT * V_m)
    s3_sl = (
        Cao * (Na_sl * Na_sl * Na_sl)
        + KmCao * (Na_sl * Na_sl * Na_sl)
        + (Nao * Nao * Nao) * Ca_sl
        + KmCai
        * (Nao * Nao * Nao)
        * (1 + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
        + (KmNao * KmNao * KmNao) * (1 + Ca_sl / KmCai) * Ca_sl
    )
    I_ncx_junc = (
        Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_junc + s1_junc)
        * Ka_junc
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc)
    )
    I_ncx_sl = (
        IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_sl + s1_sl)
        * Fsl
        * Ka_sl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl)
    )

    # Expressions for the I_PCa component
    I_pca_junc = (
        Fjunc
        * IbarSLCaP
        * np.power(Q10SLCaP, Qpow)
        * np.power(Ca_j, 1.6)
        / (np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6))
    )
    I_pca_sl = (
        IbarSLCaP
        * np.power(Q10SLCaP, Qpow)
        * np.power(Ca_sl, 1.6)
        * Fsl
        / (np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6))
    )

    # Expressions for the I_CaBK component
    I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m)
    I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl

    # Expressions for the SR Fluxes component
    kCaSR = MaxSR - (MaxSR - MinSR) / (1 + np.power(ec50SR / Ca_sr, 2.5))
    koSRCa = koCa / kCaSR
    kiSRCa = kiCa * kCaSR
    RI = 1 - Ry_Ri - Ry_Ro - Ry_Rr
    dRy_Rr_dt = (
        kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa
    )
    dRy_Rr_dt_linearized = -kim - (Ca_j * Ca_j) * koSRCa - Ca_j * kiSRCa
    states[13] = (
        np.where(
            np.abs(dRy_Rr_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dRy_Rr_dt_linearized))
            * dRy_Rr_dt
            / dRy_Rr_dt_linearized,
            dt * dRy_Rr_dt,
        )
        + Ry_Rr
    )
    dRy_Ro_dt = (
        kim * Ry_Ri
        - kom * Ry_Ro
        + (Ca_j * Ca_j) * Ry_Rr * koSRCa
        - Ca_j * Ry_Ro * kiSRCa
    )
    dRy_Ro_dt_linearized = -kom - Ca_j * kiSRCa
    states[14] = (
        np.where(
            np.abs(dRy_Ro_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dRy_Ro_dt_linearized))
            * dRy_Ro_dt
            / dRy_Ro_dt_linearized,
            dt * dRy_Ro_dt,
        )
        + Ry_Ro
    )
    dRy_Ri_dt = (
        -kim * Ry_Ri - kom * Ry_Ri + (Ca_j * Ca_j) * RI * koSRCa + Ca_j * Ry_Ro * kiSRCa
    )
    dRy_Ri_dt_linearized = -kim - kom - (Ca_j * Ca_j) * koSRCa
    states[15] = (
        np.where(
            np.abs(dRy_Ri_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dRy_Ri_dt_linearized))
            * dRy_Ri_dt
            / dRy_Ri_dt_linearized,
            dt * dRy_Ri_dt,
        )
        + Ry_Ri
    )
    J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro
    J_serca = (
        Vmax_SRCaP
        * np.power(Q10SRCaP, Qpow)
        * (np.power(Ca_i / Kmf, hillSRCaP) - np.power(Ca_sr / Kmr, hillSRCaP))
        / (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
    )
    J_SRleak = 5.348e-06 * Ca_sr - 5.348e-06 * Ca_j

    # Expressions for the Na Buffers component
    dNa_Bj_dt = -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j
    dNa_Bj_dt_linearized = -koff_na - kon_na * Na_j
    states[16] = Na_Bj + np.where(
        np.abs(dNa_Bj_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dNa_Bj_dt_linearized)) * dNa_Bj_dt / dNa_Bj_dt_linearized,
        dt * dNa_Bj_dt,
    )
    dNa_Bsl_dt = -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl
    dNa_Bsl_dt_linearized = -koff_na - kon_na * Na_sl
    states[17] = Na_Bsl + np.where(
        np.abs(dNa_Bsl_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dNa_Bsl_dt_linearized))
        * dNa_Bsl_dt
        / dNa_Bsl_dt_linearized,
        dt * dNa_Bsl_dt,
    )

    # Expressions for the Cytosolic Ca Buffers component
    dTn_CL_dt = -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    dTn_CL_dt_linearized = -koff_tncl - kon_tncl * Ca_i
    states[18] = (
        np.where(
            np.abs(dTn_CL_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dTn_CL_dt_linearized))
            * dTn_CL_dt
            / dTn_CL_dt_linearized,
            dt * dTn_CL_dt,
        )
        + Tn_CL
    )
    dTn_CHc_dt = (
        -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
    )
    dTn_CHc_dt_linearized = -koff_tnchca - kon_tnchca * Ca_i
    states[19] = (
        np.where(
            np.abs(dTn_CHc_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dTn_CHc_dt_linearized))
            * dTn_CHc_dt
            / dTn_CHc_dt_linearized,
            dt * dTn_CHc_dt,
        )
        + Tn_CHc
    )
    dTn_CHm_dt = -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (
        Bmax_TnChigh - Tn_CHc - Tn_CHm
    )
    dTn_CHm_dt_linearized = -koff_tnchmg - Mgi * kon_tnchmg
    states[20] = (
        np.where(
            np.abs(dTn_CHm_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dTn_CHm_dt_linearized))
            * dTn_CHm_dt
            / dTn_CHm_dt_linearized,
            dt * dTn_CHm_dt,
        )
        + Tn_CHm
    )
    dCaM_dt = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i
    dCaM_dt_linearized = -koff_cam - kon_cam * Ca_i
    states[21] = CaM + np.where(
        np.abs(dCaM_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dCaM_dt_linearized)) * dCaM_dt / dCaM_dt_linearized,
        dt * dCaM_dt,
    )
    dMyo_c_dt = -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
    dMyo_c_dt_linearized = -koff_myoca - kon_myoca * Ca_i
    states[22] = Myo_c + np.where(
        np.abs(dMyo_c_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dMyo_c_dt_linearized)) * dMyo_c_dt / dMyo_c_dt_linearized,
        dt * dMyo_c_dt,
    )
    dMyo_m_dt = -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
    dMyo_m_dt_linearized = -koff_myomg - Mgi * kon_myomg
    states[23] = Myo_m + np.where(
        np.abs(dMyo_m_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dMyo_m_dt_linearized)) * dMyo_m_dt / dMyo_m_dt_linearized,
        dt * dMyo_m_dt,
    )
    dSRB_dt = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i
    dSRB_dt_linearized = -koff_sr - kon_sr * Ca_i
    states[24] = (
        np.where(
            np.abs(dSRB_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dSRB_dt_linearized)) * dSRB_dt / dSRB_dt_linearized,
            dt * dSRB_dt,
        )
        + SRB
    )
    J_CaB_cytosol = (
        -koff_cam * CaM
        - koff_myoca * Myo_c
        - koff_myomg * Myo_m
        - koff_sr * SRB
        - koff_tnchca * Tn_CHc
        - koff_tnchmg * Tn_CHm
        - koff_tncl * Tn_CL
        + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
        + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
        + kon_cam * (Bmax_CaM - CaM) * Ca_i
        + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
        + kon_sr * (Bmax_SR - SRB) * Ca_i
        + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
        + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i
    )

    # Expressions for the Junctional and SL Ca Buffers component
    Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl
    Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc
    Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl
    Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc
    dSLL_j_dt = -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j
    dSLL_j_dt_linearized = -koff_sll - kon_sll * Ca_j
    states[25] = (
        np.where(
            np.abs(dSLL_j_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dSLL_j_dt_linearized))
            * dSLL_j_dt
            / dSLL_j_dt_linearized,
            dt * dSLL_j_dt,
        )
        + SLL_j
    )
    dSLL_sl_dt = -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl
    dSLL_sl_dt_linearized = -koff_sll - kon_sll * Ca_sl
    states[26] = (
        np.where(
            np.abs(dSLL_sl_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dSLL_sl_dt_linearized))
            * dSLL_sl_dt
            / dSLL_sl_dt_linearized,
            dt * dSLL_sl_dt,
        )
        + SLL_sl
    )
    dSLH_j_dt = -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
    dSLH_j_dt_linearized = -koff_slh - kon_slh * Ca_j
    states[27] = (
        np.where(
            np.abs(dSLH_j_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dSLH_j_dt_linearized))
            * dSLH_j_dt
            / dSLH_j_dt_linearized,
            dt * dSLH_j_dt,
        )
        + SLH_j
    )
    dSLH_sl_dt = -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
    dSLH_sl_dt_linearized = -koff_slh - kon_slh * Ca_sl
    states[28] = (
        np.where(
            np.abs(dSLH_sl_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dSLH_sl_dt_linearized))
            * dSLH_sl_dt
            / dSLH_sl_dt_linearized,
            dt * dSLH_sl_dt,
        )
        + SLH_sl
    )
    J_CaB_junction = (
        -koff_slh * SLH_j
        - koff_sll * SLL_j
        + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
        + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j
    )
    J_CaB_sl = (
        -koff_slh * SLH_sl
        - koff_sll * SLL_sl
        + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
        + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl
    )

    # Expressions for the SR Ca Concentrations component
    Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr
    dCsqn_b_dt = -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
    dCsqn_b_dt_linearized = -koff_csqn - kon_csqn * Ca_sr
    states[29] = Csqn_b + np.where(
        np.abs(dCsqn_b_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dCsqn_b_dt_linearized))
        * dCsqn_b_dt
        / dCsqn_b_dt_linearized,
        dt * dCsqn_b_dt,
    )
    dCa_sr_dt = (
        -J_SRCarel
        + koff_csqn * Csqn_b
        - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
        - J_SRleak * Vmyo / Vsr
        + J_serca
    )
    dJ_serca_dCa_sr = -Vmax_SRCaP * hillSRCaP * np.power(Q10SRCaP, Qpow) * np.power(
        Ca_sr / Kmr, hillSRCaP
    ) / (
        (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP)) * Ca_sr
    ) - Vmax_SRCaP * hillSRCaP * np.power(
        Q10SRCaP, Qpow
    ) * np.power(
        Ca_sr / Kmr, hillSRCaP
    ) * (
        np.power(Ca_i / Kmf, hillSRCaP) - np.power(Ca_sr / Kmr, hillSRCaP)
    ) / (
        (
            (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
            * (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
        )
        * Ca_sr
    )
    dCa_sr_dt_linearized = (
        -kon_csqn * (-Csqn_b + Bmax_Csqn)
        - ks * Ry_Ro
        - 5.348e-06 * Vmyo / Vsr
        + dJ_serca_dCa_sr
    )
    states[30] = Ca_sr + np.where(
        np.abs(dCa_sr_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dCa_sr_dt_linearized)) * dCa_sr_dt / dCa_sr_dt_linearized,
        dt * dCa_sr_dt,
    )

    # Expressions for the Na Concentrations component
    I_Na_tot_junc = (
        3 * I_nak_junc + 3 * I_ncx_junc + I_CaNa_junc + I_Na_junc + I_nabk_junc
    )
    I_Na_tot_sl = 3 * I_nak_sl + 3 * I_ncx_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl
    dNa_j_dt = (
        -dNa_Bj_dt
        + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
        - Cmem * I_Na_tot_junc / (Frdy * Vjunc)
    )
    dI_CaNa_junc_dibarna_j = (
        0.45 * Fjunc_CaL * np.power(Q10CaL, Qpow) * (1 + fcaCaj - f_Ca_Bj) * d * f
    )
    ds3_junc_dNa_j = (
        3 * Cao * (Na_j * Na_j)
        + 3 * KmCao * (Na_j * Na_j)
        + 3 * KmCai * (Nao * Nao * Nao) * (Na_j * Na_j) / (KmNai * KmNai * KmNai)
    )
    dena_junc_dNa_j = -1 / (FoRT * Na_j)
    ds1_junc_dNa_j = 3 * Cao * (Na_j * Na_j) * np.exp(nu * FoRT * V_m)
    dI_nabk_junc_dena_junc = -Fjunc * GNaB
    dibarna_j_dNa_j = (
        0.75
        * Frdy
        * GCaL
        * pNa
        * FoRT
        * V_m
        * np.exp(FoRT * V_m)
        / (-1 + np.exp(FoRT * V_m))
    )
    dI_Na_junc_dena_junc = -Fjunc * GNa * (m * m * m) * h * j
    dI_ncx_junc_ds1_junc = (
        Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * Ka_junc
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc)
    )
    dI_nak_junc_dNa_j = (
        4
        * Fjunc
        * IbarNaK
        * Ko
        * np.power(KmNaip, 4)
        * fnak
        / (
            (
                (1 + np.power(KmNaip, 4) / np.power(Na_j, 4))
                * (1 + np.power(KmNaip, 4) / np.power(Na_j, 4))
            )
            * (KmKo + Ko)
            * np.power(Na_j, 5)
        )
    )
    dI_ncx_junc_ds3_junc = (
        -Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_junc + s1_junc)
        * Ka_junc
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * (s3_junc * s3_junc))
    )
    dNa_j_dt_linearized = -J_na_juncsl / Vjunc - Cmem * (
        3 * dI_nak_junc_dNa_j
        + dI_CaNa_junc_dibarna_j * dibarna_j_dNa_j
        + dI_Na_junc_dena_junc * dena_junc_dNa_j
        + dI_nabk_junc_dena_junc * dena_junc_dNa_j
        + 3 * dI_ncx_junc_ds1_junc * ds1_junc_dNa_j
        + 3 * dI_ncx_junc_ds3_junc * ds3_junc_dNa_j
    ) / (Frdy * Vjunc)
    states[31] = Na_j + np.where(
        np.abs(dNa_j_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dNa_j_dt_linearized)) * dNa_j_dt / dNa_j_dt_linearized,
        dt * dNa_j_dt,
    )
    dNa_sl_dt = (
        -dNa_Bsl_dt
        + J_na_juncsl * (-Na_sl + Na_j) / Vsl
        + J_na_slmyo * (-Na_sl + Na_i) / Vsl
        - Cmem * I_Na_tot_sl / (Frdy * Vsl)
    )
    dI_nabk_sl_dena_sl = -GNaB * Fsl
    dI_CaNa_sl_dibarna_sl = (
        0.45 * np.power(Q10CaL, Qpow) * (1 + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f
    )
    dibarna_sl_dNa_sl = (
        0.75
        * Frdy
        * GCaL
        * pNa
        * FoRT
        * V_m
        * np.exp(FoRT * V_m)
        / (-1 + np.exp(FoRT * V_m))
    )
    dI_Na_sl_dena_sl = -GNa * (m * m * m) * Fsl * h * j
    dI_ncx_sl_ds1_sl = (
        IbarNCX
        * np.power(Q10NCX, Qpow)
        * Fsl
        * Ka_sl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl)
    )
    ds3_sl_dNa_sl = (
        3 * Cao * (Na_sl * Na_sl)
        + 3 * KmCao * (Na_sl * Na_sl)
        + 3 * KmCai * (Nao * Nao * Nao) * (Na_sl * Na_sl) / (KmNai * KmNai * KmNai)
    )
    dI_nak_sl_dNa_sl = (
        4
        * IbarNaK
        * Ko
        * np.power(KmNaip, 4)
        * Fsl
        * fnak
        / (
            (
                (1 + np.power(KmNaip, 4) / np.power(Na_sl, 4))
                * (1 + np.power(KmNaip, 4) / np.power(Na_sl, 4))
            )
            * (KmKo + Ko)
            * np.power(Na_sl, 5)
        )
    )
    dI_ncx_sl_ds3_sl = (
        -IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_sl + s1_sl)
        * Fsl
        * Ka_sl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * (s3_sl * s3_sl))
    )
    dena_sl_dNa_sl = -1 / (FoRT * Na_sl)
    ds1_sl_dNa_sl = 3 * Cao * (Na_sl * Na_sl) * np.exp(nu * FoRT * V_m)
    dNa_sl_dt_linearized = (
        -J_na_juncsl / Vsl
        - J_na_slmyo / Vsl
        - Cmem
        * (
            3 * dI_nak_sl_dNa_sl
            + dI_CaNa_sl_dibarna_sl * dibarna_sl_dNa_sl
            + dI_Na_sl_dena_sl * dena_sl_dNa_sl
            + dI_nabk_sl_dena_sl * dena_sl_dNa_sl
            + 3 * dI_ncx_sl_ds1_sl * ds1_sl_dNa_sl
            + 3 * dI_ncx_sl_ds3_sl * ds3_sl_dNa_sl
        )
        / (Frdy * Vsl)
    )
    states[32] = Na_sl + np.where(
        np.abs(dNa_sl_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dNa_sl_dt_linearized)) * dNa_sl_dt / dNa_sl_dt_linearized,
        dt * dNa_sl_dt,
    )
    dNa_i_dt = J_na_slmyo * (-Na_i + Na_sl) / Vmyo
    dNa_i_dt_linearized = -J_na_slmyo / Vmyo
    states[33] = Na_i + np.where(
        np.abs(dNa_i_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dNa_i_dt_linearized)) * dNa_i_dt / dNa_i_dt_linearized,
        dt * dNa_i_dt,
    )

    # Expressions for the K Concentration component
    I_K_tot = -2 * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to
    dK_i_dt = 0
    states[34] = dt * dK_i_dt + K_i

    # Expressions for the Ca Concentrations component
    I_Ca_tot_junc = -2 * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc
    I_Ca_tot_sl = -2 * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl
    dCa_j_dt = (
        -J_CaB_junction
        + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
        + J_SRCarel * Vsr / Vjunc
        + J_SRleak * Vmyo / Vjunc
        - Cmem * I_Ca_tot_junc / (2 * Frdy * Vjunc)
    )
    dI_Ca_junc_dibarca_j = (
        0.45 * Fjunc_CaL * np.power(Q10CaL, Qpow) * (1 + fcaCaj - f_Ca_Bj) * d * f
    )
    dI_cabk_junc_deca_junc = -Fjunc * GCaB
    dibarca_j_dCa_j = (
        1.364
        * Frdy
        * GCaL
        * pCa
        * FoRT
        * V_m
        * np.exp(2 * FoRT * V_m)
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    dKa_junc_dCa_j = (
        2
        * (Kdact * Kdact)
        / (
            (
                (1 + (Kdact * Kdact) / (Ca_j * Ca_j))
                * (1 + (Kdact * Kdact) / (Ca_j * Ca_j))
            )
            * (Ca_j * Ca_j * Ca_j)
        )
    )
    dJ_CaB_junction_dCa_j = kon_slh * (-SLH_j + Bmax_SLhighj) + kon_sll * (
        -SLL_j + Bmax_SLlowj
    )
    dI_ncx_junc_ds2_junc = (
        -Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * Ka_junc
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc)
    )
    ds2_junc_dCa_j = (Nao * Nao * Nao) * np.exp((-1 + nu) * FoRT * V_m)
    ds3_junc_dCa_j = (
        (Nao * Nao * Nao)
        + (KmNao * KmNao * KmNao) * (1 + Ca_j / KmCai)
        + (KmNao * KmNao * KmNao) * Ca_j / KmCai
    )
    dJ_SRCarel_dCa_j = -ks * Ry_Ro
    dI_pca_junc_dCa_j = 1.6 * Fjunc * IbarSLCaP * np.power(Q10SLCaP, Qpow) * np.power(
        Ca_j, 0.6000000000000001
    ) / (
        np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6)
    ) - 1.6 * Fjunc * IbarSLCaP * np.power(
        Q10SLCaP, Qpow
    ) * np.power(
        Ca_j, 2.2
    ) / (
        (np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6))
        * (np.power(KmPCa, 1.6) + np.power(Ca_j, 1.6))
    )
    dI_ncx_junc_dKa_junc = (
        Fjunc
        * IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_junc + s1_junc)
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc)
    )
    deca_junc_dCa_j = -1 / (2 * Ca_j * FoRT)
    dCa_j_dt_linearized = (
        -dJ_CaB_junction_dCa_j
        - J_ca_juncsl / Vjunc
        - 5.348e-06 * Vmyo / Vjunc
        + Vsr * dJ_SRCarel_dCa_j / Vjunc
        - Cmem
        * (
            dI_Ca_junc_dibarca_j * dibarca_j_dCa_j
            + dI_cabk_junc_deca_junc * deca_junc_dCa_j
            - 2 * dI_ncx_junc_dKa_junc * dKa_junc_dCa_j
            - 2 * dI_ncx_junc_ds2_junc * ds2_junc_dCa_j
            - 2 * dI_ncx_junc_ds3_junc * ds3_junc_dCa_j
            + dI_pca_junc_dCa_j
        )
        / (2 * Frdy * Vjunc)
    )
    states[35] = Ca_j + np.where(
        np.abs(dCa_j_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dCa_j_dt_linearized)) * dCa_j_dt / dCa_j_dt_linearized,
        dt * dCa_j_dt,
    )
    dCa_sl_dt = (
        -J_CaB_sl
        + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
        + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
        - Cmem * I_Ca_tot_sl / (2 * Frdy * Vsl)
    )
    dI_Ca_sl_dibarca_sl = (
        0.45 * np.power(Q10CaL, Qpow) * (1 + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f
    )
    ds2_sl_dCa_sl = (Nao * Nao * Nao) * np.exp((-1 + nu) * FoRT * V_m)
    dKa_sl_dCa_sl = (
        2
        * (Kdact * Kdact)
        / (
            (
                (1 + (Kdact * Kdact) / (Ca_sl * Ca_sl))
                * (1 + (Kdact * Kdact) / (Ca_sl * Ca_sl))
            )
            * (Ca_sl * Ca_sl * Ca_sl)
        )
    )
    ds3_sl_dCa_sl = (
        (Nao * Nao * Nao)
        + (KmNao * KmNao * KmNao) * (1 + Ca_sl / KmCai)
        + (KmNao * KmNao * KmNao) * Ca_sl / KmCai
    )
    dI_ncx_sl_dKa_sl = (
        IbarNCX
        * np.power(Q10NCX, Qpow)
        * (-s2_sl + s1_sl)
        * Fsl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl)
    )
    dibarca_sl_dCa_sl = (
        1.364
        * Frdy
        * GCaL
        * pCa
        * FoRT
        * V_m
        * np.exp(2 * FoRT * V_m)
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    dI_pca_sl_dCa_sl = 1.6 * IbarSLCaP * np.power(Q10SLCaP, Qpow) * np.power(
        Ca_sl, 0.6000000000000001
    ) * Fsl / (
        np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6)
    ) - 1.6 * IbarSLCaP * np.power(
        Q10SLCaP, Qpow
    ) * np.power(
        Ca_sl, 2.2
    ) * Fsl / (
        (np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6))
        * (np.power(KmPCa, 1.6) + np.power(Ca_sl, 1.6))
    )
    dJ_CaB_sl_dCa_sl = kon_slh * (-SLH_sl + Bmax_SLhighsl) + kon_sll * (
        -SLL_sl + Bmax_SLlowsl
    )
    deca_sl_dCa_sl = -1 / (2 * Ca_sl * FoRT)
    dI_ncx_sl_ds2_sl = (
        -IbarNCX
        * np.power(Q10NCX, Qpow)
        * Fsl
        * Ka_sl
        / ((1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl)
    )
    dI_cabk_sl_deca_sl = -GCaB * Fsl
    dCa_sl_dt_linearized = (
        -dJ_CaB_sl_dCa_sl
        - J_ca_juncsl / Vsl
        - J_ca_slmyo / Vsl
        - Cmem
        * (
            dI_Ca_sl_dibarca_sl * dibarca_sl_dCa_sl
            + dI_cabk_sl_deca_sl * deca_sl_dCa_sl
            - 2 * dI_ncx_sl_dKa_sl * dKa_sl_dCa_sl
            - 2 * dI_ncx_sl_ds2_sl * ds2_sl_dCa_sl
            - 2 * dI_ncx_sl_ds3_sl * ds3_sl_dCa_sl
            + dI_pca_sl_dCa_sl
        )
        / (2 * Frdy * Vsl)
    )
    states[36] = Ca_sl + np.where(
        np.abs(dCa_sl_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dCa_sl_dt_linearized)) * dCa_sl_dt / dCa_sl_dt_linearized,
        dt * dCa_sl_dt,
    )
    dCa_i_dt = (
        -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo
    )
    dJ_serca_dCa_i = Vmax_SRCaP * hillSRCaP * np.power(Q10SRCaP, Qpow) * np.power(
        Ca_i / Kmf, hillSRCaP
    ) / (
        (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP)) * Ca_i
    ) - Vmax_SRCaP * hillSRCaP * np.power(
        Q10SRCaP, Qpow
    ) * np.power(
        Ca_i / Kmf, hillSRCaP
    ) * (
        np.power(Ca_i / Kmf, hillSRCaP) - np.power(Ca_sr / Kmr, hillSRCaP)
    ) / (
        (
            (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
            * (1 + np.power(Ca_i / Kmf, hillSRCaP) + np.power(Ca_sr / Kmr, hillSRCaP))
        )
        * Ca_i
    )
    dJ_CaB_cytosol_dCa_i = (
        kon_cam * (Bmax_CaM - CaM)
        + kon_myoca * (Bmax_myosin - Myo_c - Myo_m)
        + kon_sr * (Bmax_SR - SRB)
        + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
        + kon_tncl * (Bmax_TnClow - Tn_CL)
    )
    dCa_i_dt_linearized = (
        -dJ_CaB_cytosol_dCa_i - J_ca_slmyo / Vmyo - Vsr * dJ_serca_dCa_i / Vmyo
    )
    states[37] = Ca_i + np.where(
        np.abs(dCa_i_dt_linearized) > 1e-08,
        (-1.0 + np.exp(dt * dCa_i_dt_linearized)) * dCa_i_dt / dCa_i_dt_linearized,
        dt * dCa_i_dt,
    )

    # Expressions for the Membrane potential component
    i_Stim = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
            t - stim_period * np.floor(t / stim_period) >= stim_start,
        ),
        -stim_amplitude,
        0,
    )
    I_Na_tot = I_Na_tot_junc + I_Na_tot_sl
    I_Cl_tot = I_ClCa + I_Clbk
    I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl
    I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot
    dV_m_dt = -I_tot - i_Stim
    dI_nak_sl_dfnak = (
        IbarNaK
        * Ko
        * Fsl
        / ((1 + np.power(KmNaip, 4) / np.power(Na_sl, 4)) * (KmKo + Ko))
    )
    dibarca_j_dV_m = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_j * np.exp(2 * FoRT * V_m))
        * FoRT
        / (-1 + np.exp(2 * FoRT * V_m))
        - 8
        * Frdy
        * GCaL
        * pCa
        * (FoRT * FoRT)
        * (-0.341 * Cao + 0.341 * Ca_j * np.exp(2 * FoRT * V_m))
        * V_m
        * np.exp(2 * FoRT * V_m)
        / ((-1 + np.exp(2 * FoRT * V_m)) * (-1 + np.exp(2 * FoRT * V_m)))
        + 2.728
        * Frdy
        * GCaL
        * pCa
        * (FoRT * FoRT)
        * Ca_j
        * V_m
        * np.exp(2 * FoRT * V_m)
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    dI_Na_sl_dV_m = GNa * (m * m * m) * Fsl * h * j
    ds1_junc_dV_m = Cao * nu * (Na_j * Na_j * Na_j) * FoRT * np.exp(nu * FoRT * V_m)
    ds2_junc_dV_m = (
        (Nao * Nao * Nao) * (-1 + nu) * Ca_j * FoRT * np.exp((-1 + nu) * FoRT * V_m)
    )
    dI_ncx_junc_dV_m = Fjunc * IbarNCX * np.power(Q10NCX, Qpow) * (
        -ds2_junc_dV_m + ds1_junc_dV_m
    ) * Ka_junc / (
        (1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_junc
    ) - Fjunc * IbarNCX * ksat * np.power(
        Q10NCX, Qpow
    ) * (
        -1 + nu
    ) * (
        -s2_junc + s1_junc
    ) * FoRT * Ka_junc * np.exp(
        (-1 + nu) * FoRT * V_m
    ) / (
        (
            (1 + ksat * np.exp((-1 + nu) * FoRT * V_m))
            * (1 + ksat * np.exp((-1 + nu) * FoRT * V_m))
        )
        * s3_junc
    )
    dibark_dV_m = (
        Frdy
        * GCaL
        * pK
        * (-0.75 * Ko + 0.75 * K_i * np.exp(FoRT * V_m))
        * FoRT
        / (-1 + np.exp(FoRT * V_m))
        - Frdy
        * GCaL
        * pK
        * (FoRT * FoRT)
        * (-0.75 * Ko + 0.75 * K_i * np.exp(FoRT * V_m))
        * V_m
        * np.exp(FoRT * V_m)
        / ((-1 + np.exp(FoRT * V_m)) * (-1 + np.exp(FoRT * V_m)))
        + 0.75
        * Frdy
        * GCaL
        * pK
        * (FoRT * FoRT)
        * K_i
        * V_m
        * np.exp(FoRT * V_m)
        / (-1 + np.exp(FoRT * V_m))
    )
    dI_kp_sl_dkp_kp = GKp * (-ek + V_m) * Fsl
    dI_CaK_dibark = (
        0.45
        * np.power(Q10CaL, Qpow)
        * (Fjunc_CaL * (1 + fcaCaj - f_Ca_Bj) + (1 + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
        * d
        * f
    )
    dI_kp_junc_dkp_kp = Fjunc * GKp * (-ek + V_m)
    dI_ClCa_sl_dV_m = GClCa * Fsl / (1 + KdClCa / Ca_sl)
    dkp_kp_dV_m = (
        298.741733340907
        * np.exp(-0.16722408026755853 * V_m)
        / (
            (1 + 1786.4755653786237 * np.exp(-0.16722408026755853 * V_m))
            * (1 + 1786.4755653786237 * np.exp(-0.16722408026755853 * V_m))
        )
    )
    drkr_dV_m = -np.exp(37 / 12 + V_m / 24) / (
        24 * ((1 + np.exp(37 / 12 + V_m / 24)) * (1 + np.exp(37 / 12 + V_m / 24)))
    )
    dI_ks_sl_dV_m = (x_ks * x_ks) * Fsl * gks_sl
    dI_tos_dV_m = GtoSlow * x_to_s * y_to_s
    dfnak_dV_m = (
        0.012450000000000001 * FoRT * np.exp(-0.1 * FoRT * V_m)
        + 0.0365 * FoRT * np.exp(-FoRT * V_m) * sigma
    ) / (
        (1 + 0.1245 * np.exp(-0.1 * FoRT * V_m) + 0.0365 * np.exp(-FoRT * V_m) * sigma)
        * (
            1
            + 0.1245 * np.exp(-0.1 * FoRT * V_m)
            + 0.0365 * np.exp(-FoRT * V_m) * sigma
        )
    )
    dbki_dV_m = (
        7.122279797276982e-18 * np.exp(0.06175 * V_m - 0.06175 * ek)
        + 0.06125396020258667 * np.exp(0.08032 * V_m - 0.08032 * ek)
    ) / (
        1 + 0.08677229415769332 * np.exp(0.5143 * ek - 0.5143 * V_m)
    ) + 0.04462699088530167 * (
        0.7626240065063081 * np.exp(0.08032 * V_m - 0.08032 * ek)
        + 1.1534056351865558e-16 * np.exp(0.06175 * V_m - 0.06175 * ek)
    ) * np.exp(
        0.5143 * ek - 0.5143 * V_m
    ) / (
        (1 + 0.08677229415769332 * np.exp(0.5143 * ek - 0.5143 * V_m))
        * (1 + 0.08677229415769332 * np.exp(0.5143 * ek - 0.5143 * V_m))
    )
    dibarca_sl_dV_m = (
        4
        * Frdy
        * GCaL
        * pCa
        * (-0.341 * Cao + 0.341 * Ca_sl * np.exp(2 * FoRT * V_m))
        * FoRT
        / (-1 + np.exp(2 * FoRT * V_m))
        - 8
        * Frdy
        * GCaL
        * pCa
        * (FoRT * FoRT)
        * (-0.341 * Cao + 0.341 * Ca_sl * np.exp(2 * FoRT * V_m))
        * V_m
        * np.exp(2 * FoRT * V_m)
        / ((-1 + np.exp(2 * FoRT * V_m)) * (-1 + np.exp(2 * FoRT * V_m)))
        + 2.728
        * Frdy
        * GCaL
        * pCa
        * (FoRT * FoRT)
        * Ca_sl
        * V_m
        * np.exp(2 * FoRT * V_m)
        / (-1 + np.exp(2 * FoRT * V_m))
    )
    dI_ks_junc_dV_m = Fjunc * (x_ks * x_ks) * gks_junc
    dI_tof_dV_m = GtoFast * x_to_f * y_to_f
    dI_ClCa_junc_dV_m = Fjunc * GClCa / (1 + KdClCa / Ca_j)
    dI_kr_drkr = (-ek + V_m) * gkr * x_kr
    daki_dV_m = (
        -1.7891395565206892e-07
        * np.exp(0.2385 * V_m - 0.2385 * ek)
        / (
            (1 + 7.35454251046446e-07 * np.exp(0.2385 * V_m - 0.2385 * ek))
            * (1 + 7.35454251046446e-07 * np.exp(0.2385 * V_m - 0.2385 * ek))
        )
    )
    dkiss_dbki = -aki / ((aki + bki) * (aki + bki))
    dkiss_daki = 1.0 / (aki + bki) - aki / ((aki + bki) * (aki + bki))
    dI_K1_dV_m = 0.4303314829119352 * GK1 * np.sqrt(
        Ko
    ) * kiss + 0.4303314829119352 * GK1 * np.sqrt(Ko) * (-ek + V_m) * (
        daki_dV_m * dkiss_daki + dbki_dV_m * dkiss_dbki
    )
    dI_kp_junc_dV_m = Fjunc * GKp * kp_kp + Fjunc * GKp * (-ek + V_m) * dkp_kp_dV_m
    dibarna_j_dV_m = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_j * np.exp(FoRT * V_m))
        * FoRT
        / (-1 + np.exp(FoRT * V_m))
        - Frdy
        * GCaL
        * pNa
        * (FoRT * FoRT)
        * (-0.75 * Nao + 0.75 * Na_j * np.exp(FoRT * V_m))
        * V_m
        * np.exp(FoRT * V_m)
        / ((-1 + np.exp(FoRT * V_m)) * (-1 + np.exp(FoRT * V_m)))
        + 0.75
        * Frdy
        * GCaL
        * pNa
        * (FoRT * FoRT)
        * Na_j
        * V_m
        * np.exp(FoRT * V_m)
        / (-1 + np.exp(FoRT * V_m))
    )
    dI_nak_junc_dfnak = (
        Fjunc
        * IbarNaK
        * Ko
        / ((1 + np.power(KmNaip, 4) / np.power(Na_j, 4)) * (KmKo + Ko))
    )
    dI_K1_dkiss = 0.4303314829119352 * GK1 * np.sqrt(Ko) * (-ek + V_m)
    ds1_sl_dV_m = Cao * nu * (Na_sl * Na_sl * Na_sl) * FoRT * np.exp(nu * FoRT * V_m)
    ds2_sl_dV_m = (
        (Nao * Nao * Nao) * (-1 + nu) * Ca_sl * FoRT * np.exp((-1 + nu) * FoRT * V_m)
    )
    dI_ncx_sl_dV_m = IbarNCX * np.power(Q10NCX, Qpow) * (
        -ds2_sl_dV_m + ds1_sl_dV_m
    ) * Fsl * Ka_sl / (
        (1 + ksat * np.exp((-1 + nu) * FoRT * V_m)) * s3_sl
    ) - IbarNCX * ksat * np.power(
        Q10NCX, Qpow
    ) * (
        -1 + nu
    ) * (
        -s2_sl + s1_sl
    ) * FoRT * Fsl * Ka_sl * np.exp(
        (-1 + nu) * FoRT * V_m
    ) / (
        (
            (1 + ksat * np.exp((-1 + nu) * FoRT * V_m))
            * (1 + ksat * np.exp((-1 + nu) * FoRT * V_m))
        )
        * s3_sl
    )
    dI_kr_dV_m = gkr * rkr * x_kr + (-ek + V_m) * drkr_dV_m * gkr * x_kr
    dI_Na_junc_dV_m = Fjunc * GNa * (m * m * m) * h * j
    dI_kp_sl_dV_m = GKp * Fsl * kp_kp + GKp * (-ek + V_m) * Fsl * dkp_kp_dV_m
    dibarna_sl_dV_m = (
        Frdy
        * GCaL
        * pNa
        * (-0.75 * Nao + 0.75 * Na_sl * np.exp(FoRT * V_m))
        * FoRT
        / (-1 + np.exp(FoRT * V_m))
        - Frdy
        * GCaL
        * pNa
        * (FoRT * FoRT)
        * (-0.75 * Nao + 0.75 * Na_sl * np.exp(FoRT * V_m))
        * V_m
        * np.exp(FoRT * V_m)
        / ((-1 + np.exp(FoRT * V_m)) * (-1 + np.exp(FoRT * V_m)))
        + 0.75
        * Frdy
        * GCaL
        * pNa
        * (FoRT * FoRT)
        * Na_sl
        * V_m
        * np.exp(FoRT * V_m)
        / (-1 + np.exp(FoRT * V_m))
    )
    dV_m_dt_linearized = (
        -GClB
        - dI_ClCa_junc_dV_m
        - dI_ClCa_sl_dV_m
        - dI_K1_dV_m
        - dI_Na_junc_dV_m
        - dI_Na_sl_dV_m
        - dI_kp_junc_dV_m
        - dI_kp_sl_dV_m
        - dI_kr_dV_m
        - dI_ks_junc_dV_m
        - dI_ks_sl_dV_m
        - dI_ncx_junc_dV_m
        - dI_ncx_sl_dV_m
        - dI_tof_dV_m
        - dI_tos_dV_m
        - Fjunc * GCaB
        - Fjunc * GNaB
        - GCaB * Fsl
        - GNaB * Fsl
        - (daki_dV_m * dkiss_daki + dbki_dV_m * dkiss_dbki) * dI_K1_dkiss
        - dI_CaK_dibark * dibark_dV_m
        - dI_CaNa_junc_dibarna_j * dibarna_j_dV_m
        - dI_CaNa_sl_dibarna_sl * dibarna_sl_dV_m
        - dI_Ca_junc_dibarca_j * dibarca_j_dV_m
        - dI_Ca_sl_dibarca_sl * dibarca_sl_dV_m
        - dI_kp_junc_dkp_kp * dkp_kp_dV_m
        - dI_kp_sl_dkp_kp * dkp_kp_dV_m
        - dI_kr_drkr * drkr_dV_m
        - dI_nak_junc_dfnak * dfnak_dV_m
        - dI_nak_sl_dfnak * dfnak_dV_m
        - dI_ncx_junc_ds1_junc * ds1_junc_dV_m
        - dI_ncx_junc_ds2_junc * ds2_junc_dV_m
        - dI_ncx_sl_ds1_sl * ds1_sl_dV_m
        - dI_ncx_sl_ds2_sl * ds2_sl_dV_m
    )
    states[38] = (
        np.where(
            np.abs(dV_m_dt_linearized) > 1e-08,
            (-1.0 + np.exp(dt * dV_m_dt_linearized)) * dV_m_dt / dV_m_dt_linearized,
            dt * dV_m_dt,
        )
        + V_m
    )

    # Return results
    return states
