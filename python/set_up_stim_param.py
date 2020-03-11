import numpy as np

def set_up_stim_param(ionic_model, G, mesh):
    "Set up parameters for a EMI model with stimulation of the correct grid points"

    overridden_parameters = {}
    if hasattr(G, 'stim_time'):
        overridden_parameters['stim_duration'] = G.stim_time
    if hasattr(G, 'stim_start_time'):
        overridden_parameters['stim_start'] = G.stim_start_time

    parameters_single = ionic_model.init_parameter_values_single(*overridden_parameters)

    # adjust stim_amplitude
    stim_amplitude_idx = ionic_model.parameter_index('stim_amplitude')
    stim_amplitude = parameters_single[stim_amplitude_idx]
    P = ionic_model.init_parameter_values_2d_from_array(parameters_single, G.nv)
    p_tmp = np.zeros(G.N, dtype=P.dtype)
    p_tmp[mesh.to_stim] = stim_amplitude*G.stim_amp
    p_tmp = p_tmp[mesh.v]
    P[stim_amplitude_idx, :] = p_tmp

    return P
