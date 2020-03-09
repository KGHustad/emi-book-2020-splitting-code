function [states, varargout] = grandi_init_states()
  % % Default state values for ODE model: grandi_pasqualini_bers_2010
  % % ---------------------------------------------------------------
  % %
  % % states = grandi_init_states();
  % % [states, states_names] = grandi_init_states();

  % --- Default initial state values --- 
  states = zeros(39, 1);

  % --- I_Na ---
  states(1) = 0.00379308741444; % m;
  states(2) = 0.626221949492; % h;
  states(3) = 0.62455357249; % j;

  % --- I_Kr ---
  states(4) = 0.0210022533039; % x_kr;

  % --- I_Ks ---
  states(5) = 0.00428016666259; % x_ks;

  % --- I_to ---
  states(6) = 0.000440445885643; % x_to_s;
  states(7) = 0.785115828275; % y_to_s;
  states(8) = 0.000440438103759; % x_to_f;
  states(9) = 0.999995844039; % y_to_f;

  % --- I_Ca ---
  states(10) = 2.92407183949e-06; % d;
  states(11) = 0.995135796704; % f;
  states(12) = 0.0246760872106; % f_Ca_Bj;
  states(13) = 0.0152723084239; % f_Ca_Bsl;

  % --- SR Fluxes ---
  states(14) = 0.890806040818; % Ry_Rr;
  states(15) = 7.40481128854e-07; % Ry_Ro;
  states(16) = 9.07666168961e-08; % Ry_Ri;

  % --- Na Buffers ---
  states(17) = 3.45437733033; % Na_Bj;
  states(18) = 0.753740951478; % Na_Bsl;

  % --- Cytosolic Ca Buffers ---
  states(19) = 0.00893455096919; % Tn_CL;
  states(20) = 0.117412025937; % Tn_CHc;
  states(21) = 0.0106160166693; % Tn_CHm;
  states(22) = 0.000295573424135; % CaM;
  states(23) = 0.00192322252438; % Myo_c;
  states(24) = 0.137560495023; % Myo_m;
  states(25) = 0.00217360235649; % SRB;

  % --- Junctional and SL Ca Buffers ---
  states(26) = 0.0074052452168; % SLL_j;
  states(27) = 0.00990339304377; % SLL_sl;
  states(28) = 0.0735890020284; % SLH_j;
  states(29) = 0.114583623437; % SLH_sl;

  % --- SR Ca Concentrations ---
  states(30) = 1.19723145924; % Csqn_b;
  states(31) = 0.554760499828; % Ca_sr;

  % --- Na Concentrations ---
  states(32) = 8.40537012593; % Na_j;
  states(33) = 8.40491910001; % Na_sl;
  states(34) = 8.40513364345; % Na_i;

  % --- K Concentration ---
  states(35) = 120; % K_i;

  % --- Ca Concentrations ---
  states(36) = 0.000175882395147; % Ca_j;
  states(37) = 0.000106779509977; % Ca_sl;
  states(38) = 8.72509677797e-05; % Ca_i;

  % --- Membrane potential ---
  states(39) = -81.4552030513; % V_m;

  if nargout == 2

    % --- State names --- 
    state_names = cell(39, 1);

    % --- I_Na ---
    state_names{1} = 'm';
    state_names{2} = 'h';
    state_names{3} = 'j';

    % --- I_Kr ---
    state_names{4} = 'x_kr';

    % --- I_Ks ---
    state_names{5} = 'x_ks';

    % --- I_to ---
    state_names{6} = 'x_to_s';
    state_names{7} = 'y_to_s';
    state_names{8} = 'x_to_f';
    state_names{9} = 'y_to_f';

    % --- I_Ca ---
    state_names{10} = 'd';
    state_names{11} = 'f';
    state_names{12} = 'f_Ca_Bj';
    state_names{13} = 'f_Ca_Bsl';

    % --- SR Fluxes ---
    state_names{14} = 'Ry_Rr';
    state_names{15} = 'Ry_Ro';
    state_names{16} = 'Ry_Ri';

    % --- Na Buffers ---
    state_names{17} = 'Na_Bj';
    state_names{18} = 'Na_Bsl';

    % --- Cytosolic Ca Buffers ---
    state_names{19} = 'Tn_CL';
    state_names{20} = 'Tn_CHc';
    state_names{21} = 'Tn_CHm';
    state_names{22} = 'CaM';
    state_names{23} = 'Myo_c';
    state_names{24} = 'Myo_m';
    state_names{25} = 'SRB';

    % --- Junctional and SL Ca Buffers ---
    state_names{26} = 'SLL_j';
    state_names{27} = 'SLL_sl';
    state_names{28} = 'SLH_j';
    state_names{29} = 'SLH_sl';

    % --- SR Ca Concentrations ---
    state_names{30} = 'Csqn_b';
    state_names{31} = 'Ca_sr';

    % --- Na Concentrations ---
    state_names{32} = 'Na_j';
    state_names{33} = 'Na_sl';
    state_names{34} = 'Na_i';

    % --- K Concentration ---
    state_names{35} = 'K_i';

    % --- Ca Concentrations ---
    state_names{36} = 'Ca_j';
    state_names{37} = 'Ca_sl';
    state_names{38} = 'Ca_i';

    % --- Membrane potential ---
    state_names{39} = 'V_m';
    varargout(1) = {state_names};
  end
end