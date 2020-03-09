function [parameters, varargout] = grandi_init_parameters()
  % % Default parameter values for ODE model: grandi_pasqualini_bers_2010
  % % -------------------------------------------------------------------
  % %
  % % parameters = grandi_init_parameters();
  % % [parameters, parameters_names] = grandi_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(107, 1);

  % --- Geometry ---
  parameters(1) = 0.11; % Fjunc;
  parameters(2) = 0.9; % Fjunc_CaL;
  parameters(3) = 100; % cellLength;
  parameters(4) = 10.25; % cellRadius;
  parameters(5) = 1; %use_I_ion% distJuncSL;
  parameters(6) = 0; %GRest distSLcyto;
  parameters(7) = 0.16; % junctionLength;
  parameters(8) = 0.015; % junctionRadius;

  % --- I_Na ---
  parameters(9) = 23; % GNa;

  % --- I_NaBK ---
  parameters(10) = 0.000597; % GNaB;

  % --- I_NaK ---
  parameters(11) = 1.8; % IbarNaK;
  parameters(12) = 1.5; % KmKo;
  parameters(13) = 11; % KmNaip;
  parameters(14) = 0.35; % g_K1 %%%1.39; % Q10KmNai;
  parameters(15) = 1.63; % Q10NaK;

  % --- I_Kp ---
  parameters(16) = 0.002; % gkp;

  % --- I_Ks ---
  parameters(17) = 0.01833; % pNaK;

  % --- I_to ---
  parameters(18) = 1; % epi;

  % --- I_ClCa ---
  parameters(19) = 0.009; % GClB;
  parameters(20) = 0.0548125; % GClCa;
  parameters(21) = 0.1; % KdClCa;

  % --- I_Ca ---
  parameters(22) = 1.8; % Q10CaL;
  parameters(23) = 0.00027; % pCa;
  parameters(24) = 1.35e-07; % pK;
  parameters(25) = 7.5e-09; % pNa;

  % --- I_NCX ---
  parameters(26) = 4.5; % IbarNCX;
  parameters(27) = 0.00015; % Kdact;
  parameters(28) = 0.00359; % KmCai;
  parameters(29) = 1.3; % KmCao;
  parameters(30) = 12.29; % KmNai;
  parameters(31) = 87.5; % KmNao;
  parameters(32) = 1.57; % Q10NCX;
  parameters(33) = 0.32; % ksat;
  parameters(34) = 0.27; % nu;

  % --- I_PCa ---
  parameters(35) = 0.0673; % IbarSLCaP;
  parameters(36) = 0.0005; % KmPCa;
  parameters(37) = 2.35; % Q10SLCaP;

  % --- I_CaBK ---
  parameters(38) = 0.0005513; % GCaB;

  % --- SR Fluxes ---
  parameters(39) = 0.000246; % Kmf;
  parameters(40) = 1.7; % Kmr;
  parameters(41) = 15; % MaxSR;
  parameters(42) = 1; % MinSR;
  parameters(43) = 2.6; % Q10SRCaP;
  parameters(44) = 0.0053114; % Vmax_SRCaP;
  parameters(45) = 0.45; % ec50SR;
  parameters(46) = 1.787; % hillSRCaP;
  parameters(47) = 0.5; % kiCa;
  parameters(48) = 0.005; % kim;
  parameters(49) = 10; % koCa;
  parameters(50) = 0.06; % kom;
  parameters(51) = 25; % ks;

  % --- Na Buffers ---
  parameters(52) = 7.561; % Bmax_Naj;
  parameters(53) = 1.65; % Bmax_Nasl;
  parameters(54) = 0.001; % koff_na;
  parameters(55) = 0.0001; % kon_na;

  % --- Cytosolic Ca Buffers ---
  parameters(56) = 0.024; % Bmax_CaM;
  parameters(57) = 0.0171; % Bmax_SR;
  parameters(58) = 0.14; % Bmax_TnChigh;
  parameters(59) = 0.07; % Bmax_TnClow;
  parameters(60) = 0.14; % Bmax_myosin;
  parameters(61) = 0.238; % koff_cam;
  parameters(62) = 0.00046; % koff_myoca;
  parameters(63) = 5.7e-05; % koff_myomg;
  parameters(64) = 0.06; % koff_sr;
  parameters(65) = 3.2e-05; % koff_tnchca;
  parameters(66) = 0.00333; % koff_tnchmg;
  parameters(67) = 0.0196; % koff_tncl;
  parameters(68) = 34; % kon_cam;
  parameters(69) = 13.8; % kon_myoca;
  parameters(70) = 0.0157; % kon_myomg;
  parameters(71) = 100; % kon_sr;
  parameters(72) = 2.37; % kon_tnchca;
  parameters(73) = 0.003; % kon_tnchmg;
  parameters(74) = 32.7; % kon_tncl;

  % --- Junctional and SL Ca Buffers ---
  parameters(75) = 0.000165; % Bmax_SLhighj0;
  parameters(76) = 0.0134; % Bmax_SLhighsl0;
  parameters(77) = 0.00046; % Bmax_SLlowj0;
  parameters(78) = 0.0374; % Bmax_SLlowsl0;
  parameters(79) = 0.03; % koff_slh;
  parameters(80) = 1.3; % koff_sll;
  parameters(81) = 100; % kon_slh;
  parameters(82) = 100; % kon_sll;

  % --- SR Ca Concentrations ---
  parameters(83) = 0.14; % Bmax_Csqn0;
  parameters(84) = 1.64e-06; % DcaJuncSL;
  parameters(85) = 1.22e-06; % DcaSLcyto;
  parameters(86) = 8.2413e-13; % J_ca_juncsl;
  parameters(87) = 3.7243e-12; % J_ca_slmyo;
  parameters(88) = 65; % koff_csqn;
  parameters(89) = 100; % kon_csqn;

  % --- Na Concentrations ---
  parameters(90) = 1.09e-05; % DnaJuncSL;
  parameters(91) = 1.79e-05; % DnaSLcyto;
  parameters(92) = 1.8313e-14; % J_na_juncsl;
  parameters(93) = 1.6386e-12; % J_na_slmyo;
  parameters(94) = 140; % Nao;

  % --- K Concentration ---
  parameters(95) = 5.4; % Ko;

  % --- Ca Concentrations ---
  parameters(96) = 1.8; % Cao;

  % --- Cl Concentrations ---
  parameters(97) = 15; % Cli;
  parameters(98) = 150; % Clo;

  % --- Mg Concentrations ---
  parameters(99) = 1; % Mgi;

  % --- Membrane potential ---
  parameters(100) = 1.381e-10; % Cmem;
  parameters(101) = 96485; % Frdy;
  parameters(102) = 8314; % R;
  parameters(103) = 310; % Temp;
  parameters(104) = 40.0; % stim_amplitude;
  parameters(105) = 1.0; % stim_duration;
  parameters(106) = 1000.0; % stim_period;
  parameters(107) = 0.0; % stim_start;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(107, 1);

    % --- Geometry ---
    parameter_names{1} = 'Fjunc';
    parameter_names{2} = 'Fjunc_CaL';
    parameter_names{3} = 'cellLength';
    parameter_names{4} = 'cellRadius';
    parameter_names{5} = 'distJuncSL';
    parameter_names{6} = 'distSLcyto';
    parameter_names{7} = 'junctionLength';
    parameter_names{8} = 'junctionRadius';

    % --- I_Na ---
    parameter_names{9} = 'GNa';

    % --- I_NaBK ---
    parameter_names{10} = 'GNaB';

    % --- I_NaK ---
    parameter_names{11} = 'IbarNaK';
    parameter_names{12} = 'KmKo';
    parameter_names{13} = 'KmNaip';
    parameter_names{14} = 'g_K1';
    parameter_names{15} = 'Q10NaK';

    % --- I_Kp ---
    parameter_names{16} = 'gkp';

    % --- I_Ks ---
    parameter_names{17} = 'pNaK';

    % --- I_to ---
    parameter_names{18} = 'epi';

    % --- I_ClCa ---
    parameter_names{19} = 'GClB';
    parameter_names{20} = 'GClCa';
    parameter_names{21} = 'KdClCa';

    % --- I_Ca ---
    parameter_names{22} = 'Q10CaL';
    parameter_names{23} = 'pCa';
    parameter_names{24} = 'pK';
    parameter_names{25} = 'pNa';

    % --- I_NCX ---
    parameter_names{26} = 'IbarNCX';
    parameter_names{27} = 'Kdact';
    parameter_names{28} = 'KmCai';
    parameter_names{29} = 'KmCao';
    parameter_names{30} = 'KmNai';
    parameter_names{31} = 'KmNao';
    parameter_names{32} = 'Q10NCX';
    parameter_names{33} = 'ksat';
    parameter_names{34} = 'nu';

    % --- I_PCa ---
    parameter_names{35} = 'IbarSLCaP';
    parameter_names{36} = 'KmPCa';
    parameter_names{37} = 'Q10SLCaP';

    % --- I_CaBK ---
    parameter_names{38} = 'GCaB';

    % --- SR Fluxes ---
    parameter_names{39} = 'Kmf';
    parameter_names{40} = 'Kmr';
    parameter_names{41} = 'MaxSR';
    parameter_names{42} = 'MinSR';
    parameter_names{43} = 'Q10SRCaP';
    parameter_names{44} = 'Vmax_SRCaP';
    parameter_names{45} = 'ec50SR';
    parameter_names{46} = 'hillSRCaP';
    parameter_names{47} = 'kiCa';
    parameter_names{48} = 'kim';
    parameter_names{49} = 'koCa';
    parameter_names{50} = 'kom';
    parameter_names{51} = 'ks';

    % --- Na Buffers ---
    parameter_names{52} = 'Bmax_Naj';
    parameter_names{53} = 'Bmax_Nasl';
    parameter_names{54} = 'koff_na';
    parameter_names{55} = 'kon_na';

    % --- Cytosolic Ca Buffers ---
    parameter_names{56} = 'Bmax_CaM';
    parameter_names{57} = 'Bmax_SR';
    parameter_names{58} = 'Bmax_TnChigh';
    parameter_names{59} = 'Bmax_TnClow';
    parameter_names{60} = 'Bmax_myosin';
    parameter_names{61} = 'koff_cam';
    parameter_names{62} = 'koff_myoca';
    parameter_names{63} = 'koff_myomg';
    parameter_names{64} = 'koff_sr';
    parameter_names{65} = 'koff_tnchca';
    parameter_names{66} = 'koff_tnchmg';
    parameter_names{67} = 'koff_tncl';
    parameter_names{68} = 'kon_cam';
    parameter_names{69} = 'kon_myoca';
    parameter_names{70} = 'kon_myomg';
    parameter_names{71} = 'kon_sr';
    parameter_names{72} = 'kon_tnchca';
    parameter_names{73} = 'kon_tnchmg';
    parameter_names{74} = 'kon_tncl';

    % --- Junctional and SL Ca Buffers ---
    parameter_names{75} = 'Bmax_SLhighj0';
    parameter_names{76} = 'Bmax_SLhighsl0';
    parameter_names{77} = 'Bmax_SLlowj0';
    parameter_names{78} = 'Bmax_SLlowsl0';
    parameter_names{79} = 'koff_slh';
    parameter_names{80} = 'koff_sll';
    parameter_names{81} = 'kon_slh';
    parameter_names{82} = 'kon_sll';

    % --- SR Ca Concentrations ---
    parameter_names{83} = 'Bmax_Csqn0';
    parameter_names{84} = 'DcaJuncSL';
    parameter_names{85} = 'DcaSLcyto';
    parameter_names{86} = 'J_ca_juncsl';
    parameter_names{87} = 'J_ca_slmyo';
    parameter_names{88} = 'koff_csqn';
    parameter_names{89} = 'kon_csqn';

    % --- Na Concentrations ---
    parameter_names{90} = 'DnaJuncSL';
    parameter_names{91} = 'DnaSLcyto';
    parameter_names{92} = 'J_na_juncsl';
    parameter_names{93} = 'J_na_slmyo';
    parameter_names{94} = 'Nao';

    % --- K Concentration ---
    parameter_names{95} = 'Ko';

    % --- Ca Concentrations ---
    parameter_names{96} = 'Cao';

    % --- Cl Concentrations ---
    parameter_names{97} = 'Cli';
    parameter_names{98} = 'Clo';

    % --- Mg Concentrations ---
    parameter_names{99} = 'Mgi';

    % --- Membrane potential ---
    parameter_names{100} = 'Cmem';
    parameter_names{101} = 'Frdy';
    parameter_names{102} = 'R';
    parameter_names{103} = 'Temp';
    parameter_names{104} = 'stim_amplitude';
    parameter_names{105} = 'stim_duration';
    parameter_names{106} = 'stim_period';
    parameter_names{107} = 'stim_start';
    varargout(1) = {parameter_names};
  end
end