function [values] = grandi_rhs(t, states, parameters)
  % Compute the right hand side of the grandi_pasqualini_bers_2010 ODE

  % Assign states
  m=states(1,:); h=states(2,:); j=states(3,:); x_kr=states(4,:); x_ks=states(5,:);...
    x_to_s=states(6,:); y_to_s=states(7,:); x_to_f=states(8,:); y_to_f=states(9,:);...
    d=states(10,:); f=states(11,:); f_Ca_Bj=states(12,:); f_Ca_Bsl=states(13,:);...
    Ry_Rr=states(14,:); Ry_Ro=states(15,:); Ry_Ri=states(16,:); Na_Bj=states(17,:);...
    Na_Bsl=states(18,:); Tn_CL=states(19,:); Tn_CHc=states(20,:);...
    Tn_CHm=states(21,:); CaM=states(22,:); Myo_c=states(23,:); Myo_m=states(24,:);...
    SRB=states(25,:); SLL_j=states(26,:); SLL_sl=states(27,:); SLH_j=states(28,:);...
    SLH_sl=states(29,:); Csqn_b=states(30,:); Ca_sr=states(31,:); Na_j=states(32,:);...
    Na_sl=states(33,:); Na_i=states(34,:); K_i=states(35,:); Ca_j=states(36,:);...
    Ca_sl=states(37,:); Ca_i=states(38,:); V_m=states(39,:);

  % Assign parameters
  Fjunc=parameters(1,:); Fjunc_CaL=parameters(2,:); cellLength=parameters(3,:);...
    cellRadius=parameters(4,:); GNa=parameters(9,:); GNaB=parameters(10,:);...
    IbarNaK=parameters(11,:); KmKo=parameters(12,:); KmNaip=parameters(13,:);...
    g_K1=parameters(14,:); gkp=parameters(16,:); pNaK=parameters(17,:); epi=parameters(18,:);...
    GClB=parameters(19,:); GClCa=parameters(20,:); KdClCa=parameters(21,:);...
    Q10CaL=parameters(22,:); pCa=parameters(23,:); pK=parameters(24,:);...
    pNa=parameters(25,:); IbarNCX=parameters(26,:); Kdact=parameters(27,:);...
    KmCai=parameters(28,:); KmCao=parameters(29,:); KmNai=parameters(30,:);...
    KmNao=parameters(31,:); Q10NCX=parameters(32,:); ksat=parameters(33,:);...
    nu=parameters(34,:); IbarSLCaP=parameters(35,:); KmPCa=parameters(36,:);...
    Q10SLCaP=parameters(37,:); GCaB=parameters(38,:); Kmf=parameters(39,:);...
    Kmr=parameters(40,:); MaxSR=parameters(41,:); MinSR=parameters(42,:);...
    Q10SRCaP=parameters(43,:); Vmax_SRCaP=parameters(44,:);...
    ec50SR=parameters(45,:); hillSRCaP=parameters(46,:); kiCa=parameters(47,:);...
    kim=parameters(48,:); koCa=parameters(49,:); kom=parameters(50,:);...
    ks=parameters(51,:); Bmax_Naj=parameters(52,:); Bmax_Nasl=parameters(53,:);...
    koff_na=parameters(54,:); kon_na=parameters(55,:); Bmax_CaM=parameters(56,:);...
    Bmax_SR=parameters(57,:); Bmax_TnChigh=parameters(58,:);...
    Bmax_TnClow=parameters(59,:); Bmax_myosin=parameters(60,:);...
    koff_cam=parameters(61,:); koff_myoca=parameters(62,:);...
    koff_myomg=parameters(63,:); koff_sr=parameters(64,:);...
    koff_tnchca=parameters(65,:); koff_tnchmg=parameters(66,:);...
    koff_tncl=parameters(67,:); kon_cam=parameters(68,:);...
    kon_myoca=parameters(69,:); kon_myomg=parameters(70,:);...
    kon_sr=parameters(71,:); kon_tnchca=parameters(72,:);...
    kon_tnchmg=parameters(73,:); kon_tncl=parameters(74,:);...
    Bmax_SLhighj0=parameters(75,:); Bmax_SLhighsl0=parameters(76,:);...
    Bmax_SLlowj0=parameters(77,:); Bmax_SLlowsl0=parameters(78,:);...
    koff_slh=parameters(79,:); koff_sll=parameters(80,:); kon_slh=parameters(81,:);...
    kon_sll=parameters(82,:); Bmax_Csqn0=parameters(83,:);...
    J_ca_juncsl=parameters(86,:); J_ca_slmyo=parameters(87,:);...
    koff_csqn=parameters(88,:); kon_csqn=parameters(89,:);...
    J_na_juncsl=parameters(92,:); J_na_slmyo=parameters(93,:);...
    Nao=parameters(94,:); Ko=parameters(95,:); Cao=parameters(96,:);...
    Cli=parameters(97,:); Clo=parameters(98,:); Mgi=parameters(99,:);...
    Cmem=parameters(100,:); Frdy=parameters(101,:); R=parameters(102,:);...
    Temp=parameters(103,:); stim_amplitude=parameters(104,:);...
    stim_duration=parameters(105,:); stim_period=parameters(106,:);...
    stim_start=parameters(107,:);

  % Init return args
  values = zeros(size(states));

  % Expressions for the Geometry component
  Vcell = 1e-15.*pi*cellLength.*cellRadius.^2;
  Vmyo = 0.65*Vcell;
  Vsr = 0.035*Vcell;
  Vsl = 0.02*Vcell;  % /10
  Vjunc = 0.000539*Vcell;
  Fsl = 1 - Fjunc;
  Fsl_CaL = 1 - Fjunc_CaL;

  % Expressions for the Reversal potentials component
  FoRT = Frdy./(R.*Temp);
  ena_junc = log(Nao./Na_j)./FoRT;
  ena_sl = log(Nao./Na_sl)./FoRT;
  ek = log(Ko./K_i)./FoRT;
  eca_junc = log(Cao./Ca_j)./(2.*FoRT);
  eca_sl = log(Cao./Ca_sl)./(2.*FoRT);
  ecl = log(Cli./Clo)./FoRT;
  Qpow = -31 + Temp./10;

  % Expressions for the I_Na component
  mss = (1 + 0.00184221158117.*exp(-0.110741971207.*V_m)).^(-2);
  taum = 0.1292.*exp(-(2.94658944659 + 0.0643500643501.*V_m).^2) +...
    0.06487.*exp(-(-0.0943466353678 + 0.0195618153365.*V_m).^2);
  ah = ((V_m >= -40).*(0) + ~(V_m >=...
    -40).*(4.43126792958e-07.*exp(-0.147058823529.*V_m)));
  bh = ((V_m >= -40).*(0.77./(0.13 +...
    0.0497581410839.*exp(-0.0900900900901.*V_m))) + ~(V_m >=...
    -40).*(310000.0.*exp(0.3485.*V_m) + 2.7.*exp(0.079.*V_m)));
  tauh = 1.0./(ah + bh);
  hss = (1 + 15212.5932857.*exp(0.134589502019.*V_m)).^(-2);
  aj = ((V_m >= -40).*(0) + ~(V_m >= -40).*((37.78 +...
    V_m).*(-25428.0.*exp(0.2444.*V_m) - 6.948e-06.*exp(-0.04391.*V_m))./(1 +...
    50262745826.0.*exp(0.311.*V_m))));
  bj = ((V_m >= -40).*(0.6.*exp(0.057.*V_m)./(1 + 0.0407622039784.*exp(-0.1.*V_m)))...
    + ~(V_m >= -40).*(0.02424.*exp(-0.01052.*V_m)./(1 +...
    0.0039608683399.*exp(-0.1378.*V_m))));
  tauj = 1.0./(aj + bj);
  jss = (1 + 15212.5932857.*exp(0.134589502019.*V_m)).^(-2);
  values(1,:) = (-m + mss)./taum;
  values(2,:) = (-h + hss)./tauh;
  values(3,:) = (-j + jss)./tauj;
  I_Na_junc = Fjunc.*GNa.*m.^3.*(-ena_junc + V_m).*h.*j;
  I_Na_sl = GNa.*m.^3.*(-ena_sl + V_m).*Fsl.*h.*j;

  % Expressions for the I_NaBK component
  I_nabk_junc = Fjunc.*GNaB.*(-ena_junc + V_m);
  I_nabk_sl = GNaB.*(-ena_sl + V_m).*Fsl;

  % Expressions for the I_NaK component
  sigma = -1/7 + exp(0.0148588410104.*Nao)/7;
  fnak = 1.0./(1 + 0.1245.*exp(-0.1.*FoRT.*V_m) + 0.0365.*exp(-FoRT.*V_m).*sigma);
  I_nak_junc = Fjunc.*IbarNaK.*Ko.*fnak./((1 + KmNaip.^4./Na_j.^4).*(KmKo + Ko));
  I_nak_sl = IbarNaK.*Ko.*Fsl.*fnak./((1 + KmNaip.^4./Na_sl.^4).*(KmKo + Ko));
  I_nak = I_nak_junc + I_nak_sl;

  % Expressions for the I_Kr component
  gkr = 0.0150616019019.*sqrt(Ko);
  xrss = 1.0./(1 + exp(-2 - V_m./5));
  tauxr = 230./(1 + exp(2 + V_m./20)) + 3300./((1 + exp(-22./9 - V_m./9)).*(1 +...
    exp(11./9 + V_m./9)));
  values(4,:) = (-x_kr + xrss)./tauxr;
  rkr = 1.0./(1 + exp(37./12 + V_m./24));
  I_kr = (-ek + V_m).*gkr.*rkr.*x_kr;

  % Expressions for the I_Kp component
  kp_kp = 1.0./(1 + 1786.47556538.*exp(-0.167224080268.*V_m));
  I_kp_junc = Fjunc.*gkp.*(-ek + V_m).*kp_kp;
  I_kp_sl = gkp.*(-ek + V_m).*Fsl.*kp_kp;
  I_kp = I_kp_junc + I_kp_sl;

  % Expressions for the I_Ks component
  eks = log((Ko + Nao.*pNaK)./(pNaK.*Na_i + K_i))./FoRT;
  gks_junc = 0.0035;
  gks_sl = 0.0035;
  xsss = 1.0./(1 + 0.765928338365.*exp(-0.0701754385965.*V_m));
  tauxs = 990.1./(1 + 0.841540408868.*exp(-0.070821529745.*V_m));
  values(5,:) = (-x_ks + xsss)./tauxs;
  I_ks_junc = Fjunc.*gks_junc.*x_ks.^2.*(-eks + V_m);
  I_ks_sl = gks_sl.*x_ks.^2.*(-eks + V_m).*Fsl;
  I_ks = I_ks_junc + I_ks_sl;

  % Expressions for the I_to component
  GtoSlow = ((epi == 1).*(0.0156) + ~(epi == 1).*(0.037596));
  GtoFast = ((epi == 1).*(0.1144) + ~(epi == 1).*(0.001404));
  xtoss = 1.0./(1 + exp(19./13 - V_m./13));
  ytoss = 1.0./(1 + 49.4024491055.*exp(V_m./5));
  tauxtos = 0.5 + 9./(1 + exp(1./5 + V_m./15));
  tauytos = 30 + 800./(1 + exp(6 + V_m./10));
  values(6,:) = (-x_to_s + xtoss)./tauxtos;
  values(7,:) = (-y_to_s + ytoss)./tauytos;
  I_tos = (-ek + V_m).*GtoSlow.*x_to_s.*y_to_s;
  tauxtof = 0.5 + 8.5.*exp(-(9./10 + V_m./50).^2);
  tauytof = 7 + 85.*exp(-(40 + V_m).^2./220);
  values(8,:) = (-x_to_f + xtoss)./tauxtof;
  values(9,:) = (-y_to_f + ytoss)./tauytof;
  I_tof = (-ek + V_m).*GtoFast.*x_to_f.*y_to_f;
  I_to = I_tof + I_tos;

  % Expressions for the I_Ki component
  aki = 1.02./(1 + 7.35454251046e-07.*exp(0.2385.*V_m - 0.2385.*ek));
  bki = (0.762624006506.*exp(0.08032.*V_m - 0.08032.*ek) +...
    1.15340563519e-16.*exp(0.06175.*V_m - 0.06175.*ek))./(1 +...
    0.0867722941577.*exp(0.5143.*ek - 0.5143.*V_m));
  kiss = aki./(aki + bki);
  I_ki = g_K1.*sqrt(Ko/5.4).*(-ek + V_m).*kiss;

  % Expressions for the I_ClCa component
  I_ClCa_junc = Fjunc.*GClCa.*(-ecl + V_m)./(1 + KdClCa./Ca_j);
  I_ClCa_sl = GClCa.*(-ecl + V_m).*Fsl./(1 + KdClCa./Ca_sl);
  I_ClCa = I_ClCa_junc + I_ClCa_sl;
  I_Clbk = GClB.*(-ecl + V_m);

  % Expressions for the I_Ca component
  fss = 1.0./(1 + exp(35./9 + V_m./9)) + 0.6./(1 + exp(5./2 - V_m./20));
  dss = 1.0./(1 + exp(-5./6 - V_m./6));
  taud = (1 - exp(-5./6 - V_m./6)).*dss./(0.175 + 0.035.*V_m);
  tauf = 1.0./(0.02 + 0.0197.*exp(-(0.48865 + 0.0337.*V_m).^2));
  values(10,:) = (-d + dss)./taud;
  values(11,:) = (-f + fss)./tauf;
  values(12,:) = -0.0119.*f_Ca_Bj + 1.7.*(1 - f_Ca_Bj).*Ca_j;
  values(13,:) = -0.0119.*f_Ca_Bsl + 1.7.*(1 - f_Ca_Bsl).*Ca_sl;
  fcaCaMSL = 0;
  fcaCaj = 0;
  ibarca_j = 4.*Frdy.*pCa.*(-0.341.*Cao +...
    0.341.*Ca_j.*exp(2.*FoRT.*V_m)).*FoRT.*V_m./(-1 + exp(2.*FoRT.*V_m));
  ibarca_sl = 4.*Frdy.*pCa.*(-0.341.*Cao +...
    0.341.*Ca_sl.*exp(2.*FoRT.*V_m)).*FoRT.*V_m./(-1 + exp(2.*FoRT.*V_m));
  ibark = Frdy.*pK.*(-0.75.*Ko + 0.75.*K_i.*exp(FoRT.*V_m)).*FoRT.*V_m./(-1 +...
    exp(FoRT.*V_m));
  ibarna_j = Frdy.*pNa.*(-0.75.*Nao + 0.75.*Na_j.*exp(FoRT.*V_m)).*FoRT.*V_m./(-1 +...
    exp(FoRT.*V_m));
  ibarna_sl = Frdy.*pNa.*(-0.75.*Nao + 0.75.*Na_sl.*exp(FoRT.*V_m)).*FoRT.*V_m./(-1 +...
    exp(FoRT.*V_m));
  I_Ca_junc = 0.45.*Fjunc_CaL.*Q10CaL.^Qpow.*(1 + fcaCaj - f_Ca_Bj).*d.*f.*ibarca_j;
  I_Ca_sl = 0.45.*Q10CaL.^Qpow.*(1 + fcaCaMSL - f_Ca_Bsl).*Fsl_CaL.*d.*f.*ibarca_sl;
  I_CaK = 0.45.*Q10CaL.^Qpow.*(Fjunc_CaL.*(1 + fcaCaj - f_Ca_Bj) + (1 + fcaCaMSL...
    - f_Ca_Bsl).*Fsl_CaL).*d.*f.*ibark;
  I_CaNa_junc = 0.45.*Fjunc_CaL.*Q10CaL.^Qpow.*(1 + fcaCaj - f_Ca_Bj).*d.*f.*ibarna_j;
  I_CaNa_sl = 0.45.*Q10CaL.^Qpow.*(1 + fcaCaMSL - f_Ca_Bsl).*Fsl_CaL.*d.*f.*ibarna_sl;

  % Expressions for the I_NCX component
  Ka_junc = 1.0./(1 + Kdact.^2./Ca_j.^2);
  Ka_sl = 1.0./(1 + Kdact.^2./Ca_sl.^2);
  s1_junc = Cao.*Na_j.^3.*exp(nu.*FoRT.*V_m);
  s1_sl = Cao.*Na_sl.^3.*exp(nu.*FoRT.*V_m);
  s2_junc = Nao.^3.*Ca_j.*exp((-1 + nu).*FoRT.*V_m);
  s3_junc = Cao.*Na_j.^3 + KmCao.*Na_j.^3 + Nao.^3.*Ca_j + KmCai.*Nao.^3.*(1 +...
    Na_j.^3./KmNai.^3) + KmNao.^3.*(1 + Ca_j./KmCai).*Ca_j;
  s2_sl = Nao.^3.*Ca_sl.*exp((-1 + nu).*FoRT.*V_m);
  s3_sl = Cao.*Na_sl.^3 + KmCao.*Na_sl.^3 + Nao.^3.*Ca_sl + KmCai.*Nao.^3.*(1 +...
    Na_sl.^3./KmNai.^3) + KmNao.^3.*(1 + Ca_sl./KmCai).*Ca_sl;
  I_ncx_junc = Fjunc.*IbarNCX.*Q10NCX.^Qpow.*(-s2_junc + s1_junc).*Ka_junc./((1 +...
    ksat.*exp((-1 + nu).*FoRT.*V_m)).*s3_junc);
  I_ncx_sl = IbarNCX.*Q10NCX.^Qpow.*(-s2_sl + s1_sl).*Fsl.*Ka_sl./((1 +...
    ksat.*exp((-1 + nu).*FoRT.*V_m)).*s3_sl);

  % Expressions for the I_PCa component
  I_pca_junc = Fjunc.*IbarSLCaP.*Q10SLCaP.^Qpow.*Ca_j.^1.6./(KmPCa.^1.6 + Ca_j.^1.6);
  I_pca_sl = IbarSLCaP.*Q10SLCaP.^Qpow.*Ca_sl.^1.6.*Fsl./(KmPCa.^1.6 + Ca_sl.^1.6);

  % Expressions for the I_CaBK component
  I_cabk_junc = Fjunc.*GCaB.*(-eca_junc + V_m);
  I_cabk_sl = GCaB.*(-eca_sl + V_m).*Fsl;

  % Expressions for the SR Fluxes component
  kCaSR = MaxSR - (MaxSR - MinSR)./(1 + (ec50SR./Ca_sr).^2.5);
  koSRCa = koCa./kCaSR;
  kiSRCa = kiCa.*kCaSR;
  RI = 1 - Ry_Ri - Ry_Ro - Ry_Rr;
  values(14,:) = kim.*RI + kom.*Ry_Ro - Ca_j.^2.*Ry_Rr.*koSRCa - Ca_j.*Ry_Rr.*kiSRCa;
  values(15,:) = kim.*Ry_Ri - kom.*Ry_Ro + Ca_j.^2.*Ry_Rr.*koSRCa - Ca_j.*Ry_Ro.*kiSRCa;
  values(16,:) = -kim.*Ry_Ri - kom.*Ry_Ri + Ca_j.^2.*RI.*koSRCa + Ca_j.*Ry_Ro.*kiSRCa;
  J_SRCarel = ks.*(-Ca_j + Ca_sr).*Ry_Ro;
  J_serca = Vmax_SRCaP.*Q10SRCaP.^Qpow.*((Ca_i./Kmf).^hillSRCaP -...
    (Ca_sr./Kmr).^hillSRCaP)./(1 + (Ca_i./Kmf).^hillSRCaP + (Ca_sr./Kmr).^hillSRCaP);
  J_SRleak = 5.348e-06.*Ca_sr - 5.348e-06.*Ca_j;

  % Expressions for the Na Buffers component
  values(17,:) = -koff_na.*Na_Bj + kon_na.*(Bmax_Naj - Na_Bj).*Na_j;
  values(18,:) = -koff_na.*Na_Bsl + kon_na.*(Bmax_Nasl - Na_Bsl).*Na_sl;

  % Expressions for the Cytosolic Ca Buffers component
  values(19,:) = -koff_tncl.*Tn_CL + kon_tncl.*(Bmax_TnClow - Tn_CL).*Ca_i;
  values(20,:) = -koff_tnchca.*Tn_CHc + kon_tnchca.*(Bmax_TnChigh - Tn_CHc -...
    Tn_CHm).*Ca_i;
  values(21,:) = -koff_tnchmg.*Tn_CHm + Mgi.*kon_tnchmg.*(Bmax_TnChigh - Tn_CHc -...
    Tn_CHm);
  values(22,:) = -koff_cam.*CaM + kon_cam.*(Bmax_CaM - CaM).*Ca_i;
  values(23,:) = -koff_myoca.*Myo_c + kon_myoca.*(Bmax_myosin - Myo_c -...
    Myo_m).*Ca_i;
  values(24,:) = -koff_myomg.*Myo_m + Mgi.*kon_myomg.*(Bmax_myosin - Myo_c - Myo_m);
  values(25,:) = -koff_sr.*SRB + kon_sr.*(Bmax_SR - SRB).*Ca_i;
  J_CaB_cytosol = -koff_cam.*CaM - koff_myoca.*Myo_c - koff_myomg.*Myo_m -...
    koff_sr.*SRB - koff_tnchca.*Tn_CHc - koff_tnchmg.*Tn_CHm - koff_tncl.*Tn_CL +...
    Mgi.*kon_myomg.*(Bmax_myosin - Myo_c - Myo_m) +...
    Mgi.*kon_tnchmg.*(Bmax_TnChigh - Tn_CHc - Tn_CHm) + kon_cam.*(Bmax_CaM -...
    CaM).*Ca_i + kon_myoca.*(Bmax_myosin - Myo_c - Myo_m).*Ca_i +...
    kon_sr.*(Bmax_SR - SRB).*Ca_i + kon_tnchca.*(Bmax_TnChigh - Tn_CHc -...
    Tn_CHm).*Ca_i + kon_tncl.*(Bmax_TnClow - Tn_CL).*Ca_i;

  % Expressions for the Junctional and SL Ca Buffers component
  Bmax_SLlowsl = Bmax_SLlowsl0.*Vmyo./Vsl;
  Bmax_SLlowj = Bmax_SLlowj0.*Vmyo./Vjunc;
  Bmax_SLhighsl = Bmax_SLhighsl0.*Vmyo./Vsl;
  Bmax_SLhighj = Bmax_SLhighj0.*Vmyo./Vjunc;
  values(26,:) = -koff_sll.*SLL_j + kon_sll.*(-SLL_j + Bmax_SLlowj).*Ca_j;
  values(27,:) = -koff_sll.*SLL_sl + kon_sll.*(-SLL_sl + Bmax_SLlowsl).*Ca_sl;
  values(28,:) = -koff_slh.*SLH_j + kon_slh.*(-SLH_j + Bmax_SLhighj).*Ca_j;
  values(29,:) = -koff_slh.*SLH_sl + kon_slh.*(-SLH_sl + Bmax_SLhighsl).*Ca_sl;
  J_CaB_junction = -koff_slh.*SLH_j - koff_sll.*SLL_j + kon_slh.*(-SLH_j +...
    Bmax_SLhighj).*Ca_j + kon_sll.*(-SLL_j + Bmax_SLlowj).*Ca_j;
  J_CaB_sl = -koff_slh.*SLH_sl - koff_sll.*SLL_sl + kon_slh.*(-SLH_sl +...
    Bmax_SLhighsl).*Ca_sl + kon_sll.*(-SLL_sl + Bmax_SLlowsl).*Ca_sl;

  % Expressions for the SR Ca Concentrations component
  Bmax_Csqn = Bmax_Csqn0.*Vmyo./Vsr;
  values(30,:) = -koff_csqn.*Csqn_b + kon_csqn.*(-Csqn_b + Bmax_Csqn).*Ca_sr;
  values(31,:) = -J_SRCarel + koff_csqn.*Csqn_b - kon_csqn.*(-Csqn_b +...
    Bmax_Csqn).*Ca_sr - J_SRleak.*Vmyo./Vsr + J_serca;

  % Expressions for the Na Concentrations component
  I_Na_tot_junc = 3.*I_nak_junc + 3.*I_ncx_junc + I_CaNa_junc + I_Na_junc +...
    I_nabk_junc;
  I_Na_tot_sl = 3.*I_nak_sl + 3.*I_ncx_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl;
  values(32,:) = -values(17,:) + J_na_juncsl.*(-Na_j + Na_sl)./Vjunc -...
    Cmem.*I_Na_tot_junc./(Frdy.*Vjunc);
  values(33,:) = -values(18,:) + J_na_juncsl.*(-Na_sl + Na_j)./Vsl +...
    J_na_slmyo.*(-Na_sl + Na_i)./Vsl - Cmem.*I_Na_tot_sl./(Frdy.*Vsl);
  values(34,:) = J_na_slmyo.*(-Na_i + Na_sl)./Vmyo;

  % Expressions for the K Concentration component
  I_K_tot = -2.*I_nak + I_CaK + I_ki + I_kp + I_kr + I_ks + I_to;
  values(35,:) = 0;

  % Expressions for the Ca Concentrations component
  I_Ca_tot_junc = -2.*I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc;
  I_Ca_tot_sl = -2.*I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl;
  values(36,:) = -J_CaB_junction + J_ca_juncsl.*(-Ca_j + Ca_sl)./Vjunc +...
    J_SRCarel.*Vsr./Vjunc + J_SRleak.*Vmyo./Vjunc -...
    Cmem.*I_Ca_tot_junc./(2.*Frdy.*Vjunc);
  values(37,:) = -J_CaB_sl + J_ca_juncsl.*(-Ca_sl + Ca_j)./Vsl +...
    J_ca_slmyo.*(-Ca_sl + Ca_i)./Vsl - Cmem.*I_Ca_tot_sl./(2.*Frdy.*Vsl);
  values(38,:) = -J_CaB_cytosol + J_ca_slmyo.*(-Ca_i + Ca_sl)./Vmyo -...
    J_serca.*Vsr./Vmyo;

  % Expressions for the Membrane potential component
  i_Stim = ((t - stim_period.*floor(t./stim_period) <= stim_duration +...
    stim_start & t - stim_period.*floor(t./stim_period) >=...
    stim_start).*(-stim_amplitude));
  I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
  I_Cl_tot = I_ClCa + I_Clbk;
  I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
  I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot;
  values(39,:) = -I_tot - i_Stim;
end