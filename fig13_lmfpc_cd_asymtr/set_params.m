function [q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, delta_phi, kvalue] = set_params
% Set Plasma Parameters for iSH
%  set_params sets the charge and mass needed by the Lorentz Force Law
%  and the dimensionless plasma parameters
%     vtic = v_ti / c = 1e-4 from PLUME
%     bi = Ion Plasma beta
%     be = Electron Plasma beta
  
% Set species parameters (Ions: q=+1, m=mp/me; Electrons: q=-1, m=1)
q = 1.;
m = 1.;
mime = 1836.;  % Set mi/me
tite = 1.; % set Ti/Te

% from PLUME: vtp = 1.E-4
vtic = 1.e-4;

% RSR
em_eps = 0.02;

field_choice = -653;

delta_phi = 0. * pi / 4.;

% filename = "eigenMode.mat";

% Wave period
switch field_choice
    % ============= KAWs ================
    case {-41, -45, -47, -451}
        bi = 1.;
        waveT = 110.72;
        kvalue = 0.05;
    case {-49, -491}
        bi = 1.;
        waveT = 110.72;
        kvalue = -0.05;
    case {-495, -4951}
        bi = 1.;
        waveT = 110.72;
        kvalue = [-0.05, 0.05];
    case -43
        bi = 1.;
        waveT = 55.61;
        kvalue = 0.1;
    % ============= ICWs ================
    case {-51, -53, -511}
        bi = 1.;
        waveT = 20.1657;
        kvalue = 0.8;
    case -55
        bi = 1.;
        waveT = 25.78599;
        kvalue = 0.4;
    case -59
        bi = 1.;
        waveT = 22.44201;
        kvalue = 0.6;
    case -61
        bi = 1.;
        waveT = 23.94790;
        kvalue = 0.5;
    case {-63, -6300}
        bi = 1.;
        waveT = 23.57534;
        kvalue = 0.525;
    case -65
        bi = 1.;
        waveT = 23.57534;
        kvalue = -0.525;
    case {-653, -6531, -6532, -6533}
        bi = 1.;
        waveT = 23.57534;
        kvalue = [-0.525, 0.525];
    case -6301
        bi = 0.1;
        waveT = 10.61342774;
        kvalue = 0.525;
    case -6303
        bi = 0.3;
        waveT = 14.52614851;
        kvalue = 0.525;        
    case -633
        bi = 3.;
        waveT = 41.67231726;
        kvalue = 0.525;         
    case -6310
        bi = 10.;
        waveT = 90.66208508;
        kvalue = 0.525; 
end

% Physical Initial Time
t_init = - 4. * waveT;

% Plasma beta_i and beta_e
be = bi * mime / tite;

% EM field choice
% -35: ExB Drift Fields
% -----------------------------------------------------
% -37: single KAW, 
%      from PLUME calculation,
%      betap = 1, kperp = 1, kpar = 1e-3, mi/me = 1836
%      period: T = 2pi/0.001137 = 5526.11, unit: t \Omega_0
% -----------------------------------------------------
% -39: single KAW, TBF
% -----------------------------------------------------
% -41: single KAW,
%      from PLUME calculation,
%      betap = 1, kperp = 1, kpar = 0.05, mi/me = 1836
%      large B0, no re-scale to PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -43: single KAW,
%      from PLUME calculation,
%      betap = 1, kperp = 1, kpar = 0.1, mi/me = 1836
%      period: T = 2pi/1.129786e-01 = 55.61, unit: t \Omega_0
% -----------------------------------------------------
% -45: single KAW,
%      from PLUME calculation,
%      betap = 1, kperp = 1, kpar = 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -451: single KAW, with window function
%      from PLUME calculation,
%      betap = 1, kperp = 1, kpar = 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -47: single KAW,
%      from PLUME calculation, only Epar is retained
%      betap = 1, kperp = 1, kpar = 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -49: single KAW,
%      from PLUME calculation, Negative kpar
%      betap = 1, kperp = 1, kpar = - 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -491: single KAW, with window function
%      from PLUME calculation, Negative kpar
%      betap = 1, kperp = 1, kpar = - 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -495: two KAWs, no window function
%      from PLUME calculation, Positive and Negative kpar
%      betap = 1, kperp = 1, kpar = \mp 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -4951: two KAWs, with window function
%      from PLUME calculation, Positive and Negative kpar
%      betap = 1, kperp = 1, kpar = \mp 0.05, mi/me = 1836
%      re-scaled PLUME outputs
%      period: T = 2pi/0.056748 = 110.72, unit: t \Omega_0
% -----------------------------------------------------
% -51: single ICW,
%      from PLUME calculation,
%      betap = 1, kperp = 0.01, kpar = 0.8, mi/me = 1836
%      no time re-scale
%      re-scaled PLUME outputs, with Window Function
%      period: T = 2pi/0.3116 = 20.1657, unit: t \Omega_0
% -----------------------------------------------------
% -53: single ICW,
%      from PLUME calculation,
%      betap = 1, kperp = 0.01, kpar = 0.8, mi/me = 1836
%      no time re-scale
%      re-scaled PLUME outputs, No Window Function
%      period: T = 2pi/0.3116 = 20.1657, unit: t \Omega_0
% -----------------------------------------------------
% -511: single ICW, 
%       same field as -51
%       verify my reading data from file code
% -----------------------------------------------------
% -55: single ICW, kperp = 0.01, kpar = 0.4, betap = 1
%      gamma/omega = -0.088, T = 25.785993267766635
% -----------------------------------------------------
% -57: single ICW, kperp = 0.01, kpar = XXX, betap = 1
% -----------------------------------------------------
% -59: single ICW, kperp = 0.01, kpar = 0.6, betap = 1
%      gamma/omega = -0.45132964, T = 22.442007381318877
% -----------------------------------------------------
% -61: single ICW, kperp = 0.01, kpar = 0.5, betap = 1
%      gamma/omega = -0.25421113, T = 23.947895167415304
% -----------------------------------------------------
% -63: single ICW, kperp = 0.01, kpar = 0.525, betap = 1
%      gamma/omega = -0.29831079, T = 23.57533569259684
% -----------------------------------------------------
% -6300: No Win, single ICW, kperp = 0.01, kpar = 0.525, betap = 1
%      gamma/omega = -0.29831079, T = 23.57533569259684
% -----------------------------------------------------
% -65: single ICW, kperp = 0.01, kpar = -0.525, betap = 1
%      gamma/omega = -0.29831079, T = 23.57533569259684
% -----------------------------------------------------
% -653: two ICWs, kperp = 0.01, kpar = \mp 0.525, betap = 1
%      gamma/omega = -0.48779923, T = 23.57533569259684
% -----------------------------------------------------
% -6531: Non-Physical, two ICWs, kperp = 0.01, kpar = \mp 0.525, betap = 1
%      gamma/omega = -0.29831079, T = 23.57533569259684
%      test shut down random EM field components
% -----------------------------------------------------
% -6532: two ICWs, kperp = 0.01, kpar = \mp 0.525, betap = 1
%      gamma/omega = -0.29831079, T = 23.57533569259684
%      test uneven wave amplitude, 2 kpar_pos; 1 kpar_neg
% -----------------------------------------------------
% -6301: single ICW, kperp = 0.01, kpar = 0.525, betap = 0.1
%      gamma/omega = -0.48779923, T = 10.61342774
% -----------------------------------------------------
% -6303: single ICW, kperp = 0.01, kpar = 0.525, betap = 0.3
%      gamma/omega = -0.37363259, T = 14.52614851
% -----------------------------------------------------
% -633: single ICW, kperp = 0.01, kpar = 0.525, betap = 3
%      gamma/omega = -0.27533507, T = 41.67231726
% -----------------------------------------------------
% -6310: single ICW, kperp = 0.01, kpar = 0.525, betap = 10
%      gamma/omega = -0.29255560, T = 90.66208508

