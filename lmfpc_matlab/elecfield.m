function E = elecfield(t, x, choice)
% ELECTRIC Calculates electric field E(t,x)
% Input values are scalar time t and vector position x(3)= (x,y,z)
% Output values are column vector E(3)=(Ex,Ey,Ez)

% FIELD Geometry choice (egeom):
% See set_params.m

% Persistent variables for data interpolation (egeom=-511)
  % persistent initialized initialized_fields eigenMode
  % persistent q m mime tite bi be vtic field_choice em_eps waveT t_init filename

  % if isempty(initialized)
    [q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, delta_phi, kpardi] = set_params;
    % initialized = 1;
  % end

  % Check if `choice` argument is provided
  if nargin == 3
    egeom = choice;
  else
    egeom = field_choice;
  end

% Initialize column vector
  E=zeros(3,1);
			
switch(egeom)
%=====================================================================
    case -35
        E = [0.; 0.1; 0.];
    case -37
        k = [1.; 0.; 0.001];
        E_tilde = em_eps * [1. + 0.i; -0.000014 + 0.000527i; -0.000223 + 0.000054i];
        % omega = 0.001137 - 0.000033i;
        omega = 0.001137;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
 
    case -41
        k = [1.; 0.; 0.05];
        E_tilde = [1. + 0.i; -0.000644 + 0.026270i; -0.011231 + 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
    
    case -43
        k = [1.; 0.; 0.1];
        E_tilde = [1. + 0.i; -1.098791e-03 + 5.210673e-02i; -2.287604e-02 + 5.622485e-03i];
        % omega = 1.129786e-01 - 3.214358e-03i;
        omega = 1.129786e-01;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

    case -45 % no window function
        k = [1.; 0.; 0.05];
        E_tilde = em_eps * [1. + 0.i; -0.000644 + 0.026270i; -0.011231 + 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);
        
        % delta_phi = 0;
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

    case -451 % single KAW positive kpar with window function
        k = [1.; 0.; 0.05];
        E_tilde = em_eps * [1. + 0.i; -0.000644 + 0.026270i; -0.011231 + 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
    
    case -47
        k = [1.; 0.; 0.05];
        E_tilde = em_eps * [1. + 0.i; -0.000644 + 0.026270i; -0.011231 + 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
        E(1) = 0.;
        E(2) = 0.;

    case -49
        delta_phi = 0.;
        k = [1.; 0.; -0.05];
        E_tilde = em_eps * [1. + 0.i; -0.000644 + 0.026270i; 0.011231 - 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);
        
        % delta_phi = 0;
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));


    case -491
        delta_phi = 0.;
        k = [1.; 0.; -0.05];
        E_tilde = em_eps * [1. + 0.i; -0.000644 + 0.026270i; 0.011231 - 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);
        
        % delta_phi = 0;
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
    
    case -495
        k1 = [1.; 0.; -0.05];
        E_tilde1 = em_eps * [1. + 0.i; -0.000644 + 0.026270i; 0.011231 - 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega1 = 0.056748;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        E1 = E_tilde1 * exp_term1;

        k2 = [1.; 0.; 0.05];
        E_tilde2 = em_eps * [1. + 0.i; -0.000644 + 0.026270i; -0.011231 + 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega2 = 0.056748;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        E2 = E_tilde2 * exp_term2;

        E = E1 + E2;
        E = real(E);
    
    
    case -4951
        k1 = [1.; 0.; -0.05];
        E_tilde1 = em_eps * [1. + 0.i; -0.000644 + 0.026270i; 0.011231 - 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega1 = 0.056748;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % E1 = 0.5 * (E_tilde1 * exp_term1 + conj(E_tilde1) * conj(exp_term1));
        E1 = E_tilde1 * exp_term1;

        k2 = [1.; 0.; 0.05];
        E_tilde2 = em_eps * [1. + 0.i; -0.000644 + 0.026270i; -0.011231 + 0.002723i];
        % omega = 0.056748 - 0.001630i;
        omega2 = 0.056748;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % E2 = 0.5 * (E_tilde2 * exp_term2 + conj(E_tilde2) * conj(exp_term2));
        E2 = E_tilde2 * exp_term2;

        E = E1 + E2;
    
        dt = t - t_init;

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
        E = real(E);
    
    case -51
        k = [0.01; 0.; 0.8];
        E_tilde = em_eps * [1. + 0.i; 8.942908e-05 - 9.999534e-01i; -3.822730e-03 + 5.121012e-03i];
        % omega = 3.115772e-01 - 2.729087e-01i;
        omega = 3.115772e-01;
        k_dot_x = dot(k, x);

        % period = 2. * pi / omega;
        dt = t - t_init;
        
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

    case -53 % no window function
        k = [0.01; 0.; 0.8];
        E_tilde = em_eps * [1. + 0.i; 8.942908e-05 - 9.999534e-01i; -3.822730e-03 + 5.121012e-03i];
        % omega = 3.115772e-01 - 2.729087e-01i;
        omega = 3.115772e-01;
        k_dot_x = dot(k, x);
        
        exp_term = exp(1i * (k_dot_x - omega * t));
        
        E = 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

    case -511
        if isempty(initialized_fields)
            load(filename);
            initialized_fields = 1;
        end
        
        k = [eigenMode.kperp; 0.; eigenMode.kpar];
        E_tilde = em_eps * [eigenMode.exr + eigenMode.exi*1i; eigenMode.eyr + eigenMode.eyi*1i; eigenMode.ezr + eigenMode.ezi*1i];
        % omega = 3.115772e-01 - 2.729087e-01i;
        omega = eigenMode.w;

        k_dot_x = dot(k, x);
        
        dt = t - t_init;
        
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));        
    case -55
        k = [0.01; 0.; 0.4];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00080315 - 0.9989638 * 1i; -0.00253034 + 0.00569913 * 1i];
        % omega = 0.2436666 - 0.02148634 *1i;
        omega = 0.2436666;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
     case -59
        k = [0.01; 0.; 0.6];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00023585 - 0.99979840 * 1i; -0.00325777 + 0.00569910 * 1i];
        % omega = 0.27997430 - 0.12636070 *1i;
        omega = 0.27997430;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
     case -61
        k = [0.01; 0.; 0.5];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00041470 - 0.99956460 * 1i; -0.00286969 + 0.00579388 * 1i];
        % omega = 0.26236900 - 0.06669712 *1i;
        omega = 0.26236900;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

     case -63
        k = [0.01; 0.; 0.525];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; -0.00296036 + 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega = 0.26651520;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
     
    case -6300
        k = [0.01; 0.; 0.525];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; -0.00296036 + 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega = 0.26651520;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            % smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
     
    case -65
        delta_phi = 0.;
        k = [0.01; 0.; -0.525];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; 0.00296036 - 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega = 0.26651520;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
    
    case -653
        k1 = [0.01; 0.; -0.525];
        E_tilde1 = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; 0.00296036 - 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % E1 = 0.5 * (E_tilde1 * exp_term1 + conj(E_tilde1) * conj(exp_term1));
        E1 = E_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        E_tilde2 = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; -0.00296036 + 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % E2 = 0.5 * (E_tilde2 * exp_term2 + conj(E_tilde2) * conj(exp_term2));
        E2 = E_tilde2 * exp_term2;

        E = E1 + E2;
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
        
        E = real(E);

    case -6531 % !!!WARNING: Non-physical fields
        k1 = [0.01; 0.; -0.525];
        E_tilde1 = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; 0.00296036 - 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % E1 = 0.5 * (E_tilde1 * exp_term1 + conj(E_tilde1) * conj(exp_term1));
        E1 = E_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        E_tilde2 = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; -0.00296036 + 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % E2 = 0.5 * (E_tilde2 * exp_term2 + conj(E_tilde2) * conj(exp_term2));
        E2 = E_tilde2 * exp_term2;

        E = E1 + E2;
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
        E = real(E);
        E(3) = 0.;

    case -6532 % pos k mode has an amplitude twice of that of neg k mode
        k1 = [0.01; 0.; -0.525];
        E_tilde1 = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; 0.00296036 - 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % E1 = 0.5 * (E_tilde1 * exp_term1 + conj(E_tilde1) * conj(exp_term1));
        E1 = E_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        E_tilde2 = 2 * em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; -0.00296036 + 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % E2 = 0.5 * (E_tilde2 * exp_term2 + conj(E_tilde2) * conj(exp_term2));
        E2 = E_tilde2 * exp_term2;

        E = E1 + E2;
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
        
        E = real(E);

    case -6533 % neg k mode has an amplitude twice of that of pos k mode
        k1 = [0.01; 0.; -0.525];
        E_tilde1 = 2 * em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; 0.00296036 - 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % E1 = 0.5 * (E_tilde1 * exp_term1 + conj(E_tilde1) * conj(exp_term1));
        E1 = E_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        E_tilde2 = em_eps * [1. + 0. * 1i; 0.00036167 - 0.99963800 * 1i; -0.00296036 + 0.00578877 * 1i];
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % E2 = 0.5 * (E_tilde2 * exp_term2 + conj(E_tilde2) * conj(exp_term2));
        E2 = E_tilde2 * exp_term2;

        E = E1 + E2;
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        E = smooth_factor * E;
        
        E = real(E);

     case -6301
        k = [0.01; 0.; 0.525];
        E_tilde = em_eps * [1. + 0. * 1i; -0.00021820 - 0.99968420 * 1i; -0.00241984 + 0.02332500 * 1i];
        % omega = 0.59200340 - 0.28877880 *1i;
        omega = 0.59200340;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));
    
    case -6303
        k = [0.01; 0.; 0.525];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00013122 - 0.99963530 * 1i; -0.00438972 + 0.00906996 * 1i];
        % omega = 0.43254310 - 0.16161220 *1i;
        omega = 0.43254310;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

    case -633
        k = [0.01; 0.; 0.525];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00051026 - 0.99971560 * 1i; -0.00227139 + 0.00486228 * 1i];
        % omega = 0.15077600 - 0.04151392 *1i;
        omega = 0.15077600;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

    case -6310
        k = [0.01; 0.; 0.525];
        E_tilde = em_eps * [1. + 0. * 1i; 0.00057961 - 0.99979260 * 1i; -0.00188862 + 0.00446625 * 1i];
        % omega = 0.06930334 - 0.02027508 *1i;
        omega = 0.06930334;
        k_dot_x = dot(k, x);

        dt = t - t_init;      
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end
        
        E = smooth_factor * 0.5 * (E_tilde * exp_term + conj(E_tilde) * conj(exp_term));

end
