function B = magfield(t, x, choice)
% MAGNETIC Calculates magnetic field B(t,x)
% Input values are scalar time t and vector position x(3)= (x,y,z) (de normalized)
% Output values are column vector B(3)=(Bx,By,Bz)
%
% FIELD Geometry choice (geom):
% See set_params.m

  % persistent initialized initialized_fields eigenMode
  % persistent q m mime tite bi be vtic field_choice em_eps waveT t_init filename

  % if isempty(initialized)
    [q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, delta_phi, kvalue] = set_params;
    % initialized = 1;
  % end

  % Check if `choice` argument is provided
  if nargin == 3
    geom = choice;
  else
    geom = field_choice;
  end

% Initialize column vector
  B = zeros(3,1);

switch(geom)
%=====================================================================
    case -35
        % testing case
        B(1) = 0.;
        B(2) = 0.;
        B(3) = 1.;
    case -37
        k = [1.; 0.; 0.001];
        B_tilde_raw = vtic * [0.252270 - 4.627666i; 10767.230000 - 163.727200i; -252.270200 + 4627.666000i];
        % omega = 0.001137 - 0.000033i;
        B_tilde = em_eps * B_tilde_raw;
        omega = 0.001137;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B(3) = B(3) + abs(B_tilde_raw(2)) ;

    case -41
        k = [1.; 0.; 0.05];
        B_tilde_raw = vtic * [12.318400 - 231.109100i; 10794.910000 - 169.658700i; -246.368000 + 4622.182000i];
        B_tilde = 1 * B_tilde_raw;
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B(3) = B(3) + abs(B_tilde_raw(2)) / em_eps;

    case -43
        k = [1.; 0.; 0.1];
        B_tilde = [2.282905e+01 - 4.605594e+02i; 1.088140e+04 - 1.880723e+02i; -2.282905e+02 + 4.605594e+03i];
        B_tilde = vtic * B_tilde;
        % omega = 1.129786e-01 - 3.214358e-03i;
        omega = 1.129786e-01;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B(3) = B(3) + vtic * sqrt(1.088140e+04^2 + 1.880723e+02^2) * B0By;

    case -45 % no window function
        k = [1.; 0.; 0.05];
        B_tilde_raw = vtic * [12.318400 - 231.109100i; 10794.910000 - 169.658700i; -246.368000 + 4622.182000i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        % B(3) = B(3) + abs(B_tilde_raw(2));
        B(3) = B(3) + 1.;

    case -451 % single KAW positive kpar with window function
        k = [1.; 0.; 0.05];
        B_tilde_raw = vtic * [12.318400 - 231.109100i; 10794.910000 - 169.658700i; -246.368000 + 4622.182000i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        % B(3) = B(3) + abs(B_tilde_raw(2));
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -47
        k = [1.; 0.; 0.05];
        B_tilde_raw = vtic * [12.318400 - 231.109100i; 10794.910000 - 169.658700i; -246.368000 + 4622.182000i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        % B(3) = B(3) + abs(B_tilde_raw(2));
        B(3) = B(3) + 1.;
        B(1) = 0.;
        B(2) = 0.;

    case -49 % No window function
        delta_phi = 0.;
        k = [1.; 0.; -0.05];
        B_tilde_raw = vtic * [-12.318400 + 231.109100i; -10794.910000 + 169.658700i; -246.368000 + 4622.182000i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        % B(3) = B(3) + abs(B_tilde_raw(2));
        B(3) = B(3) + 1.;

    case -491 % single KAW negative kpar with window function
        delta_phi = 0.;
        k = [1.; 0.; -0.05];
        B_tilde_raw = vtic * [-12.318400 + 231.109100i; -10794.910000 + 169.658700i; -246.368000 + 4622.182000i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.056748 - 0.001630i;
        omega = 0.056748;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        
        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -495
        k1 = [1.; 0.; -0.05];
        B_tilde_raw1 = vtic * [-12.318400 + 231.109100i; -10794.910000 + 169.658700i; -246.368000 + 4622.182000i];
        B_tilde1 = em_eps * B_tilde_raw1;
        % omega = 0.056748 - 0.001630i;
        omega1 = 0.056748;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % B1 = 0.5 * (B_tilde1 * exp_term1 + conj(B_tilde1) * conj(exp_term1));
        B1 = B_tilde1 * exp_term1;

        k2 = [1.; 0.; 0.05];
        B_tilde_raw2 = vtic * [12.318400 - 231.109100i; 10794.910000 - 169.658700i; -246.368000 + 4622.182000i];
        B_tilde2 = em_eps * B_tilde_raw2;
        % omega = 0.056748 - 0.001630i;
        omega2 = 0.056748;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % B2 = 0.5 * (B_tilde2 * exp_term2 + conj(B_tilde2) * conj(exp_term2));
        B2 = B_tilde2 * exp_term2;
        
        B = B1 + B2;
        B(3) = B(3) + 1.;
        B = real(B);

    case -4951
        k1 = [1.; 0.; -0.05];
        B_tilde_raw1 = vtic * [-12.318400 + 231.109100i; -10794.910000 + 169.658700i; -246.368000 + 4622.182000i];
        B_tilde1 = em_eps * B_tilde_raw1;
        % omega = 0.056748 - 0.001630i;
        omega1 = 0.056748;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % B1 = 0.5 * (B_tilde1 * exp_term1 + conj(B_tilde1) * conj(exp_term1));
        B1 = B_tilde1 * exp_term1;

        k2 = [1.; 0.; 0.05];
        B_tilde_raw2 = vtic * [12.318400 - 231.109100i; 10794.910000 - 169.658700i; -246.368000 + 4622.182000i];
        B_tilde2 = em_eps * B_tilde_raw2;
        % omega = 0.056748 - 0.001630i;
        omega2 = 0.056748;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % B2 = 0.5 * (B_tilde2 * exp_term2 + conj(B_tilde2) * conj(exp_term2));
        B2 = B_tilde2 * exp_term2;
        
        B = B1 + B2;

        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end
        B = smooth_factor * B;
        B(3) = B(3) + 1.;
        B = real(B);

     case -51
        k = [0.01; 0.; 0.8];
        B_tilde_raw = vtic * [-1.278282e+04 + 1.459139e+04 * 1i; 1.459472e+04 + 1.278179e+04 * 1i; 1.590841e+02 - 1.815920e+02 * 1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 3.115772e-01 - 2.729087e-01 * 1i;
        omega = 3.115772e-01;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));
        dt = t - t_init;

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

     case -53
        k = [0.01; 0.; 0.8];
        B_tilde_raw = vtic * [-1.278282e+04 + 1.459139e+04i; 1.459472e+04 + 1.278179e+04i; 1.590841e+02 - 1.815920e+02i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 3.115772e-01 - 2.729087e-01i;
        omega = 3.115772e-01;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B(3) = B(3) + 1.;
    
    case -511
        if isempty(initialized_fields)
            load(filename);
            initialized_fields = 1;
        end

        k = [eigenMode.kperp; 0.; eigenMode.kpar];
        B_tilde_raw = vtic * [eigenMode.bxr + eigenMode.bxi*1i; eigenMode.byr + eigenMode.byi*1i; eigenMode.bzr + eigenMode.bzi*1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 3.115772e-01 - 2.729087e-01i;
        omega = eigenMode.w;
        k_dot_x = dot(k, x);

        exp_term = exp(1i * (k_dot_x - omega * t));
        dt = t - t_init;

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

     case -55
        k = [0.01; 0.; 0.4];
        B_tilde_raw = vtic * [-1441.113 + 16194.18 * 1i; 16213.37 + 1427.344 * 1i; 36.19913 - 406.7795 * 1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.2436666 - 0.02148634 * 1i;
        omega = 0.2436666;
        k_dot_x = dot(k, x);
        
        dt = t - t_init;
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -59
        k = [0.01; 0.; 0.6];
        B_tilde_raw = vtic * [-8072.301 + 17874.35 * 1i; 17881.59 + 8068.456 * 1i; 133.9669 - 296.6403 * 1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.27997430 - 0.12636070 *1i;
        omega = 0.27997430;
        k_dot_x = dot(k, x);
        
        dt = t - t_init;
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -61
        k = [0.01; 0.; 0.5];
        B_tilde_raw = vtic * [-4566.727 + 17933.15 * 1i; 17944.41 + 4559.459 * 1i; 91.11818 - 357.8133 * 1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.26236900 - 0.06669712 *1i;
        omega = 0.26236900;
        k_dot_x = dot(k, x);
        
        dt = t - t_init;
        exp_term = exp(1i * (k_dot_x - omega * t));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -63
        k = [0.01; 0.; 0.525];
        B_tilde_raw = vtic * [-5398.741 + 18073.83 * 1i; 18083.95 + 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega = 0.26651520;
        k_dot_x = dot(k, x);
        
        dt = t - t_init;
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -6300
        k = [0.01; 0.; 0.525];
        B_tilde_raw = vtic * [-5398.741 + 18073.83 * 1i; 18083.95 + 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde = em_eps * B_tilde_raw;
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

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -65
        delta_phi = 0.;
        k = [0.01; 0.; -0.525];
        B_tilde_raw = vtic * [5398.741 - 18073.83 * 1i; -18083.95 - 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde = em_eps * B_tilde_raw;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega = 0.26651520;
        k_dot_x = dot(k, x);
        
        dt = t - t_init;
        exp_term = exp(1i * (k_dot_x - omega * t + delta_phi));

        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -653
        k1 = [0.01; 0.; -0.525];
        B_tilde_raw1 = vtic * [5398.741 - 18073.83 * 1i; -18083.95 - 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde1 = em_eps * B_tilde_raw1;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % B1 = 0.5 * (B_tilde1 * exp_term1 + conj(B_tilde1) * conj(exp_term1));
        B1 = B_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        B_tilde_raw2 = vtic * [-5398.741 + 18073.83 * 1i; 18083.95 + 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde2 = em_eps * B_tilde_raw2;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % B2 = 0.5 * (B_tilde2 * exp_term2 + conj(B_tilde2) * conj(exp_term2));
        B2 = B_tilde2 * exp_term2;
        
        B = B1 + B2;

        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = smooth_factor * B;
        B(3) = B(3) + 1.;
        B = real(B);

    case -6531 % !!!WARNING: Non-physical fields
        k1 = [0.01; 0.; -0.525];
        B_tilde_raw1 = vtic * [5398.741 - 18073.83 * 1i; -18083.95 - 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde1 = em_eps * B_tilde_raw1;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % B1 = 0.5 * (B_tilde1 * exp_term1 + conj(B_tilde1) * conj(exp_term1));
        B1 = B_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        B_tilde_raw2 = vtic * [-5398.741 + 18073.83 * 1i; 18083.95 + 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde2 = em_eps * B_tilde_raw2;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % B2 = 0.5 * (B_tilde2 * exp_term2 + conj(B_tilde2) * conj(exp_term2));
        B2 = B_tilde2 * exp_term2;
        
        B = B1 + B2;

        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = smooth_factor * B;
        B(3) = B(3) + 1.;
        B = real(B);
    
    case -6532 % pos k mode has an amplitude twice of that of neg k mode
        k1 = [0.01; 0.; -0.525];
        B_tilde_raw1 = vtic * [5398.741 - 18073.83 * 1i; -18083.95 - 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde1 = em_eps * B_tilde_raw1;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % B1 = 0.5 * (B_tilde1 * exp_term1 + conj(B_tilde1) * conj(exp_term1));
        B1 = B_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        B_tilde_raw2 = vtic * [-5398.741 + 18073.83 * 1i; 18083.95 + 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde2 = 2 * em_eps * B_tilde_raw2;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % B2 = 0.5 * (B_tilde2 * exp_term2 + conj(B_tilde2) * conj(exp_term2));
        B2 = B_tilde2 * exp_term2;
        
        B = B1 + B2;

        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = smooth_factor * B;
        B(3) = B(3) + 1.;
        B = real(B);

    case -6533 % neg k mode has an amplitude twice of that of pos k mode
        k1 = [0.01; 0.; -0.525];
        B_tilde_raw1 = vtic * [5398.741 - 18073.83 * 1i; -18083.95 - 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde1 = 2 * em_eps * B_tilde_raw1;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega1 = 0.26651520;
        k_dot_x1 = dot(k1, x);

        exp_term1 = exp(1i * (k_dot_x1 - omega1 * t));

        % B1 = 0.5 * (B_tilde1 * exp_term1 + conj(B_tilde1) * conj(exp_term1));
        B1 = B_tilde1 * exp_term1;

        k2 = [0.01; 0.; 0.525];
        B_tilde_raw2 = vtic * [-5398.741 + 18073.83 * 1i; 18083.95 + 5392.464 * 1i; 102.8709 - 344.3898 * 1i];
        B_tilde2 = em_eps * B_tilde_raw2;
        % omega = 0.26651520 - 0.07950436 *1i;
        omega2 = 0.26651520;
        k_dot_x2 = dot(k2, x);

        exp_term2 = exp(1i * (k_dot_x2 - omega2 * t + delta_phi));

        % B2 = 0.5 * (B_tilde2 * exp_term2 + conj(B_tilde2) * conj(exp_term2));
        B2 = B_tilde2 * exp_term2;
        
        B = B1 + B2;

        dt = t - t_init;
        if dt < waveT
            smooth_factor = 0.5 * (1 - cos(pi * dt / waveT));
            % smooth_factor = 1.;
        else
            smooth_factor = 1.;
        end

        B = smooth_factor * B;
        B(3) = B(3) + 1.;
        B = real(B);

    case -6301
        k = [0.01; 0.; 0.525];
        B_tilde_raw = vtic * [-3490.462 + 7159.495 * 1i; 7162.877 + 3490.106 * 1i; 66.50938 - 136.4214 * 1i];
        B_tilde = em_eps * B_tilde_raw;
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

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;
    
    case -6303
        k = [0.01; 0.; 0.525];
        B_tilde_raw = vtic * [-3977.924 + 10642.36 * 1i; 10648.34 + 3976.471 * 1i; 75.79778 - 202.786 * 1i];
        B_tilde = em_eps * B_tilde_raw;
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

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -633
        k = [0.01; 0.; 0.525];
        B_tilde_raw = vtic * [-8922.261 + 32340.59 * 1i; 32356.56 + 8905.672 * 1i; 170.0102 - 616.2372 * 1i];
        B_tilde = em_eps * B_tilde_raw;
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

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;

    case -6310
        k = [0.01; 0.; 0.525];
        B_tilde_raw = vtic * [-20443.66 + 69729.53 * 1i; 69760.07 + 20402.26 * 1i; 389.5459 - 1328.669 * 1i];
        B_tilde = em_eps * B_tilde_raw;
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

        B = 0.5 * (B_tilde * exp_term + conj(B_tilde) * conj(exp_term));
        B = smooth_factor * B;
        B(3) = B(3) + 1.;
end
