function dy = lorentz4(t,y)
% LORENTZ Calculates particle motion in given E and B fields according to
%    dimensionless equations for single particle motion
% Dimensionless equations:
%    dx'/dt' =v'
%    dv'/dt'= E' + v' x B'
% Input values are
%      vector y(6)= (x,y,z,vx,vy,vz) (Column vector)(1:6)
%      t=time (scalar)
%      q'= q/e
%      m'= m/m_e
% The position y(1:3) and time t are used to obtain E and B at particle position
%     using functions magfield and elecfield
% Output values are
%      dy(6)=(dx/dt,dy/dt,dz/dt,dvx/dt,dvy/dt,dvz/dt) (Column vector)(1:6)

% Key parameters to define Plasma

  % persistent initialized
  % persistent q m mime tite bi be vtic field_choice em_eps waveT t_init filename

  % if isempty(initialized)
    [q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, delta_phi, kvalue] = set_params;
    % initialized = 1;
  % end
% t
% Get electric and magnetic fields
E0 = elecfield(t, y(1:3));
B0 = magfield(t, y(1:3));
  
% Initialize dy solution as column vector
dy = zeros(6,1) ;  %Column vector with 6 values

% follow the form of SPM dimensionless equations
dy(1:3) = y(4:6);
dy(4) = q / m * (E0(1) + y(5) * B0(3) - y(6) * B0(2));
dy(5) = q / m * (E0(2) + y(6) * B0(1) - y(4) * B0(3));
dy(6) = q / m * (E0(3) + y(4) * B0(2) - y(5) * B0(1));


