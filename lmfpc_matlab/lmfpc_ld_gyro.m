% Determine Distribution Using Liouville-Mapping: ion Landau damping
% SPM: Driver for Single Particle Motion in given E and B fields
% Dimensionless equations:
%   dx'/dt' = v'
%   dv'/dt'= E' + v' x B'
% Dimensionless choice (CGS units):
%   x' = x / rho_i
%   v' = v / v_ti
%   t' = Omega_i * t
%   E' = E / (v_ti * B0 / c)
%   B' = B / B0
% PLUME Normalizations (CGS units):
%   E~ = E / (v_ti * B0 / c)
%   B~ = B / (v_ti * B0 / c)
%   vtp = v_ti / c = 1E-4
% To apply PLUME outputs to this iSH SPM Solver
%   Use E~ directly
%   B' = B~ * (v_ti / c)
% Numerical Method choices 
%    Adaptive Runge-Kutta 45 method with RelTol specified
% Input x0, v0, B0, E0,
% Compute x and v as a function of time

close all
clear
clc

format long

disp('=======================================================================');
disp("You've started lmfpc_ld_gyro.");
% This version only calculates f, df/dz, cez
% This version starts with vperp and theta grids.

% Key parameters to define Plasma and turbulence
[q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, filename, delta_phi, kpardi] = set_params;
field_choice;
em_eps;
waveT;
t_init;

% Setting Position to Compute f(vx,vy)----------------------
% Set x limits
% xval = (0.1:-1.0:0.1)';
xval = 0.1;
yval = 0.1;
zval = 0.1;
nx = length(xval);

% Create (vperp, theta, vz) grid ================================================
% Create Velocity Space Grid (v_ti units)

vlength = 3.;
nvperp = 9; % Number of sampled points along vperp axis
ntheta = 1; % Number of sampled theta
nvz = 62; % Number of sampled points along vz axis

% Compute vmin_perp, vmax_perp, dv_perp
dv_perp = vlength / nvperp;
vmin_perp = 0.;
vmax_perp = vlength;

% Compute dtheta
dtheta = 2. * pi / ntheta;

% Compute vmin_z, vmax_z, dv_z
dv_z = 2. * vlength / nvz;
vmin_z = -1. * vlength - dv_z / 2.;
vmax_z = vlength + dv_z / 2.;

% Set up (vperp, theta, vz) grid
vperp_values = (vmin_perp:dv_perp:vmax_perp)
% theta_values = (0+dtheta:dtheta:2.*pi)
theta_values = 2. * pi
vz_values = (vmin_z:dv_z:vmax_z)

nvperp = length(vperp_values);
ntheta = length(theta_values);
nvz = length(vz_values);

% Set Time grids
% ti = -333.; % has been specified in set_params
dt_final = waveT / 40.;
t_final = (0.:dt_final:1.*waveT)' % The correlation interval
nt = length(t_final);

% END Create (vperp, theta, vz) grid ================================================

% Initialize f, dfdvz, cez, cez_PerpPar
f = zeros(nx, nvperp, ntheta, nvz, nt);
dfdvz = zeros(nvperp, ntheta, nvz);
cez = zeros(nx, nvperp, ntheta, nvz, nt);
cez_PerpPar = zeros(nx, nvperp, ntheta, nvz, nt);

% Set tolerance
options = odeset('RelTol', 1.0e-3);

% Some parameters that will be used later
adjust_lim = 0;
vperp_cutoff = 0;

% Set function handle to Lorentz Force Law function
lorhandle = @lorentz4;

% Get the current date and time for filename suffix
current_time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
time_suffix = char(current_time);

save_figure = true; % if true, save figures as png files
fpcfigure_filename = sprintf('./plots/iSHLDV12_%s.png', time_suffix);

save_mat_data = true; % if true, save all workspace 
mat_data_name = sprintf("./data/iSHLDV12_%s.mat", time_suffix);

% pause
disp("Initiation done.");
disp("Now start the numerical calculation process for each velocity space grid point.")

% LOOP OVER x-positions for distribution function calculation
% This is the core code of this script.
% It calculates f(xval, vx, vy, vz), df/dv, FPC, j dot E, reduced FPC
for it = 1:nt
    fprintf('t_final= %3.2f\n', t_final(it));
    for ix = 1:nx %________________________________________________________________
        % Loop over each point in the (vperp, theta, vz) grid
        for ivperp = 1:nvperp
            tic; % Start timing the time loop
            for itheta = 1:ntheta
                for ivz = 1:nvz
                    vx0 = vperp_values(ivperp) * cos(theta_values(itheta));
                    vy0 = vperp_values(ivperp) * sin(theta_values(itheta));
                    vz0 = vz_values(ivz);
                    x0 = [xval(ix), yval(ix), zval(ix)]';  % Initial position (Column vector)
                    v0 = [vx0, vy0, vz0]';  % Initial velocity (Column vector)
                    % Convert x0, v0 to y0(1:6) initial conditions
                    y0 = [x0', v0']; % Row vector
                    % I may declutter these two transposes later.
        
                    % Calling ode45 for Runge-Kutta with option for setting error tolerance
                    % [t4, y4] is the solution set
                    % t4: the time array, a column vector, each element
                    %     represents a time slice
                    % y4: the solution matrix, contains 6 columns x, y, z, vx,
                    %     vy, vz, each row corresponds to one time slice
                    % nsteps: how many timeslices we have
                    tspan = [t_final(it), t_init]';
                    [t4,y4] = ode45(lorhandle, tspan, y0, options);
                    nsteps = size(t4,1);
        
                    % y4(end) gives x and v at the "initial" time (because of integrating backwards)
                    % here, acrodding to Liouville's theorem, the f value calculated from 
                    % the "initial" time, which is a Maxwellian, has been assigned to 
                    % the phase-space grid we want, at (x, vx, vy, vz)
                    % need to change the expression of f for iSH
    
                    f(ix, ivperp, itheta, ivz, it) = exp(- y4(end,4)^2. - y4(end,5)^2.- y4(end,6)^2.);
        
                end % vz loop 
            end  % theta loop

            fprintf('vperp = %3.2f theta = %3.2f vz =  %3.2f \n', vperp_values(ivperp), theta_values(itheta),  vz_values(ivz));
            % Stop timing and output the computation time
            elapsedTime = toc;
            fprintf('Time taken for computing f at this vperp loop is\n %f seconds\n', elapsedTime);

            % if (mod(ivperp,4)==0)
            %     [vperp_values(ivperp)] % This output serves as a progress bar.
            % end
    
        end % vperp loop
    
        disp('Numerical calculation at one time point is done.')
        disp('Now start computing df/dvz, cez, and cez_PerpPar.')
    
        % ======================================================================
        % Compute vz derivatives of f at each x and t
        % df/dvz has been initialized
        for i = 1:nvperp
            for j = 1:ntheta
                for k = 2:nvz-1
                    dfdvz(i, j, k) = (f(ix, i, j, k+1, it) - f(ix, i, j, k-1, it)) / (vz_values(k+1) - vz_values(k-1));
                end
                % End points
                dfdvz(i, j, 1) = (f(ix, i, j, 2, it) - f(ix, i, j, 1, it)) / (vz_values(2) - vz_values(1));
                dfdvz(i, j, nvz) = (f(ix, i, j, nvz, it) - f(ix, i, j, nvz-1, it)) / (vz_values(nvz) - vz_values(nvz-1));
            end
        end    
        %======================================================================
        % Compute FPC at each x and t
        E0 = elecfield(t_final(it), [xval(ix), yval(ix), zval(ix)]);
        for i = 1:nvperp
            for j = 1:ntheta
                for k = 1:nvz
                    if (0==1) % Total v^2
	                    cex(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(1)*dfdvx(i,j,k);
	                    cey(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(2)*dfdvy(i,j,k);
	                    cez(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(3)*dfdvz(i,j,k);
                    else % Just component vy^2 or vx^2
                        % Note that here only the code below else is executed
	                    cez(ix, i, j, k, it) = - q * vz_values(k)^2 / 2. * E0(3) * dfdvz(i, j, k);
                    end
                end
            end
        end
    
        % Compute Alternate FPC Cprime and J dot E
        % but it does not compute anything because 0 neq 1
        if (0==1) %----------------------------------------
            for i = 1:nvx
                for j = 1:nvy
                    for k = 1:nvz
                          cpex(ix,i,j,k)=q*vx_values(i)*E0(1)*f(ix,i,j,k);
                          cpey(ix,i,j,k)=q*vy_values(j)*E0(2)*f(ix,i,j,k);
                          cpez(ix,i,j,k)=q*vz_values(k)*E0(3)*f(ix,i,j,k);
                    end
                end
            end
        end %----------------------------------------
    
        %======================================================================
    
        % Calculate cez_PerpPar
        % cez_PerpPar has been initialized
        cez_PerpPar = sum(cez, 3) * dtheta / (2. * pi);
    end %end for ix=1:nx %____________________________________________________
end % end for it = 1:nt

avg_cez_PerpPar = sum(cez_PerpPar(:, :, :, :, 1:size(cez_PerpPar, 5)-1), 5) * dt_final / (t_final(end) - t_final(1));
    
disp("Necessary quantities computation done.");
disp("Now start plotting routine.")
% ===============================================================
% iSH_FPC_Plotter % Plotting routine is called for every xval

% Create the meshgrid for plotting
[VZ, VPERP] = meshgrid(vz_values, vperp_values);

% Plot cez_PerpPar
tmpf(:,:) = squeeze(avg_cez_PerpPar);

% Set Figure Size
% scrsz = get(0,'ScreenSize');

if save_mat_data
    save(mat_data_name)
end

h10 = figure('Position', [1, 1, 1400, 700], 'Visible','on');
[~, h] = contourf(VZ, VPERP, tmpf, 50);
hold on;
xline(1.135, 'LineStyle','--', ...
    'LineWidth', 3);
xline(-1.135, 'LineStyle','--', 'LineWidth', 3);

set(h,'edgecolor','none');
set(gca,'FontSize',16, ...
'FontName','TimesNewRoman', ...
'FontWeight','normal', ...
'LineWidth',2)
colorbar('FontSize',16,'FontName','TimesNewRoman');

% xlim([vzmin - dv / 2, vzmax + dv / 2]);
% ylim([min(vperp_values) - dv / 2, max(vperp_values) + dv / 2]);

% Sets the colormap limits such that the center of the color bar is 0
cL = caxis;  
caxis([-max(abs(cL)) max(abs(cL))]); 

daspect([1 1 1]);
grid on;

title('$C_{E_z}(v_\parallel, v_\perp)$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

xlabel('$v_\parallel/v_{ti}$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

ylabel('$v_\perp/v_{ti}$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

colormap(bluewhitered);

if save_figure
    % Save the figure as a PNG file
    print(h10, fpcfigure_filename, '-dpng', '-r150');  % 150 dpi resolution

end


% Create variable sda for Shock-Drift Acceleration flag (to be set by hand)
% sda=zeros(nx,1);
% sda2=zeros(nx,1);
% sda3=zeros(nx,1);

% Compute sums over velocity-space for C_Ex, C_Ey, and C_Ez
% This should give us the energy density (Capital W in our notation)??
% I changed the following variable name from E_cex to Ws_cex
Ws_cez = sum(cez,[2 3 4]);














