% Determine Distribution Using Liouville-Mapping: ion Stochastic Heating
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

% format long

disp('=======================================================================');
disp("You've started iSH_FPC_Calculator_V2_1_SingleICWCD_VxVy.");
% This version only calculates f, df/dz, cez
% This version starts with vperp and theta grids.

% Get the MATLAB version string
matlab_version = version('-release'); % Returns '2023a', '2022b', etc.

% Display the version
% fprintf('MATLAB Version: %s\n', matlab_version);

% Key parameters to define Plasma and turbulence
[q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, filename, delta_phi, kvalue] = set_params;

% Setting Position to Compute f(vx,vy)----------------------
% Set x limits
% xval = (0.1:-1.0:0.1)';
xval = 0.1;
yval = 0.1;
zval = 0.1;
nx = length(xval);

% Create (vx, vy, vz) grid ================================================
% Create Velocity Space Grid (v_ti units)

vlength = 4.;
% Compute max and min ends of velocity space axes
dvx = 2. * vlength / 39.;
dvy = 2. * vlength / 39.;
dvz = 2. * vlength / 39.;

% Set up (vx, vy, vz) grid
vx_values = (-vlength:dvx:vlength)
vy_values = (-vlength:dvy:vlength)
vz_values = (-vlength:dvz:vlength)
% vz_values = -1.7734;

nvx = length(vx_values);
nvy = length(vy_values);
nvz = length(vz_values);

% Set t_final slices
% ti = -333.; % now t_init has been specified in set_params
dt_final = waveT / 40.;
t_final = (0.:dt_final:2.*waveT)' % The correlation interval
nt_final = length(t_final);

% END Create (vperp, theta, vz) grid ================================================

% Initialize f, dfdvz, cez, cez_PerpPar
f = zeros(nx, nvx, nvy, nvz, nt_final);

dfdvx = zeros(nx, nvx, nvy, nvz, nt_final);
dfdvy = zeros(nx, nvx, nvy, nvz, nt_final);
dfdvz = zeros(nx, nvx, nvy, nvz, nt_final);

cex = zeros(nx, nvx, nvy, nvz, nt_final);
cey = zeros(nx, nvx, nvy, nvz, nt_final);
cez = zeros(nx, nvx, nvy, nvz, nt_final);

cex_VxVy = zeros(nx, nvx, nvy, nt_final);
cey_VxVy = zeros(nx, nvx, nvy, nt_final);
cez_VxVy = zeros(nx, nvx, nvy, nt_final);

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

save_figure = false; % if true, save figures as png files
fpcfigure_filename = sprintf('./plots_iSH/iSHCDV21_%s.png', time_suffix);

save_mat_data = false; % if true, save all workspace 
mat_data_name = sprintf("./data_iSH/iSHCDV21_%s.mat", time_suffix);

% pause
disp("Initiation done.");
disp("Now start the numerical calculation process for each velocity space grid point.")

% LOOP OVER x-positions for distribution function calculation
% This is the core code of this script.
% It calculates f(xval, vx, vy, vz), df/dv, FPC, j dot E, reduced FPC
for it_final = 1:nt_final
    fprintf('t_final= %3.2f\n', t_final(it_final));
    for ix = 1:nx %________________________________________________________________
        % Loop over each point in the (vx, vy, vz) grid
        for ivz = 1:nvz
            tic; % Start timing the time loop
            for ivy = 1:nvy
                for ivx = 1:nvx
                    vx0 = vx_values(ivx);
                    vy0 = vy_values(ivy);
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
                    tspan = [t_final(it_final), t_init]';
                    [t4,y4] = ode45(lorhandle, tspan, y0, options);
                    nsteps = size(t4,1);
        
                    % y4(end) gives x and v at the "initial" time (because of integrating backwards)
                    % here, acrodding to Liouville's theorem, the f value calculated from 
                    % the "initial" time, which is a Maxwellian, has been assigned to 
                    % the phase-space grid we want, at (x, vx, vy, vz)
                    % need to change the expression of f for iSH
    
                    f(ix, ivx, ivy, ivz, it_final) = exp(- y4(end,4)^2. - y4(end,5)^2.- y4(end,6)^2.);
        
                end % vx loop 
            end  % vy loop

            fprintf('vx = %3.2f vy = %3.2f vz =  %3.2f \n', vx_values(ivx), vy_values(ivy),  vz_values(ivz));
            % Stop timing and output the computation time
            elapsedTime = toc;
            fprintf('Time taken for computing f at this vz loop is\n %f seconds\n', elapsedTime);

            % if (mod(ivperp,4)==0)
            %     [vx_values(ivperp)] % This output serves as a progress bar.
            % end
    
        end % vz loop
    end % x loop
end % t_final loop

disp('Numerical calculation at one time point is done.')
disp('Now start computing df/dvx, df/dvy, df/dvz,')
disp('cex, cey, cez, cex_VxVy, cey_VxVy, cez_VxVy.')

for it_final = 1:nt_final
    for ix = 1:nx
        E0 = elecfield(t_final(it_final), [xval(ix), yval(ix), zval(ix)]);    
        % ======================================================================
        % Compute vx derivatives of f at each x and t
        % df/dvx has been initialized
        for k = 1:nvz
            for j = 1:nvy
                for i = 2:nvx-1
                    dfdvx(ix, i, j, k, it_final) = (f(ix, i+1, j, k, it_final) - f(ix, i-1, j, k, it_final)) / (vx_values(i+1) - vx_values(i-1));
                end
                % End points
                dfdvx(ix, 1, j, k, it_final) = (f(ix, 2, j, k, it_final) - f(ix, 1, j, k, it_final)) / (vx_values(2) - vx_values(1));
                dfdvx(ix, nvx, j, k, it_final) = (f(ix, nvx, j, k, it_final) - f(ix, nvx-1, j, k, it_final)) / (vx_values(nvx) - vx_values(nvx-1));
            end
        end    
        %======================================================================
        % Compute vy derivatives of f at each x and t
        % df/dvy has been initialized
        for i = 1:nvx
            for k = 1:nvz
                for j = 2:nvy-1
                    dfdvy(ix, i, j, k, it_final) = (f(ix, i, j+1, k, it_final) - f(ix, i, j-1, k, it_final)) / (vy_values(j+1) - vy_values(j-1));
                end
                % End points
                dfdvy(ix, i, 1, k, it_final) = (f(ix, i, 2, k, it_final) - f(ix, i, 1, k, it_final)) / (vy_values(2) - vy_values(1));
                dfdvy(ix, i, nvy, k, it_final) = (f(ix, i, nvy, k, it_final) - f(ix, i, nvy-1, k, it_final)) / (vy_values(nvy) - vy_values(nvy-1));
            end
        end    
        %======================================================================
        % Compute vz derivatives of f at each x and t
        % df/dvz has been initialized
        if nvz < 3
            % continue
        else
            for i = 1:nvx
                for j = 1:nvy
                    for k = 2:nvz-1
                        dfdvz(ix, i, j, k, it_final) = (f(ix, i, j, k+1, it_final) - f(ix, i, j, k-1, it_final)) / (vz_values(k+1) - vz_values(k-1));
                    end
                    % End points
                    dfdvz(ix, i, j, 1, it_final) = (f(ix, i, j, 2, it_final) - f(ix, i, j, 1, it_final)) / (vz_values(2) - vz_values(1));
                    dfdvz(ix, i, j, nvz, it_final) = (f(ix, i, j, nvz, it_final) - f(ix, i, j, nvz-1, it_final)) / (vz_values(nvz) - vz_values(nvz-1));
                end
            end
        end
        %======================================================================
        
        %======================================================================
        % Compute Cex at each x and t     
        for i = 1:nvx
            for j = 1:nvy
                for k = 1:nvz
                    if (0==1) % Total v^2
	                    cex(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(1)*dfdvx(i,j,k);
	                    cey(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(2)*dfdvy(i,j,k);
	                    cez(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(3)*dfdvz(i,j,k);
                    else % Just component vy^2 or vx^2
                        % Note that here only the code below else is executed
	                    cex(ix, i, j, k, it_final) = - q * vx_values(i)^2 / 2. * E0(1) * dfdvx(ix, i, j, k, it_final);
                    end
                end
            end
        end

        %======================================================================
        % Compute Cey at each x and t
        for i = 1:nvx
            for j = 1:nvy
                for k = 1:nvz
                    if (0==1) % Total v^2
	                    cex(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(1)*dfdvx(i,j,k);
	                    cey(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(2)*dfdvy(i,j,k);
	                    cez(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(3)*dfdvz(i,j,k);
                    else % Just component vy^2 or vx^2
                        % Note that here only the code below else is executed
	                    cey(ix, i, j, k, it_final) = - q * vy_values(j)^2 / 2. * E0(2) * dfdvy(ix, i, j, k, it_final);
                    end
                end
            end
        end
        %======================================================================
        % Compute Cez at each x and t
        if nvz < 3
            % continue
        else
            for i = 1:nvx
                for j = 1:nvy
                    for k = 1:nvz
                        if (0==1) % Total v^2
                            cex(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(1)*dfdvx(i,j,k);
                            cey(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(2)*dfdvy(i,j,k);
                            cez(ix,i,j,k)=-q*(vx_values(i)^2+vy_values(j)^2+vz_values(k)^2)/2.*E0(3)*dfdvz(i,j,k);
                        else % Just component vy^2 or vx^2
                            % Note that here only the code below else is executed
                            cez(ix, i, j, k, it_final) = - q * vz_values(k)^2 / 2. * E0(3) * dfdvz(ix, i, j, k, it_final);
                        end
                    end
                end
            end
        end
        %======================================================================
    
        % Calculate cex_VxVy, cey_VxVy, cez_VxVy
        % cez_VxVy has been initialized
        cex_VxVy = sum(cex, 4) * dvz;
        cey_VxVy = sum(cey, 4) * dvz;
        if nvz < 3
            % continue
        else
            cez_VxVy = sum(cez, 4) * dvz;
        end
    end %end for ix=1:nx %____________________________________________________
end % end for it = 1:nt

tavg_cex_VxVy = sum(cex_VxVy(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
tavg_cey_VxVy = sum(cey_VxVy(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
if nvz < 3
    % continue
else
    tavg_cez_VxVy = sum(cez_VxVy(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
end
    
disp("Necessary quantities computation done.");
disp("Now start plotting routine.")
% ===============================================================
% iSH_FPC_Plotter % Plotting routine is called for every xval

% Create the meshgrid for plotting
[VX, VY] = meshgrid(vx_values, vy_values);

% Set Figure Size
% scrsz = get(0,'ScreenSize');

if save_mat_data
    save(mat_data_name)
end

h10 = figure('Position', [1, 1, 1400, 700], 'Visible','on');
% Check if the version is R2023a or earlier
if strcmp(matlab_version, '2023a') || (str2double(matlab_version(1:4)) < 2023)
	% disp('This is MATLAB R2023a or an earlier version.');
	% Continue
else
	% disp('This is MATLAB later than R2023a.');
	h10.Theme = 'Light';
end

subplot(1, 2, 1);
% Plot tavg_cex_VxVy
tmpf(:,:) = squeeze(tavg_cex_VxVy);
[~, h1] = contourf(VX, VY, transpose(tmpf), 50);
hold on;
% xline(1.135, 'LineStyle','--', 'LineWidth', 3);
% xline(-1.135, 'LineStyle','--', 'LineWidth', 3);

set(h1,'edgecolor','none');
set(gca,'FontSize',16, ...
'FontName','TimesNewRoman', ...
'FontWeight','normal', ...
'LineWidth',2)
colorbar('FontSize',16,'FontName','TimesNewRoman');

% xlim([vzmin - dv / 2, vzmax + dv / 2]);
% ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);

% Sets the colormap limits such that the center of the color bar is 0
cL = caxis;  
caxis([-max(abs(cL)) max(abs(cL))]); 

daspect([1 1 1]);
grid on;

title('$C_{E_x}(v_x, v_y)$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

xlabel('$v_x/v_{ti}$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

ylabel('$v_y/v_{ti}$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

% =========================================
subplot(1, 2, 2);
% Plot tavg_cey_VxVy
tmpf(:,:) = squeeze(tavg_cey_VxVy);
[~, h2] = contourf(VX, VY, transpose(tmpf), 50);
hold on;
% xline(1.135, 'LineStyle','--', 'LineWidth', 3);
% xline(-1.135, 'LineStyle','--', 'LineWidth', 3);

set(h2,'edgecolor','none');
set(gca,'FontSize',16, ...
'FontName','TimesNewRoman', ...
'FontWeight','normal', ...
'LineWidth',2)
colorbar('FontSize',16,'FontName','TimesNewRoman');

% xlim([vzmin - dv / 2, vzmax + dv / 2]);
% ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);

% Sets the colormap limits such that the center of the color bar is 0
cL = caxis;  
caxis([-max(abs(cL)) max(abs(cL))]); 

daspect([1 1 1]);
grid on;

title('$C_{E_y}(v_x, v_y)$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

xlabel('$v_x/v_{ti}$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',24, ...
    'FontWeight','bold');

ylabel('$v_y/v_{ti}$','Interpreter','latex', ...
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














