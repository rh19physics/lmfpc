% SPM: Driver for Single Particle Motion in given E and B fields
% For iSH projects
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
% Numerical Method choices are made by setting variable (meth):
%    Possible selections are:
%    5 Adaptive Runge-Kutta 45 method with RelTol specified
% Input x0, v0, B0, E0,
% Compute x and v as a function of time

close all;
clear;
clc;

% Get the current date and time for filename suffix
current_time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
time_suffix = char(current_time);

test_SPM_Solver = false; % if true, make sure to choose -36 EB fields.
save_figure = false; % if true, save figures as png or pdf files
figure_name = sprintf("./plots_iSH/iSHSPM_%s", time_suffix);

save_mat_data = false; % if true, save all workspace 
mat_data_name = sprintf("./data_iSH/iSHSPM_%s.mat", time_suffix);

plot_circs = true; % Plot circles corresponding to particles' initial velocity in ICW frame
th_circs = linspace(0, pi, 100);

plot_circs_0 = true; % Plot circles centered at v = 0

keep_anno = false; % If true, label the parameters

% Get screen size
screenSize = get(0, 'ScreenSize'); % Get screen dimensions [left bottom width height]
screenWidth = screenSize(3);
screenHeight = screenSize(4);

figWidth = round(screenWidth / 3); % Set figure width (3 columns)
figHeight = round(screenHeight / 2); % Set figure height (half of screen height)

% Key parameters to define Plasma 
% [q,m,mpme,vac,bi,be,ma,theta,x0shock,lramp] = set_params;
[q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, delta_phi, kvalue] = set_params;

if q > 0
    charge_label = "q > 0, ";
elseif q < 0
    charge_label = "q < 0, ";
end

fprintf('field_choice = %1.1d', field_choice);
% em_eps
% waveT
% t_init
if field_choice == -653
    center_x1 = 2 * pi / waveT / kvalue(1);
    n_pos1_mode1 = (2. * pi / waveT - 1.) / kvalue(1);
    n_neg1_mode1 = (2. * pi / waveT + 1.) / kvalue(1);
    center_x2 = 2 * pi / waveT / kvalue(2);
    n_pos1_mode2 = (2. * pi / waveT - 1.) / kvalue(2);
    n_neg1_mode2 = (2. * pi / waveT + 1.) / kvalue(2);
else
    center_x1 = 2 * pi / waveT / kvalue;
    n_pos1_mode1 = (2. * pi / waveT - 1.) / kvalue;
    n_neg1_mode1 = (2. * pi / waveT + 1.) / kvalue;
end

% Setting timeSpan and initial conditions
tspan = [t_init + 0.* waveT, 2.*waveT];
% tspan = [2.*waveT, 8.*waveT];

% ============ Initial Conditions =============
% Define multiple initial conditions for (x, y, z, vx, vy, vz)
varying_which_ic = -1;
switch varying_which_ic
    case -1 % Testing case
        angle = 0. * pi / 4.;
        rad = 1;
        init_conditions = [
            % 0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), 1.135; 
        % 0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), 0.508 + 1 * 1.905; 
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1 - 1;    
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1;
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1 + 1;
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_neg1_mode1 - 1;    
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_neg1_mode1;
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_neg1_mode1 + 1;
        ];
    case 0 % Testing case
        angle = -1. * pi / 3.;
        rad = 2.;
        init_conditions = [
        % 0.1, 0.1, 0.1, 0., 2, n_pos1_mode2 - 1;    
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode2;
        % 0.1, 0.1, 0.1, 0., 2, n_pos1_mode2 + 1;
        % 0.1, 0.1, 0.1, 0., 2, n_pos1_mode1 - 1;    
        0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1;
        % 0.1, 0.1, 0.1, 0., 2, n_pos1_mode1 + 1;
        ];
    case 4
        init_conditions = [
            0.1, 0.1, 0.1, 4, 0, -1.4;
            0.1, 0.1, 0.1, 1.4, 0, -1.4;
            0.1, 0.1, 0.1, 1.2, 0, -1.4;    
            0.1, 0.1, 0.1, 1, 0, -1.4;
            0.1, 0.1, 0.1, 0.8, 0, -1.4;
            0.1, 0.1, 0.1, 0.6, 0, -1.4;
            0.1, 0.1, 0.1, 0.4, 0, -1.4;
            0.1, 0.1, 0.1, 0.2, 0, -1.4;
            % 0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -1.4;
        %    0.1, 0.1, 0.1, 0, 1, -1.4;
        ];
    case 5
        init_conditions = [
            0.1, 0.1, 0.1, 0, 4, -1.4;
            0.1, 0.1, 0.1, 0, 1.4, -1.4;
            0.1, 0.1, 0.1, 0, 1.2, -1.4;    
            0.1, 0.1, 0.1, 0, 1, -1.4;
            0.1, 0.1, 0.1, 0, 0.8, -1.4;
            0.1, 0.1, 0.1, 0, 0.6, -1.4;
            0.1, 0.1, 0.1, 0, 0.4, -1.4;
            0.1, 0.1, 0.1, 0, 0.2, -1.4;
        ];
    
    case 6
        init_conditions = [
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -2.8;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -2.6;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -2.4;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -2.2;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -2.0;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -1.8;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -1.6;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -1.4;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -1.2;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -1.0;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, -0.8;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 0.;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 0.2;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 0.4;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 1.0;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 2.2;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 2.413;
            0.1, 0.1, 0.1, sqrt(2)/2, sqrt(2)/2, 2.6;
        ];
end

% ntspan_cases = size(tspan, 1);
nIC_cases = size(init_conditions, 1);

% Preallocate storage for results using a structure
results = struct('t', [], 'y', [], 'vz', [], 'vperp', []);

lorhandle = @lorentz4;
options = odeset('RelTol', 1.0e-10);

% Define different colors for plots
colors = lines(nIC_cases);

for i = 1:nIC_cases
    % Extract the current initial condition
    y0 = init_conditions(i, 1:6);

    % Solve the ODE
    [t4, y4] = ode45(lorhandle, tspan, y0, options);

    % Extract velocity components
    vx = y4(:, 4);
    vy = y4(:, 5);
    vz = y4(:, 6);

    % Compute vperp
    vperp = sqrt(vx.^2 + vy.^2);
    
    % Store results in the structure
    results(i).t = t4;
    results(i).y = y4;
    results(i).vz = vz;
    results(i).vperp = vperp;

end

% Solve the ODE for each initial condition
h = figure('Position', [1 1 1300 750]);
hold on; % Allows multiple plots on the same figure

% Define different colors for plots
colors = lines(nIC_cases);

% Preallocate legend entries
legend_entries = cell(nIC_cases, 1);

% TIL xline and yline are always occupying the top layer...
xline(0, 'LineStyle','-', 'LineWidth', 3, 'Color', "#b2b2b2", Layer='bottom', HandleVisibility="off");

if field_choice == -653
    xline(n_neg1_mode1, 'LineStyle','--', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");
    xline(center_x1, 'LineStyle',':', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");
    xline(n_pos1_mode1, 'LineStyle','--', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");
    xline(n_neg1_mode2, 'LineStyle','--', 'LineWidth', 3, "Color", "b", Layer='bottom', HandleVisibility="off");
    xline(center_x2, 'LineStyle',':', 'LineWidth', 3, "Color", "b", Layer='bottom', HandleVisibility="off");
    xline(n_pos1_mode2, 'LineStyle','--', 'LineWidth', 3, "Color", "b", Layer='bottom', HandleVisibility="off");
else
    xline(n_neg1_mode1, 'LineStyle','--', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");
    xline(center_x1, 'LineStyle',':', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");
    xline(n_pos1_mode1, 'LineStyle','--', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");
end

yline(1, 'LineStyle',':', 'LineWidth', 3, Layer='bottom', HandleVisibility="off");

if plot_circs_0
    nth_cirs = length(th_circs);
    radii = zeros(4);
    j = 1;
for i = [1, 2, 3, 6]
    radii(j) = sqrt(init_conditions(i, 4)^2. + init_conditions(i, 5)^2. + init_conditions(i, 6)^2.);
    j = j + 1;
end

    nR = length(radii);
    x_circs_0 = zeros(nR, nth_cirs);
    y_circs_0 = zeros(nR, nth_cirs);
    for iR = 1:nR
        for ith_cirs = 1:nth_cirs
            x_circs_0(iR, ith_cirs) = radii(iR) .* cos(th_circs(ith_cirs));
            y_circs_0(iR, ith_cirs) = radii(iR) .* sin(th_circs(ith_cirs));
        end
    end
for iR = 1:nR
    p3 = plot(x_circs_0(iR, :), y_circs_0(iR, :), 'Color', 'k', 'LineWidth', 3, "HandleVisibility", "off");
    set(p3, 'Color', [p3.Color, 0.3]);
    hold on;
end
end

for i = 1:nIC_cases
    % Extract initial conditions
    vx0 = init_conditions(i, 4);
    vy0 = init_conditions(i, 5);
    vz0 = init_conditions(i, 6);
    
    % Compute initial vperp
    vperp0 = sqrt(vx0^2 + vy0^2);

    radius = sqrt((vz0 - center_x1).^2 + (vperp0 - 0).^2);
    % for ith_circs = 1:length(th_circs)
    x_circs = radius .* cos(th_circs) + center_x1;
    y_circs = radius .* sin(th_circs);
    % end

    if plot_circs
        if i <= 3
        % p2 = plot(x_circs, y_circs, 'Color', colors(i, :), 'LineWidth', 3);
        p2 = plot(x_circs, y_circs, 'Color', 'k', 'LineWidth', 3, HandleVisibility="off");
        set(p2, 'Color', [p2.Color, 0.7]);
        hold on;
        end
    end


    % Plot trajectory
    p1 = plot(results(i).vz, results(i).vperp, 'Color', colors(i, :), 'LineWidth', 4);
    set(p1, 'Color', [p1.Color, 0.4]);

    % Plot initial condition as a dot
    scatter(vz0, vperp0, 100, colors(i, :), 'filled', 'o', "HandleVisibility", "off");
    
    % Plot final position
    if tspan(end) > tspan(1) % physically, star happens after dot
        scatter(results(i).vz(end), results(i).vperp(end), 250, colors(i, :), 'filled', "pentagram", "HandleVisibility", "off");
    elseif tspan(end) < tspan(1) % physically, triangle happens before dot
        scatter(results(i).vz(end), results(i).vperp(end), 100, colors(i, :), 'filled', "^", "HandleVisibility", "off");
    end

    % Create legend entry with initial condition values
    legend_entries{i} = sprintf('x0=[%.1f,%.1f,%.1f], v0=[%.2f,%.2f,%.2f; %.2f]', ...
        init_conditions(i, 1), init_conditions(i, 2), init_conditions(i, 3), ...
        init_conditions(i, 4), init_conditions(i, 5), init_conditions(i, 6), sqrt(init_conditions(i, 4).^2+init_conditions(i, 5).^2.));
end
box on;

% Formatting plot
xlim([-3.7, 3.7]);
ylim([0, 3.7]);
format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$');
% Add legend with initial conditions
% legend(legend_entries, 'Interpreter', 'latex', 'Location', "northeast", "FontSize", 15);

anno_label = sprintf(charge_label+"$t_i = %1.1f T, t_{span}(1) = %1.1f T, t_{span}(end) = %1.1f T, field = %d, RSR = %3.2f, \\delta \\phi = %3.2f \\pi$", ...
    t_init/waveT, tspan(1)/waveT, tspan(end)/waveT, field_choice, em_eps, delta_phi/pi);

if keep_anno
title(anno_label, 'FontSize', 18, ...
     "Interpreter", "latex", 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
end
hold off;

% ========= Optional: Save Data =========
if save_mat_data
    save(mat_data_name);
end

if save_figure
    % print(h, figure_name, '-dpng', '-r150'); 
    % Match the figure size (convert pixels to inches, or set directly)
fig_width = 10;   % in inches
fig_height = 6;   % in inches
set(h, 'PaperSize', [fig_width fig_height]);
set(h, 'PaperPosition', [0 0 fig_width fig_height]);
    print(h, figure_name, '-dpdf', '-vector'); 
end

% ============== Plot 96 panels on one figure with varing vz0 ============
if (1==0)
    % Initial conditions
xic = 0.1;
yic = 0.1;
zic = 0.1;

% Re-scale I.C.s
% rescaleICs = false;
% 
% if field_choice == -45
%     if rescaleICs
%         vx0 = em_eps * 1.;
%         vy0 = em_eps * 1.;
%         vz0 = em_eps * 0.;
%     else
%         vx0 = 1.;
%         vy0 = 1.;
%         vz0 = 0.;
%     end
% else
%     vx0 = 0.;
%     vy0 = 0.;
%     vz0 = 1.135;
% end
vx0 = sqrt(2)/2;
vy0 = sqrt(2)/2;
% Define a vector of different vz0 values
vmag = 4;
nvz0 = 96;
vz0_vec = (-vmag:vmag/(nvz0/2):vmag-vmag/(nvz0/2));
% vz0_vec = (-vmag+1:(vmag-1)/(nvz0/2):(vmag-1)-(vmag-1)/(nvz0/2));

num_conditions = length(vz0_vec);

h = figure('Position',[1, 1, screenWidth*4, screenHeight*4], 'Visible','off');

t = tiledlayout(8, 12, 'TileSpacing','compact', 'Padding','compact');

for i = 1:num_conditions
    disp(i);
    ax = nexttile;
    
    % Set initial conditions with current vz0
    x0 = [xic yic zic ]';  % Initial position (Column vector)
    v0 = [vx0 vy0 vz0_vec(i)]';
    y0 = [x0' v0'];
    
    % Solve the ODE
    lorhandle = @lorentz4;
    options = odeset('RelTol', 1.0e-10);
    [t_temp, y_temp] = ode45(lorhandle, tspan, y0, options);
    
    % Extract and compute variables for plotting
    vx_temp = y_temp(:, 4);
    vy_temp = y_temp(:, 5);
    vz_temp = y_temp(:, 6);
    vperp_temp = sqrt(vx_temp.^2 + vy_temp.^2);
    vpar_temp = vz_temp;
    
    % Plot in the current tile
    plot(ax, vpar_temp, vperp_temp, 'LineWidth', 1);
    hold on;
    yline(1);
    xline(-1.4);
    xline(0.508);
    xline(2.413);
    % title(ax, sprintf('vz0 = %.2f', vz0_vec(i)), 'FontSize', 15);
    text(ax, -vmag*0.85, vmag*1.2, strcat("$v_{z0} = ", num2str(vz0_vec(i), "%4.3f"), "$"), ...
        "Interpreter", "latex", 'FontName','TimesNewRoman','FontSize',40,'FontWeight','bold');
    xlim(ax, [-4, 4]);
    ylim(ax, [0, 4*sqrt(2)]);
    format_subplot('$v_\parallel$', '$v_\perp$');
end

% sgtitle(sprintf('Velocity Space Trajectory with vx0 = %4.3f, vy0 = %4.3f', vx0, vy0), 'FontSize', 32, ...
%     'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
% 
anno_label = sprintf("$tspan(1) = %1.1f T, v_{x0} = %4.3f, v_{y0} = %4.3f, x_0 = %3.2f, y_0 = %3.2f, z_0 = %3.2f, tspan(end) = %1.1f T, field = %d, \delta \phi = %4.3f \pi$", ...
    tspan(1)/waveT, vx0, vy0, xic, yic, zic, tspan(end)/waveT, field_choice, delta_phi/pi);
sgtitle(anno_label, 'FontSize', 40, ...
     "Interpreter", "latex", 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
% % anno_label = "WTFqwjfdlbewqulhdweoq;hdiewqo;djziwq;hdiowq;ehdiewo;qfjqw;jfiewoq;jfwelq;fjewqjfoq'wjfewoqfj'qwejfopewjqfoqw;jfepdw'jqfoewjqfopweqj'";
% 
% annotation("textbox", [0., 0.9, 1, 0.3], 'String', ...
%     {anno_label},...
%     'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 40, ...
%     'HorizontalAlignment', 'left', ...
%     'FitBoxToText','off');

% anno_label = sprintf("At $t_1 = %1.1f T$,\n$v_{x0} = %4.3f, v_{y0} = %4.3f$, \n$x_0 = %3.2f, y_0 = %3.2f, z_0 = %3.2f$, \nand we solve SPM equations to $t_2 = %1.1f T$ with field = %d", ...
%     tspan(1)/waveT, vx0, vy0, xic, yic, zic, tspan(end)/waveT, field_choice);
% 
% annotation("textbox", [0.05, 0.8, 0.8, 0.15], 'String', ...
%     {anno_label}, 'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 45, ...
%     'HorizontalAlignment', 'left', 'FitBoxToText', 'on');


figure_filename = sprintf('./plots_iSH/iSHSPMVST_%s.png', time_suffix);
% exportgraphics(h, figure_filename, 'Resolution', 300);
print(h, figure_filename, '-dpng', '-r150'); 
end
% ================= Commented Old Ways of Plotting ============
if (1==1)
% Specify method for integration
meth = 5;  % RK45
    % Initial conditions
xic = 0.1;
yic = 0.1;
zic = -1.4;

% Re-scale I.C.s
% rescaleICs = false;
% 
% if field_choice == -45
%     if rescaleICs
%         vx0 = em_eps * 1.;
%         vy0 = em_eps * 1.;
%         vz0 = em_eps * 0.;
%     else
%         vx0 = 1.;
%         vy0 = 1.;
%         vz0 = 0.;
%     end
% else
%     vx0 = 0.;
%     vy0 = 0.;
%     vz0 = 1.135;
% end

vx0 = 0.;
vy0 = 1.;
vz0 = 1.4;


x0 = [xic yic zic ]';  % Initial position (Column vector)
v0 = [vx0 vy0 vz0 ]';  % Initial velocity (Column vector)


% Analytical solutions for ExB drift:
if test_SPM_Solver
    disp("You are testing the SPM Solver by comparing the numerical solution");
    disp("with the analytical solution for ExB drift.");
    disp("Did you choose the correct field setting?");
    disp("If so, press any key to continue. If not, stop the program.");
    pause
    nn = 10000;
    dt = (tspan(2) - tspan(1)) / nn;
    t_arr = (tspan(1):dt:tspan(2))';

    q = 1;
    m = 1;
    E_an = [0. 0.1 0.]';
    B_an = [0. 0. 1.]'; 
    
    Omega_0 = q * B_an(3) / m;
    v_EcrossB = cross(E_an, B_an)/B_an(3)^2;
    rl = sqrt((v0(1) - v_EcrossB(1))^2. + (v0(2) - v_EcrossB(2))^2.) / Omega_0;
    x_an = rl * sin(Omega_0 * t_arr) + E_an(2) * B_an(3) * t_arr / B_an(3)^2;
    y_an = rl * cos(Omega_0 * t_arr) + (x0(2) - rl);
    z_an = 0. * t_arr;
else
    % matlab accept an empty if body or else body as doing nothing
end
   


switch(meth)

case 5
    % Adaptive Runge-Kutta 45 method with RelTol specified
    method = 'Adaptive Runge-Kutta 45 method with RelTol specified'
    % Set function handle to Lorentz Force Law function
    lorhandle = @lorentz4;
    % Convert x0, v0 to y0(1:6) initial conditions
    y0 = [x0' v0'];
    % Calling ode45 for Runge-Kutta with option for setting error tolerance
    options = odeset('RelTol', 1.0e-10);
    tic;
    [t4,y4] = ode45(lorhandle, tspan, y0, options);
    nsteps = size(t4,1)
    ComputationTime = toc;
    fprintf("The computation time for tspan = %d is \n %f seconds.\n", tspan(end), ComputationTime);
end


% Ensure that t4 and y4 are defined
if ~exist('t4', 'var') || ~exist('y4', 'var')
    error('t4 and y4 must be defined. Run the differential equation solver first.');
end


% Extract positions and velocities from y4
x = y4(:, 1);
y = y4(:, 2);
z = y4(:, 3);
vx = y4(:, 4);
vy = y4(:, 5);
vz = y4(:, 6);

r2d = sqrt(x.^2 + y.^2);
vperp = sqrt(vx.^2 + vy.^2);

f = exp(- vx.^2. - vy.^2.- vx.^2.);

% h = figure();
% plot(t4, f);

if save_mat_data
    % Save all workspace variables
    save(mat_data_name)
end

    % ==================================
    % Quick look at data: plot data versus time
    % figure('Color', 'w');
    % figure;
    figureHandle_VSpace = figure;
    
    % Set the figure size (in pixels)
    figureHandle_VSpace.Position = [1, 150, figWidth, figHeight];  % [left, bottom, width, height]
    t = tiledlayout(1, 1);
    % =========================================
    ax1 = nexttile;
    plot(vz, vperp, 'LineWidth', 3);
    xlim([-4, 4]);
    ylim([0, 4*sqrt(2)]);
    format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$');
    
    % Add the title for the whole figure
    sgtitle('Velocity Space Trajectory', 'FontSize', 32, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
    
    % ==================================
    % Plot 3D trajectory
    % figure('Color', 'w');
    % figure;
    figureHandle_3DVSpacetraj = figure;
    
    % Set the figure size (in pixels)
    figureHandle_3DVSpacetraj.Position = [figWidth, 150, figWidth, figHeight];  % [left, bottom, width, height]
    
    
    plot3(vx, vy, vz, 'Color', 'b', 'LineWidth', 3, 'DisplayName','SPM Solver');
    hold on;
    
    if test_SPM_Solver
        plot3(x_an, y_an, z_an, 'Color', 'b', 'LineWidth', 3, 'LineStyle','--', 'DisplayName',"Analytical Sol");
        legend();
    else
        % do nothing
    end
    
    grid on;
    % Plot Labels
    % set(gca,'Color', 'w', 'FontSize',20,'FontName','TimesNewRoman','FontWeight','bold','LineWidth',2)
    set(gca,'FontSize',20,'FontName','TimesNewRoman','FontWeight','bold','LineWidth',2)
    xlabel('$v_x$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
    ylabel('$v_y$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
    zlabel('$v_z$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
    title('3D Velocity Space Trajectory of the Particle', 'FontName','TimesNewRoman');
    view([-20.883 19.155]);
    
    % Adjust paper size and position for saving
    set(figureHandle_3DVSpacetraj, 'PaperPositionMode', 'auto');
    set(figureHandle_3DVSpacetraj, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 8, 6]);
    
    % saveas(gcf, 'iSH_SPM_3D_Trajectory.png');  % Save as PNG file
    
    % ==================================
    % Quick look at data: plot data versus time
    % figure('Color', 'w');
    % figure;
    figureHandle_QL = figure;
    
    % Set the figure size (in pixels)
    figureHandle_QL.Position = [1, 2*960, 1400, 960];  % [left, bottom, width, height]
    
    subplot(2, 3, 1);
    plot(t4, x, 'LineWidth', 3);
    format_subplot('$t \Omega_i$', '$x / \rho_i$');
    
    subplot(2, 3, 2);
    plot(t4, y, 'LineWidth', 3);
    format_subplot('$t \Omega_i$', '$y / \rho_i$');
    
    subplot(2, 3, 3);
    plot(t4, z, 'LineWidth', 3);
    format_subplot('$t \Omega_i$', '$z / \rho_i$');
    
    subplot(2, 3, 4);
    plot(t4, vx, 'LineWidth', 3);
    format_subplot('$t \Omega_i$', '$v_x / v_{ti}$');
    
    subplot(2, 3, 5);
    plot(t4, vy, 'LineWidth', 3);
    format_subplot('$t \Omega_i$', '$v_y / v_{ti}$');
    
    subplot(2, 3, 6);
    plot(t4, vz, 'LineWidth', 3);
    format_subplot('$t \Omega_i$', '$v_z / v_{ti}$');
    
    % Add the title for the whole figure
    sgtitle('Overview of Position and Velocity Over Time', 'FontSize', 32, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
    
    % ==================================
    % Plot 3D trajectory
    % figure('Color', 'w');
    % figure;
    figureHandle_3Dtraj = figure;
    
    % Set the figure size (in pixels)
    figureHandle_3Dtraj.Position = [figWidth*2, 150, figWidth, figHeight];  % [left, bottom, width, height]
    
    
    plot3(x, y, z, 'Color', '#6aa84f', 'LineWidth', 3, 'DisplayName','SPM Solver');
    hold on;
    
    if test_SPM_Solver
        plot3(x_an, y_an, z_an, 'Color', 'b', 'LineWidth', 3, 'LineStyle','--', 'DisplayName',"Analytical Sol");
        legend();
    else
        % do nothing
    end
    
    grid on;
    % Plot Labels
    % set(gca,'Color', 'w', 'FontSize',20,'FontName','TimesNewRoman','FontWeight','bold','LineWidth',2)
    set(gca,'FontSize',20,'FontName','TimesNewRoman','FontWeight','bold','LineWidth',2)
    xlabel('$x$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
    ylabel('$y$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
    zlabel('$z$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
    title('3D Trajectory of the Particle', 'FontName','TimesNewRoman');
    view([-20.883 19.155]);
    
    % Adjust paper size and position for saving
    set(figureHandle_3Dtraj, 'PaperPositionMode', 'auto');
    set(figureHandle_3Dtraj, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 8, 6]);
    
    % saveas(gcf, 'iSH_SPM_3D_Trajectory.png');  % Save as PNG file
    
    
    % ==================================
    % Plot trajectories on 2D position and velocity space
    % figure('Color', 'w');
    % figure;
    figureHandle_2Dtraj = figure;
    
    % Set the figure size (in pixels)
    figureHandle_2Dtraj.Position = [1, 1, 1400, 960];  % [left, bottom, width, height]
    
    subplot(2, 3, 1);
    plot(x, z, 'LineWidth', 3);
    format_subplot('$x / \rho_i$', '$z / \rho_i$');
    
    subplot(2, 3, 4);
    plot(x, y, 'LineWidth', 3);
    format_subplot('$x / \rho_i$', '$y / \rho_i$');
    
    subplot(2, 3, 2);
    plot(vz, vx, 'LineWidth', 3);
    format_subplot('$v_z / v_{ti}$', '$v_x /v_{ti}$');
    
    subplot(2, 3, 5);
    plot(vz, vy, 'LineWidth', 3);
    format_subplot('$v_z /v_{ti}$', '$v_y / v_{ti}$');
    
    subplot(2, 3, 3);
    plot(r2d, z, 'LineWidth', 3);
    format_subplot('$r / \rho_i$', '$z / \rho_i$');
    
    subplot(2, 3, 6);
    plot(vz, vperp, 'LineWidth', 3);
    format_subplot('$v_z / v_{ti}$', '$v_\perp / v_{ti}$');
    
    % Add the title for the whole figure
    sgtitle('SPM Trajectories Projected on 2D Plane', 'FontSize', 36, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
    
    
    % Get the MATLAB version string
    matlab_version = version('-release');  % Returns '2023a', '2022b', etc.
    
    % Display the version
    fprintf('MATLAB Version: %s\n', matlab_version);
    
    % Check if the version is R2023a or earlier
    if strcmp(matlab_version, '2023a') || (str2double(matlab_version(1:4)) < 2023)
        disp('This is MATLAB R2023a or an earlier version.');
        % Continue
    else
        disp('This is MATLAB later than R2023a.');
        figureHandle_QL.Theme = 'Light';
        figureHandle_3Dtraj.Theme = 'Light';
        figureHandle_2Dtraj.Theme = 'Light';
    end
    
    figureHandle_QL.Visible = 'off';
    figureHandle_3Dtraj.Visible = "on";
    figureHandle_2Dtraj.Visible = "off";
    
    
    if save_figure
        % Save the figures as PNG files
        figure1_filename = sprintf('./plots_iSH/iSH_SPM_QL_fields%d_tspan%d_%s.png', field_choice, tspan(end), time_suffix);
        print(figureHandle_QL, figure1_filename, '-dpng', '-r150');  % 150 dpi resolution
        
        figure2_filename = sprintf('./plots_iSH/iSH_SPM_3Dtraj_fields%d_tspan%d_%s.png', field_choice, tspan(end), time_suffix);
        print(figureHandle_3Dtraj, figure2_filename, '-dpng', '-r150'); 
        
        figure3_filename = sprintf('./plots_iSH/iSH_SPM_2Dtraj_fields%d_tspan%d_%s.png', field_choice, tspan(end), time_suffix);
        print(figureHandle_2Dtraj, figure3_filename, '-dpng', '-r150');  
    end
end



% Function that formats subplots
% ================================================================================
function format_subplot(xlabel_text, ylabel_text)
    grid on;
    set(gca, 'FontSize', 25, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold', 'LineWidth', 3);
    daspect([1 1 1]);
    xlabel(xlabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 25, 'FontWeight', 'bold');
    ylabel(ylabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 25, 'FontWeight', 'bold');
end







