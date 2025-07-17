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

save_figure = true; % if true, save figures as png or pdf files

plot_circs = true; % Plot circles corresponding to particles' initial velocity in ICW frame
plot_circs_0 = true; % Plot circles centered at v = 0, i.e. lab frame
th_circs = linspace(0, pi, 100);

keep_anno = true; % If true, label the parameters

% Cross-checking the wanted parameters and parameters from set_params.m
params_wanted = struct( ...
    'q', 1., ...
    'm', 1., ...
    'field_choice', -63, ...
    'em_eps', 0.02, ...
    't_init_multiple', -4, ... 
    'delta_phi', 0. ...
);

[q, m, mime, tite, bi, be, vtic, field_choice, em_eps, waveT, t_init, delta_phi, kvalue] = set_params;
params_set = struct( ...
    'q', q, ...
    'm', m, ...
    'field_choice', field_choice, ...
    'em_eps', em_eps, ...
    't_init_multiple', t_init/waveT, ... 
    'delta_phi', delta_phi ...
);

disp("Start cross checking parameters...");
if isequal(params_wanted, params_set)
    disp("Correct local field.");
else
    disp("Wrong local field. Please press Ctrl+C.");
    disp("Differences:");
    show_params_differences(params_wanted, params_set);
    pause;
end

if q > 0
    charge_label = "q > 0, ";
elseif q < 0
    charge_label = "q < 0, ";
end

if size(kvalue, 2) == 2
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
angle = 0. * pi / 4.;
rad = 1;
init_conditions = [
0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1 - 1;    
0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1;
0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_pos1_mode1 + 1;
0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_neg1_mode1 - 1;    
0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_neg1_mode1;
0.1, 0.1, 0.1, rad * cos(angle), rad * sin(angle), n_neg1_mode1 + 1;
];

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

if size(kvalue, 2) == 2
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

anno_label = sprintf(charge_label+"$t_i = %1.1f T, t_{span}(1) = %1.1f T, t_{span}(end) = %1.1f T, field = %d, RSR = %3.2f, \\delta \\phi = %3.2f \\pi$", ...
    t_init/waveT, tspan(1)/waveT, tspan(end)/waveT, field_choice, em_eps, delta_phi/pi);

if keep_anno
% Add legend with initial conditions
    legend(legend_entries, 'Interpreter', 'latex', 'Location', "northeast", "FontSize", 15);
    title(anno_label, 'FontSize', 18, "Interpreter", "latex", 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
end

if save_figure
    if keep_anno
        pngname = "./png/spm_icw_gyro_anno.png";
        epsname = "./eps/spm_icw_gyro_anno.eps";
        pdfname = "./pdf/spm_icw_gyro_anno.pdf";       
    else
        pngname = "./png/spm_icw_gyro.png";
        epsname = "./eps/spm_icw_gyro.eps";
        pdfname = "./pdf/spm_icw_gyro.pdf";
    end
    % Save png and eps
    print(h, pngname, '-dpng', '-r150');  % 150 dpi resolution
    print(h, epsname, '-painters','-depsc','-r150');
    
    % Save pdf
    % Match the figure size (convert pixels to inches, or set directly)
    fig_width = 10;   % in inches
    fig_height = 6;   % in inches
    set(h, 'PaperSize', [fig_width fig_height]);
    set(h, 'PaperPosition', [0 0 fig_width fig_height]);
    print(h, pdfname, '-dpdf', '-vector');
    disp('Figure saved.');
end

disp('Script finished successfully.');

% Function that formats subplots
% ================================================================================
function format_subplot(xlabel_text, ylabel_text)
    grid on;
    set(gca, 'FontSize', 25, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold', 'LineWidth', 3);
    daspect([1 1 1]);
    xlabel(xlabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 25, 'FontWeight', 'bold');
    ylabel(ylabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 25, 'FontWeight', 'bold');
end

% Display differences between fieldnames of s1 and s2
function show_params_differences(s1, s2)
    params = fieldnames(s1);
    for i = 1:numel(params)
        param_name = params{i};
        param1 = s1.(param_name);
        param2 = s2.(param_name);
        if ~isequal(param1, param2)
            fprintf("Parameter '%s' differs: \n", param_name);
            disp([' Wanted parameter: ', mat2str(param1)]);
            disp([' Current set parameter: ', mat2str(param2)]);
        end
    end
end





