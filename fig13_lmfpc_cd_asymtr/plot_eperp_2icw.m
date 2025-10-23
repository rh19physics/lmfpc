% Plot Eperp versus z for 2-ICW field -653

close all
clear
clc

format long

[~, ~, ~, ~, ~, ~, ~, field_choice, em_eps, waveT, t_init, delta_phi, kvalue] = set_params;


save_fig = true;

% Set position and time values for plotting
x_value = 0.1; 
y_value = 0.1; 
z_values = linspace(-(2.*pi/kvalue(2)), 2.*pi/kvalue(2), 10000);
t_value = t_init + 2.*waveT;

% Preallocate arrays for E values
E_values = zeros(length(z_values), 3);

% Calculate E for each t_value
for i = 1:length(z_values)
    E_values(i, :) = elecfield(t_value, [x_value, y_value, z_values(i)]);
end

% Extract E components for plotting
Ex = E_values(:, 1);
Ey = E_values(:, 2);
Ez = E_values(:, 3);

total_fig = figure('Position', [100, 100, 1000, 200], 'Theme', 'Light');
plot(z_values, sqrt(Ex.^2. + Ey.^2.), 'Color', 'k', 'LineStyle', '-', 'DisplayName', 'kpar > 0', 'LineWidth', 3);
hold on;
% xline([-2.9, -1.4, 0.1, 1.6, 3.1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
xline([-2.9, -1.4, 0.1, 1.6], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
format_subplot('$z$', '$\frac{E_\perp}{v_{ti} B_0 / c}$');
xlim([-(2.*pi/kvalue(2)), 2.*pi/kvalue(2)]);
text(-10, 0.039, "$(q)$", "Interpreter","latex", 'VerticalAlignment','top', 'FontSize',20);

if save_fig
    pngname = './png/Eperp2ICW.png';
    epsname = './eps/Eperp2ICW.eps';
    pdfname = './pdf/Eperp2ICW.pdf';

    % Save png and eps
    print(total_fig, pngname, '-dpng', '-r150');  % 150 dpi resolution
    print(total_fig, epsname, '-painters','-depsc','-r150');

    % Save pdf
    % Match the figure size (convert pixels to inches, or set directly)
    fig_width = 10;   % in inches
    fig_height = 2;   % in inches
    set(total_fig, 'PaperSize', [fig_width fig_height], 'PaperPosition', [0 0 fig_width fig_height]);
    print(total_fig, pdfname, '-dpdf', '-vector');
    disp("figure saved");
end

% ================ functions ================
% formatting plots
function format_subplot(xlabel_text, ylabel_text)
    % legend;
    grid on;
    set(gca, 'Fontsize', 20, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold', 'LineWidth', 2);
    xlabel(xlabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel(ylabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 24, 'FontWeight', 'bold');
    % title(title_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 24, 'FontWeight', 'bold');
end