close all
clear
clc

format long

save_plot = true;
save_vecplot = true;
keep_anno = true;

figure_number = 3;

switch(figure_number)
    case 1
        figure_label = ["(a)", "(b)"];
        load("../lmfpc_matlab/data/iSHLDV12_20250317_185835.mat");
    case 2
        figure_label = ["(c)", "(d)"];
        load("../lmfpc_matlab/data/iSHLDV12_20250417_103327.mat");
    case 3
        figure_label = ["(e)", "(f)"];
        load("../lmfpc_matlab/data/iSHLDV12_20250417_104902.mat");
end

% Re-calculate time average
avg_cez_PerpPar = sum(cez_PerpPar(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
tmpf = squeeze(avg_cez_PerpPar);

h10 = figure('Position', [100, 200, 1000, 800], "Visible", "off");

% Define manual positions: [left, bottom, width, height]
top_pos = [0.13, 0.43, 0.72, 0.48];   % top panel (larger height)
bot_pos = [0.13, 0.17, 0.72, 0.2];   % bottom panel

% Top axes
ax1 = axes('Position', top_pos);
[~, h] = contourf(VZ, VPERP, tmpf, 50);
hold on;
xline(1.135, 'LineStyle','--', 'LineWidth', 3);
xline(-1.135, 'LineStyle','--', 'LineWidth', 3);

xlim([min(vz_values), max(vz_values)]);

set(h,'edgecolor','none');
set(ax1, 'FontSize',24, 'FontName','TimesNewRoman', 'LineWidth',2, ...
    'XTickLabel', []);  % hide x-tick labels on top panel

text(-3, 3*0.97, figure_label(1), ...
    'Interpreter','latex', 'VerticalAlignment','top', FontSize=32);

% Annotations
first_line = sprintf("Field Choice: %d, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", field_choice, xval, yval, zval);
if size(kpardi, 2) == 1
    second_line = sprintf("Settings: $RSR = %3.2f, k_\\parallel \\rho_i = %4.3f, k_\\perp \\rho_i = 1, \\delta \\phi = %3.2f \\pi$", ...
        em_eps, kpardi, delta_phi/pi);
else
    second_line = sprintf("Settings: $RSR = %3.2f, k_{\\parallel, 1} \\rho_i = %4.3f, k_{\\parallel, 2} \\rho_i = %4.3f, k_\\perp \\rho_i = 1, \\delta \\phi = %3.2f \\pi$", ...
                    em_eps, kpardi(1), kpardi(2), delta_phi/pi);
end

third_line = sprintf('$t_i = %.0f T, t_f = (%.0f, %4.3f T; %4.3f T), (n_{v_\\perp}, n_{\\theta}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$', ...
    t_init/waveT, t_final(1), t_final(end-1)/waveT, dt_final/waveT, nvperp, ntheta, nvz);

if keep_anno
text(-3, 2, {[first_line], [second_line], [third_line]}, ...
    'Interpreter','latex', 'VerticalAlignment','top', FontSize=13);
end

title("$C_{E_z}(v_\parallel, v_\perp)$", ...
    'Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize', 32, ...
    'FontWeight','bold');
      
ylabel('$v_\perp/v_{ti}$','Interpreter','latex', ...
    'FontName','TimesNewRoman', ...
    'FontSize',32, ...
    'FontWeight','bold'); 

cb = colorbar('FontSize',24,'FontName','TimesNewRoman');
cb.Position(1) = 0.85;  % reposition colorbar to match
cb.Position(2) = top_pos(2);
cb.Position(4) = top_pos(4);  % same height as panel
% Sets the colormap limits such that the center of the color bar is 0
cL = caxis;  
caxis([-max(abs(cL)) max(abs(cL))]); 

daspect([1 1 1]);
grid on;

colormap(bluewhitered);

        % Bottom axes
        ax2 = axes('Position', bot_pos);
        tmpf_1DCez = sum(tmpf, 1) * dv_perp;
        
        fig = plot(vz_values, tmpf_1DCez, 'LineWidth', 3);
        set(gca,'FontSize',24, ...
        'FontName','TimesNewRoman', ...
        'FontWeight','normal', ...
        'LineWidth',2)

        hold on;

        text(-3, max(tmpf_1DCez)*0.97, figure_label(2), ...
            'Interpreter','latex', 'VerticalAlignment','top', FontSize=32);

        % plot(1.1131, 0., '.', 'MarkerSize', 40, 'Color','Black');
        % plot(1.8850, 0., '.', 'MarkerSize', 40, 'Color','Black');

        xline(1.135, 'LineStyle','--', ...
            'LineWidth', 3);
        xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        yline(0, 'Color','k', 'LineWidth',2);

        grid on;
        xlim([min(vz_values), max(vz_values)]);
        ylim([-max(abs(tmpf_1DCez))*1.1, max(abs(tmpf_1DCez))*1.1]);
        
        % title({[first_line], ...
        %     [second_line]}, ...
        %     'Interpreter','latex', ...
        %     'FontName','TimesNewRoman', ...
        %     'FontSize', 14, ...
        %     'FontWeight','bold');

        xlabel('$v_\parallel/v_{ti}$','Interpreter','latex', ...
            'FontName', 'TimesNewRoman', ...
            'FontSize', 32, ...
            'FontWeight', 'bold');

        ylabel('$C_{E_z}(v_\parallel)$','Interpreter','latex', ...
            'FontName', 'TimesNewRoman', ...
            'FontSize', 32, ...
            'FontWeight', 'bold');

        % ylabel('$C_{E_z}(v_\parallel) = \int C_{E_z}(v_\parallel, v_\perp) d v_\perp$','Interpreter','latex', ...
        %     'FontName', 'TimesNewRoman', ...
        %     'FontSize', 24, ...
        %     'FontWeight', 'bold');

        % Link x-axes
linkaxes([ax1, ax2], 'x');
        
        if save_plot
            if save_vecplot
                % Set paper units to inches (or points, cm, etc.)
set(h10, 'PaperUnits', 'inches');

% Match the figure size (convert pixels to inches, or set directly)
fig_width = 8;   % in inches
fig_height = 6.4;   % in inches
set(h10, 'PaperSize', [fig_width fig_height]);
set(h10, 'PaperPosition', [0 0 fig_width fig_height]);
figure_name = sprintf("iSHLDV12_%s", time_suffix);
print(h10, figure_name, '-dpdf', '-vector');
% exportgraphics(h10, figure_name, "ContentType","vector");
print(h10, figure_name, '-painters','-depsc','-r150');
            else
            figure_name = sprintf("iSHLDV12_%s", time_suffix);
            print(h10, figure_name, '-dpng', '-r150');  % 150 dpi resolution
            end
        end