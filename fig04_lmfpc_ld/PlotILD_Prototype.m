close all
clear
clc

format long

% Get the MATLAB version string
% matlab_version = version('-release');  % Returns '2023a', '2022b', etc.

% Display the version
% fprintf('MATLAB Version: %s\n', matlab_version);

data_src = 12; % 12 means data from V1.2, 13 means data from V1.3
% load("./data_iSH/iSHLDV12_20250427_171656.mat");
datapath = "./data_iSH/iSHLDV12_20250317_185835.mat";
save_plot = true;
save_vecplot = true;
keep_anno = false;

output_reversal_time = false;
input_vz = 1.25;
num_to_label = 5;

switch(data_src)

    case 12
        % load('~/Desktop/iSHLDV12_20241126_173804.mat');
        load(datapath);
        if exist("delta_phi") == 0
            delta_phi = 0.;
        end

        [~, ~, ~, ~, ~, ~, ~, field_choice_local, ~, ~, ~, ~, delta_phi_local, kpardi] = set_params;

        if field_choice == -495 || field_choice == -4951
            if field_choice_local == field_choice && delta_phi == delta_phi_local
                disp("Correct local field.");
            else
                disp("Wrong local field. Please press Ctrl+C.");
                pause;
            end
        else
            if field_choice_local == field_choice
                disp("Correct local field.");
            else
                disp("Wrong local field. Please press Ctrl+C.");
                pause;
            end
        end
        if exist("yval") 
        else
            yval = xval;
            zval = xval;
        end

        avg_cez_PerpPar = sum(cez_PerpPar(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
        tmpf = squeeze(avg_cez_PerpPar);

        h10 = figure('Position', [100, 200, 1000, 800], "Theme", "Light", "Visible", "off");
        % Check if the version is R2023a or earlier
        % if strcmp(matlab_version, '2023a') || (str2double(matlab_version(1:4)) < 2023)
        %     disp('This is MATLAB R2023a or an earlier version.');
        %     % Continue
        % else
        %     disp('This is MATLAB later than R2023a.');
        %     h10.Theme = 'Light';
        % end

        % Define manual positions: [left, bottom, width, height]
        top_pos = [0.13, 0.43, 0.72, 0.48];   % top panel (larger height)
        bot_pos = [0.13, 0.17, 0.72, 0.2];   % bottom panel

        % Top axes
        ax1 = axes('Position', top_pos);
        [~, h] = contourf(VZ, VPERP, tmpf, 50);
        hold on;
        xline(1.135, 'LineStyle','--', ...
            'LineWidth', 3);
        xline(-1.135, 'LineStyle','--', 'LineWidth', 3);

        xlim([min(vz_values), max(vz_values)]);
        
        set(h,'edgecolor','none');
        set(ax1, 'FontSize',24, 'FontName','TimesNewRoman', 'LineWidth',2, ...
    'XTickLabel', []);  % hide x-tick labels on top panel
        % set(gca, 'Color', 'white', 'Gridcolor', 'black');
        % set(gcf, 'Color', 'white');
        % colorbar('FontSize',16,'FontName','TimesNewRoman');

        text(-3, 3*0.97, "(e)", ...
            'Interpreter','latex', 'VerticalAlignment','top', FontSize=32);

        switch (field_choice)
            case -45
                first_line = "$C_{E_z}(v_\parallel, v_\perp)$ with 1-KAW and $k_\parallel > 0$";
                second_line = sprintf("No Window Function, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", xval, yval, zval);
            case -451
                first_line = "$C_{E_z}(v_\parallel, v_\perp)$ with 1-KAW and $k_\parallel > 0$";
                second_line = sprintf("Window Function, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", xval, yval, zval);
            case -49
                first_line = "$C_{E_z}(v_\parallel, v_\perp)$ with 1-KAW and $k_\parallel < 0$";
                second_line = sprintf("No Window Function, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", xval, yval, zval);
            case -491
                first_line = "$C_{E_z}(v_\parallel, v_\perp)$ with 1-KAW and $k_\parallel < 0$";
                second_line = sprintf("Window Function, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", xval, yval, zval);
            case -495
                first_line = "$C_{E_z}(v_\parallel, v_\perp)$ with 2-KAW";
                second_line = sprintf("No Window Function, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", xval, yval, zval);
            case -4951
                first_line = "$C_{E_z}(v_\parallel, v_\perp)$ with 2-KAW";
                second_line = sprintf("Window Function, $(x, y, z) = (%3.2f, %3.2f, %3.2f)$", xval, yval, zval);
        end

        switch (field_choice)
            case {-491, -49, -451, -45}
                third_line = sprintf("Settings: $RSR = %3.2f, k_\\parallel \\rho_i = %4.3f, k_\\perp \\rho_i = 1, \\delta \\phi = %3.2f \\pi$", ...
                    em_eps, kpardi, delta_phi/pi);
            case {-495, -4951} % 2-KAWs, -495: no window; -4951: with window
                third_line = sprintf("Settings: $RSR = %3.2f, k_{\\parallel, 1} \\rho_i = %4.3f, k_{\\parallel, 2} \\rho_i = %4.3f, k_\\perp \\rho_i = 1, \\delta \\phi = %3.2f \\pi$", ...
                    em_eps, kpardi(1), kpardi(2), delta_phi/pi);
        end
        fourth_line = sprintf('$t_i = %.0f T, t_f = (%.0f, %4.3f T; %4.3f T), (n_{v_\\perp}, n_{\\theta}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$', ...
            t_init/waveT, t_final(1), t_final(end-1)/waveT, dt_final/waveT, nvperp, ntheta, nvz);

        if keep_anno
        text(-3, 3, {[first_line], [second_line], [third_line], [fourth_line]}, ...
            'Interpreter','latex', 'VerticalAlignment','top', FontSize=15);
        end

        title("$C_{E_z}(v_\parallel, v_\perp)$", ...
            'Interpreter','latex', ...
            'FontName','TimesNewRoman', ...
            'FontSize', 32, ...
            'FontWeight','bold');
        % title({[first_line], ...
        %     [second_line]},'Interpreter','latex', ...
        %     'FontName','TimesNewRoman', ...
        %     'FontSize',16, ...
        %     'FontWeight','bold');
        
        % xlabel('$v_\parallel/v_{ti}$','Interpreter','latex', ...
        %     'FontName','TimesNewRoman', ...
        %     'FontSize',24, ...
        %     'FontWeight','bold');
        
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

        text(-3, max(tmpf_1DCez)*0.97, "(f)", ...
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
figure_name = sprintf("./plots_iSH/iSHLDV12_%s", time_suffix);
% print(h10, figure_name, '-dpdf', '-vector');
exportgraphics(h10, figure_name, "ContentType","vector");
print(h10, figure_name, '-painters','-depsc','-r150');
            else
            figure_name = sprintf("./plots_iSH/iSHLDV12_%s", time_suffix);
            print(h10, figure_name, '-dpng', '-r150');  % 150 dpi resolution
            end
        end
    case 13

        % load('~/Desktop/iSHLDV132_20241126_173925.mat');
        load(datapath);

        % if exist("delta_phi") == 0
        %    delta_phi = 0.;
        % end

        [~, ~, ~, ~, ~, ~, ~, field_choice_local, ~, ~, ~, ~, delta_phi_local, kpardi] = set_params;

        if field_choice == -495 || field_choice == -4951
            if field_choice_local == field_choice && delta_phi == delta_phi_local
                disp("Correct local field.");
            else
                disp("Wrong local field. Please press Ctrl+C.");
                pause;
            end
        else
            if field_choice_local == field_choice
                disp("Correct local field.");
            else
                disp("Wrong local field. Please press Ctrl+C.");
                pause;
            end
        end

        % avg_cez = squeeze(sum(cez(:, :, :, :, :, 1:end-1), 6) * dtau / (tau(end) - tau(1)));
        
        avg_cez = squeeze(sum(cez(:, :, :, :, :, 1:end-1), 6) * dtf / (tf_values(end) - tf_values(1)));
        % Create the meshgrid for plotting
        [VZ, ti] = meshgrid(vz_values, t_multiples);
        
        h10 = figure('Position', [1, 200, 1400, 900], 'Theme', 'Light');
        
        % % Check if the version is R2023a or earlier
        % if strcmp(matlab_version, '2023a') || (str2double(matlab_version(1:4)) < 2023)
        %     disp('This is MATLAB R2023a or an earlier version.');
        %     % Continue
        % else
        %     disp('This is MATLAB later than R2023a.');
        %     h10.Theme = 'Light';
        % end
        
        [~, h] = contourf(VZ, ti, avg_cez, 50);
        hold on;
        xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        % yline(3, 'LineStyle','--', 'LineWidth', 3);
        % yline(5, 'LineStyle','--', 'LineWidth', 3);
        % yline(8, 'LineStyle','--', 'LineWidth', 3);
        if output_reversal_time
            xline(input_vz, 'LineStyle','-', 'LineWidth', 1);
            % Find index of closest value in vz_values
            [~, closest_index] = min(abs(vz_values - input_vz));
            
            % Extract corresponding column from avg_cez
            avg_cez_1D = avg_cez(:, closest_index);
            
            % Define interval (comment this line if no interval is needed)
            interval = [min(t_multiples), max(t_multiples)];  % or something like [0.2, 0.5]
            
            % Call zero-crossing function
            zero_crossings = find_zero_crossings(t_multiples, avg_cez_1D, interval);
            
            % Display results
            disp('Zero crossings found at:');
            disp(zero_crossings);

            % Limit to available zero crossings
            num_to_label = min(num_to_label, length(zero_crossings));
            % Limit to the number of labeled crossings
            zc_subset = zero_crossings(1:num_to_label);
            differences = diff(zc_subset);  % differences between successive elements
            
            % --- Step 2: Output the differences ---
            disp('Differences between adjacent zero crossings:');
            disp(differences);
            text(-1, 60, 'Reversal Times:', ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 20, 'Color', 'k');
            disp_diff = flip(differences);
            for i = 1:length(differences)
                x = disp_diff(i);
                text(-0.5, 60-2.2*i, sprintf('%.3f', x), ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 20, 'Color', 'k');

            end
            
            % --- Step 3: Output the mean of the differences ---
            mean_diff = mean(differences);
            fprintf('Mean difference: %.4f\n', mean_diff);

            text(-1, 30, sprintf('Mean Reversal Time: %.3f', mean_diff), ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 20, 'Color', 'k');
            
            for i = 1:num_to_label
                x = zero_crossings(i);
                % Plot vertical dashed line or a point
                plot(input_vz, x, 'ko', 'MarkerSize', 6, 'LineWidth', 2);  % red marker at zero crossing
                % Label it
                text(input_vz, x, sprintf('%.3f', x), ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', 20, 'Color', 'k');
            end
            text(input_vz, 0, sprintf('%.3f', input_vz), ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 15, 'Color', 'k');
            
            hold off;

        end

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
        
        % daspect([1 1 1]);
        grid on;

        switch(field_choice)
            case {-49, -45}
                first_line = sprintf("$C_{E_z}(v_\\parallel, |t_i|)$, No Window Function, 1-KAW, $k_{\\parallel} \\rho_i = %4.3f$", ...
                    kpardi);
            case {-491, -451}
                first_line = sprintf("$C_{E_z}(v_\\parallel, |t_i|)$, Window Function, 1-KAW, $k_{\\parallel} \\rho_i = %4.3f$", ...
                    kpardi);
            case -495
                first_line = sprintf("$C_{E_z}(v_\\parallel, |t_i|)$, No Window Function, 2-KAW, $(x, y, z) = (%3.2f, %3.2f, %3.2f), k_{\\parallel 1} \\rho_i = %4.3f, k_{\\parallel 2} \\rho_i = %4.3f$", ...
                    xval, yval, zval, kpardi(1), kpardi(2));
            case -4951
                first_line = sprintf("$C_{E_z}(v_\\parallel, |t_i|)$, Window Function, 2-KAW, $k_{\\parallel 1} \\rho_i = %4.3f, k_{\\parallel 2} \\rho_i = %4.3f$", ...
                    kpardi(1), kpardi(2));
        end

        second_line = sprintf("RSR = %3.2f, $\\delta \\phi = %3.2f \\pi, n_{v_z} = %.0f, t_i = (%.0f, %.0f, %3.2f)T, t_f = (%.0f, %4.3f, %4.3f)T$", ...
            em_eps, delta_phi/pi, nvz, ...
            ti_values(1)/waveT, ti_values(end)/waveT, (ti_values(2)-ti_values(1))/waveT, ...
            tf_values(1)/waveT, tf_values(end-1)/waveT, (tf_values(2)-tf_values(1))/waveT);

        title({[first_line], ...
            [second_line]},'Interpreter','latex', ...
            'FontName','TimesNewRoman', ...
            'FontSize',20, ...
            'FontWeight','bold');
        
        xlabel('$v_\parallel/v_{ti}$','Interpreter','latex', ...
            'FontName','TimesNewRoman', ...
            'FontSize',24, ...
            'FontWeight','bold');
        
        ylabel('$|t_i|/T$','Interpreter','latex', ...
            'FontName','TimesNewRoman', ...
            'FontSize',24, ...
            'FontWeight','bold');
        
        colormap(bluewhitered);
        if save_plot
            if output_reversal_time
            figure_name = sprintf("./plots_iSH/iSHLDV132_%s_LabelReversalTime", time_suffix);
            else
            figure_name = sprintf("./plots_iSH/iSHLDV132_%s", time_suffix);
            end
            print(h10, figure_name, '-dpng', '-r150');  % 150 dpi resolution
        end

end


function zero_crossings = find_zero_crossings(t_multiples, avg_cez_1D, interval)
    zero_crossings = [];
    
    % Check if an interval is provided
    if nargin == 3
        mask = (t_multiples >= interval(1)) & (t_multiples <= interval(2));
        t_multiples = t_multiples(mask);
        avg_cez_1D = avg_cez_1D(mask);
    end
    
    % Loop through the avg_cez_1D values to find where sign changes occur
    for i = 1:length(avg_cez_1D)-1
        if sign(avg_cez_1D(i)) ~= sign(avg_cez_1D(i+1))
            % Linear interpolation to find a more accurate zero-crossing
            t1 = t_multiples(i);
            t2 = t_multiples(i+1);
            C1 = avg_cez_1D(i);
            C2 = avg_cez_1D(i+1);
            
            % Calculate the zero crossing using interpolation
            t_zero = t1 + (0 - C1) * (t2 - t1) / (C2 - C1);
            
            % Add the zero-crossing point to the result
            zero_crossings = [zero_crossings; t_zero];
        end
    end
end
