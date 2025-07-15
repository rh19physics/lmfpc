% This script reads iSHCDV21 data and calculate FPC quantities
% To reproduce figures in the paper, choose plotting_option = 2
close all
clear
clc

save_figure = true;

figure_number = 1;

switch(figure_number)
    case 1
        anno_labels = ["$(a)$", "$(b)$", "$(c)$", "$(d)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20241117_185132.mat");
    case 2
        anno_labels = ["$(e)$", "$(f)$", "$(g)$", "$(h)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20241117_185132.mat");
    case 3
        anno_labels = ["$(i)$", "$(j)$", "$(k)$", "$(l)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20241117_185132.mat");
end

keep_anno = false;
output_gif = false;
plot_logf = true;
logdyn = 6;
plotting_option = 2; % if 1, plot using subplot; if 2, plot using tilelayout
% if 0 or negative numbers, plot whatever you like
% -1: plot it in 3D
plot_2panels = false;
vertical_layout = false;



[~, ~, ~, ~, bi, ~, ~, field_choice_local, ~, ~, t_init_local, delta_phi_local, kpardi] = set_params;

if field_choice_local == field_choice && delta_phi == delta_phi_local && t_init == t_init_local
    disp("Correct local field.");
else
    disp("Wrong local field. Please press Ctrl+C.");
    pause;
end

if t_final(end) / waveT > 1.2
    chosen_tf_slices = [1, 10, 40, 70];
else
    chosen_tf_slices = [1, 10, 20, 30];
end

if field_choice == -653 || field_choice == -6531
    plot_circs = false;
else
    plot_circs = true;
end

% before 2024-11-28 I used xval for x, y, z
% later, I decided to choose 3 different spatial coordinates
% to make this plotting script compatible to old datasets
% I have to add:

if exist("yval") 
else
    yval = xval;
    zval = xval;
end

% =============== Extra Calcultions =================
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

f_VxVy = sum(f, 4) * dvz;

ceperp = cex + cey;
tavg_cex = sum(cex(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
tavg_cey = sum(cey(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));
tavg_ceperp = sum(ceperp(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));

tavg_cex_VxVz = sum(tavg_cex, 3) * dvy;
tavg_cex_VyVz = sum(tavg_cex, 2) * dvx;

tavg_cey_VxVz = sum(tavg_cey, 3) * dvy;
tavg_cey_VyVz = sum(tavg_cey, 2) * dvx;

tavg_ceperp_VxVz = sum(tavg_ceperp, 3) * dvy;
tavg_ceperp_VyVz = sum(tavg_ceperp, 2) * dvx;

ceperp_VxVy = sum(ceperp, 4) * dvz;
tavg_ceperp_VxVy = sum(ceperp_VxVy(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));

ceperpV2 = zeros(nx, nvx, nvy, nvz, nt_final);

for it_final = 1:nt_final
    for ix = 1:nx
        E0 = elecfield(t_final(it_final), [xval(ix), yval(ix), zval(ix)]);
        for ivx = 1:nvx
            for ivy = 1:nvy
                for ivz = 1:nvz
                    ceperpV2(ix, ivx, ivy, ivz, it_final) = - q * (vx_values(ivx)^2. ...
                        + vy_values(ivy)^2. + vz_values(ivz)^2.) / 2. * (E0(1) * dfdvx(ix, ivx, ivy, ivz, it_final) + ...
                        E0(2) * dfdvy(ix, ivx, ivy, ivz, it_final));
                end
            end
        end
    end
end

ceperpV2_VxVy = sum(ceperpV2, 4) * dvz;
tavg_ceperpV2_VxVy = sum(ceperpV2_VxVy(:, :, :, :, 1:end-1), 5) * dt_final / (t_final(end) - t_final(1));


% Compute vperp for each (vx, vy) pair
vperp_values = [];
for i = 1:nvx
    for j = 1:nvy
        vx = vx_values(i);
        vy = vy_values(j);
        vperp = sqrt(vx^2. + vy^2.);
        vperp_values = [vperp_values; vperp];
    end
end

% Get unique values of vperp and sort them from small to large
vperp_unique = unique(vperp_values);
nvperp = length(vperp_unique);

ceperp_VperpVz = zeros(nx, nvperp, nvz, nt_final);
f_VperpVz = zeros(nx, nvperp, nvz, nt_final);

% Loop over time, spatial, and velocity coordinates
for it_final = 1:nt_final
    for ix = 1:nx
        for ivz = 1:nvz
            for ivp = 1:nvperp
                % Find indices where vperp matches vperp_unique(ivp)
                matching_indices = find(abs(vperp_values - vperp_unique(ivp)) < 1e-6);

                sum_value = 0;
                for idx = 1:length(matching_indices)-1
                    [ii, jj] = ind2sub([nvx, nvy], matching_indices(idx));
                    sum_value = sum_value + ceperp(ix, jj, ii, ivz, it_final);
                end
                ceperp_VperpVz(ix, ivp, ivz, it_final) = sum_value / length(matching_indices);

                % transform f
                sum_value = 0; 
                for idx = 1:length(matching_indices)-1
                    [ii, jj] = ind2sub([nvx, nvy], matching_indices(idx));
                    sum_value = sum_value + f(ix, jj, ii, ivz, it_final);
                end
                f_VperpVz(ix, ivp, ivz, it_final) = sum_value/length(matching_indices);
            end
        end
    end
end

tavg_ceperp_VperpVz = sum(ceperp_VperpVz(:, :, :, 1:end-1), 4) * dt_final / (t_final(end) - t_final(1));
% tavg_f_VperpVz = sum(f_VperpVz(:, :, :, 1:end-1), 4) * dt_final / (t_final(end) - t_final(1));

[VZ, VPERP] = meshgrid(vz_values, vperp_unique);
[VZ_ZX, VX_ZX] = meshgrid(vz_values, vx_values);
[VZ_ZY, VY_ZY] = meshgrid(vz_values, vy_values);

% =============== Plotting =================
% Parameters
if plot_circs
    th_cirs = linspace(0, pi, 100);   % Angular range for semicircles
    nth_cirs = length(th_cirs);
    % radii = linspace(0.5, 2, 5); % Define a range of radii for the family of circles
    radii = [0.5, 1, 1.5, 2, 2.5, 3, 3.5];
    nR = length(radii);
    center_x = 2. * pi / waveT / kpardi;           % x-coordinate of circle centers
    n_pos1_mode = (2. * pi / waveT - 1.) / kpardi;
    n_neg1_mode = (2. * pi / waveT + 1.) / kpardi;
    
    x_cirs = zeros(nR);
    y_cirs = zeros(nR);
    
    for iR = 1:nR
        for ith_cirs = 1:nth_cirs
            x_cirs(iR, ith_cirs) = radii(iR) .* cos(th_cirs(ith_cirs)) + center_x;
            y_cirs(iR, ith_cirs) = radii(iR) .* sin(th_cirs(ith_cirs));
        end
    end
else
    center_x1 = 2. * pi / waveT / kpardi(1);           % x-coordinate of circle centers
    n_pos1_mode1 = (2. * pi / waveT - 1.) / kpardi(1);
    n_neg1_mode1 = (2. * pi / waveT + 1.) / kpardi(1);
    center_x2 = 2. * pi / waveT / kpardi(2);           % x-coordinate of circle centers
    n_pos1_mode2 = (2. * pi / waveT - 1.) / kpardi(2);
    n_neg1_mode2 = (2. * pi / waveT + 1.) / kpardi(2);
end

switch plotting_option
    case 2 % Use tilelayout
        scrsz = get(0, "ScreenSize");
    if plot_2panels
        hLF2 = figure('Position', [1 scrsz(3) 0.6*scrsz(3) 0.25*scrsz(3)]);
        t = tiledlayout(1, 2);
    elseif vertical_layout
        hLF2 = figure('Position', [1 scrsz(3) 0.5*scrsz(3) 2*scrsz(3)]);
        t = tiledlayout(4, 1);
    else
        % hLF2 = figure;
        hLF2 = figure('Position', [0 0 scrsz(3) 0.25*scrsz(3)]);
        % hLF2 = figure('Position', [0 0 1920 0.25*1920]);
        t = tiledlayout(1, 4);
    end
        
        
        % =========================================
        ax1 = nexttile;
        % Plot f_VperpVz at the 2nd time slices
        tmpf = squeeze(f_VperpVz);
        if plot_logf
            tmpf = log10(tmpf);
            lgmax=max(tmpf,[],'all');
            tmpf(tmpf< lgmax-logdyn)=lgmax-logdyn;
        end
        
        [~, h1] = contourf(VZ, VPERP, tmpf(:, :, end-1), 50);
        hold on;
        annotation('textbox', [0.05, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(1), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);
        annotation('textbox', [0.28, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(2), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);
        annotation('textbox', [0.54, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(3), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);
        annotation('textbox', [0.77, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(4), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);

        set(h1,'edgecolor','none');
        colorbar('FontSize',20,'FontName','TimesNewRoman');
        if plot_circs
            for iR = 1:nR
                plot(x_cirs(iR, :), y_cirs(iR, :), 'LineStyle','-', 'LineWidth', 1, 'Color','#6aa84f');
            end
        end
        if field_choice == -653 || field_choice == -6531
            xline(n_pos1_mode1, 'LineStyle','--', 'LineWidth', 2);
            xline(center_x1, 'LineStyle',':', 'LineWidth', 2);
            xline(n_pos1_mode2, 'LineStyle','--', 'LineWidth', 2);
            xline(center_x2, 'LineStyle',':', 'LineWidth', 2);
            yline(1, 'LineStyle',':', 'LineWidth', 2);
        else        
            xline(n_pos1_mode, 'LineStyle','--', 'LineWidth', 2);
            xline(center_x, 'LineStyle',':', 'LineWidth', 2);
            xline(n_neg1_mode, 'LineStyle','--', 'LineWidth', 2);
            yline(1, 'LineStyle',':', 'LineWidth', 2);
        end
        
        grid on;
        daspect([1 1 1]);
        
        title_label = sprintf("$\\log_{10}f(v_\\parallel, v_\\perp)$ at $t_f = %4.3f T$", t_final(end-1)/waveT);
        
        % format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$', title_label);
     format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$', "$\log_{10}f(v_\parallel, v_\perp)$");
        colormap(ax1, plasma);

        % =========================================
        ax2 = nexttile;
            tmpf = squeeze(tavg_ceperp_VperpVz);
            [~, h2] = contourf(VZ, VPERP, tmpf, 50);
            hold on;
            if plot_circs
                for iR = 1:nR
                    plot(x_cirs(iR, :), y_cirs(iR, :), 'LineStyle','-', 'LineWidth', 1, 'Color','#6aa84f');
                end
            end
            if field_choice == -653 || field_choice == -6531
                xline(n_pos1_mode1, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x1, 'LineStyle',':', 'LineWidth', 2);
                xline(n_pos1_mode2, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x2, 'LineStyle',':', 'LineWidth', 2);
                yline(1, 'LineStyle',':', 'LineWidth', 2);
            else        
                xline(n_pos1_mode, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x, 'LineStyle',':', 'LineWidth', 2);
                xline(n_neg1_mode, 'LineStyle','--', 'LineWidth', 2);
                yline(1, 'LineStyle',':', 'LineWidth', 2);
            end
            
            set(h2,'edgecolor','none');
            colorbar('FontSize', 20,'FontName','TimesNewRoman');
            
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis(ax2, [-max(abs(cL)) max(abs(cL))]); 
            
            daspect([1 1 1]);
            format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$', '$C_{E_\perp}(v_\parallel, v_\perp)$');
            colormap(ax2, bluewhitered);

        
        % =========================================
if plot_2panels == false        
        ax3 = nexttile;
        tmpf = squeeze(tavg_cex_VxVy);
        [~, h3] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h3,'edgecolor','none');
        colorbar('FontSize', 20,'FontName','TimesNewRoman');
        
        % Sets the colormap limits such that the center of the color bar is 0
        cL = caxis;  
        caxis(ax3, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_x}(v_x, v_y)$');
        colormap(ax3, bluewhitered);
            
        % =========================================
        ax4 = nexttile;
       tmpf = squeeze(tavg_cey_VxVy);
        [~, h4] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h4,'edgecolor','none');
        colorbar('FontSize', 20,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        % Sets the colormap limits such that the center of the color bar is 0
        cL = caxis;  
        caxis(ax4, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_y}(v_x, v_y)$');
        
        colormap(ax4, bluewhitered);
end
        % Add a super title
        if field_choice == -653 || field_choice == -6531
            params_line = sprintf("$\\mathbf{r} = (%3.2f, %3.2f, %3.2f), RSR = %3.2f, \\delta \\phi = %3.2f \\pi, k_{\\parallel, 1} \\rho_i = %4.3f, k_{\\parallel, 2} \\rho_i = %4.3f, t_i = %1.1d T, t_f = (%1.1d, %4.3f T; %4.3f T), (n_{v_x}, n_{v_y}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$", ...
                xval, yval, zval, em_eps, delta_phi/pi, kpardi(1), kpardi(2), t_init/waveT, t_final(1), t_final(end-1)/waveT, dt_final/waveT, nvx, nvy, nvz);
        else
            params_line = sprintf("$\\beta_i = %1.1f, RSR = %3.2f, k_{\\parallel} \\rho_i = %4.3f, t_i = %1.1d T, t_f = (%1.1d, %4.3f T; %4.3f T), (n_{v_x}, n_{v_y}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$", ...
                bi, em_eps, kpardi, t_init/waveT, t_final(1), t_final(end-1)/waveT, dt_final/waveT, nvx, nvy, nvz);
        end

        if keep_anno
            sgtitle(params_line, 'Interpreter', 'latex', 'FontSize', 15);
        end
        
        if save_figure
            if plot_2panels
                figure_filename = sprintf('./iSHCDV21_%s_LF22', time_suffix);
            else
                figure_filename = sprintf('./iSHCDV21_%s_LF2', time_suffix);
            end
            % Match the figure size (convert pixels to inches, or set directly)
            fig_width = 18;   % in inches 18
            fig_height = 4.5;   % in inches 4.5
            set(hLF2, 'PaperSize', [fig_width fig_height], 'PaperPosition', [0 0 fig_width fig_height]);
            print(hLF2, figure_filename, '-dpdf', '-vector');
            print(hLF2, figure_filename, '-painters','-depsc','-r150');
        end

    case 1 % Use subplot
        h10 = figure('Position', [1, 1, 1500, 1500], 'Visible','on', 'Theme', 'Light');
        
        ax1 = subplot(3, 3, 1);
        % Plot tavg_cex_VxVy
        tmpf(:,:) = squeeze(tavg_cex_VxVy);
        [~, h1] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h1,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        % Sets the colormap limits such that the center of the color bar is 0
        cL = caxis;  
        caxis(ax1, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_x}(v_x, v_y)$');
        colormap(ax1, bluewhitered);
        
        % =========================================
        ax2 = subplot(3, 3, 2);
        % Plot tavg_cey_VxVy
        tmpf(:,:) = squeeze(tavg_cey_VxVy);
        [~, h2] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h2,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        % Sets the colormap limits such that the center of the color bar is 0
        cL = caxis;  
        caxis(ax2, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_y}(v_x, v_y)$');
        
        colormap(ax2, bluewhitered);
        %======================================
        ax3 = subplot(3, 3, 3);
        % Plot tavg_cex_VxVy
        % plot_f = squeeze(tavg_cex_VxVy) + squeeze(tavg_cey_VxVy);
        plot_f = squeeze(tavg_ceperp_VxVy);
        [~, h3] = contourf(VX, VY, transpose(plot_f), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h3,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        % Sets the colormap limits such that the center of the color bar is 0
        cL = caxis;  
        caxis(ax3, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_\perp}(v_x, v_y)$');
        colormap(ax3, bluewhitered);
        
        %======================================
        ax4 = subplot(3, 3, 4);
        % Plot ceperpV2_VxVy
        tmpf = squeeze(tavg_ceperpV2_VxVy);
        [~, h4] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h4,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        % Sets the colormap limits such that the center of the color bar is 0
        cL = caxis;  
        caxis(ax4, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_\perp}^{(v^2)}(v_x, v_y)$');
        colormap(ax4, bluewhitered);
        
        %======================================
        ax5 = subplot(3, 3, 5);
        if nvz == 1
            % continue
            tmpf = squeeze(tavg_ceperp_VperpVz);
            plot(vperp_unique, tmpf, 'LineWidth', 3);
            hold on;
            xline(1, 'LineStyle',':', 'LineWidth', 2);
            format_subplot('$v_\perp/v_{ti}$', '$C_{E_\perp}$', '1D $C_{E_\perp}(v_\perp)$ Plot');
            vz_slice_label = sprintf("$v_z = %3.2f$", vz_values);
            text(0.7, 0.85, vz_slice_label, 'FontSize', 16,'Interpreter','latex','FontWeight','bold', 'Color','k', 'Units', 'normalized');
            ylim([- max(abs(tmpf)), max(abs(tmpf))]);
        else
            % Plot ceperpV2_VxVy
            tmpf = squeeze(tavg_ceperp_VperpVz);
            [~, h4] = contourf(VZ, VPERP, tmpf, 50);
            hold on;
            if plot_circs
                for iR = 1:nR
                    plot(x_cirs(iR, :), y_cirs(iR, :), 'LineStyle','-', 'LineWidth', 1, 'Color','#6aa84f');
                end
            end
            if field_choice == -653 || field_choice == -6531
                xline(n_pos1_mode1, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x1, 'LineStyle',':', 'LineWidth', 2);
                xline(n_pos1_mode2, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x2, 'LineStyle',':', 'LineWidth', 2);
                yline(1, 'LineStyle',':', 'LineWidth', 2);
            else        
                xline(n_pos1_mode, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x, 'LineStyle',':', 'LineWidth', 2);
                xline(n_neg1_mode, 'LineStyle','--', 'LineWidth', 2);
                yline(1, 'LineStyle',':', 'LineWidth', 2);
            end
            
            set(h4,'edgecolor','none');
            colorbar('FontSize',16,'FontName','TimesNewRoman');
            
            % xlim([-3,3]);
            % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
            
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis(ax5, [-max(abs(cL)) max(abs(cL))]); 
            
            daspect([1 1 1]);
            format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$', '$C_{E_\perp}(v_\parallel, v_\perp)$');
            colormap(ax5, bluewhitered);
        end
        % =========================================
        ax6 = subplot(3, 3, 6);
        % Plot f_VxVy at the 1st time slices
        tmpf = squeeze(f_VxVy);
        
        
        if plot_logf
            tmpf = log10(tmpf);
            lgmax=max(tmpf,[],'all');
            tmpf(tmpf< lgmax-logdyn)=lgmax-logdyn;
        end
        
        [~, h3] = contourf(VX, VY, transpose(tmpf(:, :, chosen_tf_slices(1))), 50);
        hold on;
        
        set(h3,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        if plot_logf
            colormap(ax6, plasma);
            title_label = sprintf("$\\log_{10} f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(1))/waveT);
        else
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis(ax6, [0, max(abs(cL))]);
            % colormap
            colormap(ax6, bluewhitered);
            title_label = sprintf("$f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(1))/waveT);
        end
        
        daspect(ax6, [1 1 1]);
        grid on;
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', title_label);
        
        % =========================================
        ax7 = subplot(3, 3, 7);
        % Plot f_VxVy at the 1st time slices
        tmpf = squeeze(f_VxVy);
        
        
        if plot_logf
            tmpf = log10(tmpf);
            lgmax=max(tmpf,[],'all');
            tmpf(tmpf< lgmax-logdyn)=lgmax-logdyn;
        end
        
        [~, h3] = contourf(VX, VY, transpose(tmpf(:, :, chosen_tf_slices(2))), 50);
        hold on;
        
        set(h3,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        if plot_logf
            colormap(ax7, plasma);
            title_label = sprintf("$\\log_{10} f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(2))/waveT);
        else
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis(ax7, [0, max(abs(cL))]);
            % colormap
            colormap(ax7, bluewhitered);
            title_label = sprintf("$f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(2))/waveT);
        end
        
        
        daspect([1 1 1]);
        grid on;
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', title_label);
        % =========================================
        ax8 = subplot(3, 3, 8);
        % Plot f_VxVy at the 1st time slices
        tmpf = squeeze(f_VxVy);
        
        
        if plot_logf
            tmpf = log10(tmpf);
            lgmax=max(tmpf,[],'all');
            tmpf(tmpf< lgmax-logdyn)=lgmax-logdyn;
        end
        
        [~, h3] = contourf(VX, VY, transpose(tmpf(:, :, chosen_tf_slices(3))), 50);
        hold on;
        
        set(h3,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        if plot_logf
            colormap(ax8, plasma);
            title_label = sprintf("$\\log_{10} f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(3))/waveT);
        else
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis(ax8, [0, max(abs(cL))]);
            % colormap
            colormap(ax8, bluewhitered);
            title_label = sprintf("$f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(3))/waveT);
        end
        
        daspect([1 1 1]);
        grid on;
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', title_label);
        
        % =========================================
        ax9 = subplot(3, 3, 9);
        % Plot f_VxVy at the 1st time slices
        tmpf = squeeze(f_VxVy);
        
        if plot_logf
            tmpf = log10(tmpf);
            lgmax=max(tmpf,[],'all');
            tmpf(tmpf< lgmax-logdyn)=lgmax-logdyn;
        end
        
        [~, h3] = contourf(VX, VY, transpose(tmpf(:, :, chosen_tf_slices(4))), 50);
        hold on;
        
        set(h3,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        
        % xlim([vzmin - dv / 2, vzmax + dv / 2]);
        % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
        
        
        if plot_logf
            colormap(ax9, plasma);
            title_label = sprintf("$\\log_{10}f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(4))/waveT);
        else
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis(ax9, [0, max(abs(cL))]);
            % colormap
            colormap(ax9, bluewhitered);
            title_label = sprintf("$f(v_x, v_y)$ at $t_f = $ %3.2f T", t_final(chosen_tf_slices(4))/waveT);
        end
        
        daspect([1 1 1]);
        grid on;
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', title_label);
        
        
        % Add a super title
        if field_choice == -653 || field_choice == -6531
            params_line = sprintf("$\\mathbf{r} = (%3.2f, %3.2f, %3.2f), RSR = %3.2f, \\delta \\phi = %3.2f \\pi, k_{\\parallel, 1} \\rho_i = %4.3f, k_{\\parallel, 2} \\rho_i = %4.3f, t_i = %1.1d T, t_f = (%1.1d, %1.1d T; %4.3f T), (n_{v_x}, n_{v_y}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$", ...
                xval, yval, zval, em_eps, delta_phi/pi, kpardi(1), kpardi(2), t_init/waveT, t_final(1), t_final(end)/waveT, dt_final/waveT, nvx, nvy, nvz);
        else
            params_line = sprintf("$\\beta_i = %1.1f, RSR = %3.2f, k_{\\parallel} \\rho_i = %4.3f, t_i = %1.1d T, t_f = (%1.1d, %1.1d T; %4.3f T), (n_{v_x}, n_{v_y}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$", ...
                bi, em_eps, kpardi, t_init/waveT, t_final(1), t_final(end)/waveT, dt_final/waveT, nvx, nvy, nvz);
        end
        
        sgtitle("iCD Signature Overview"  + newline + params_line, ...
                'Interpreter', 'latex', 'FontSize', 16);
        
        if save_figure
            figure_filename = sprintf('./iSHCDV21_%s.png', time_suffix);
            print(h10, figure_filename, '-dpng', '-r150');  % 150 dpi resolution
        end

    case 0 % testing case
        hLF0 = figure;
        t = tiledlayout(3, 3, "TileSpacing","compact");
        

        % ========== Cex =============================
        ax1 = nexttile;
        % Plot tavg_cex_VxVy
        tmpf = squeeze(tavg_cex_VxVy);
        [~, h1] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h1,'edgecolor','none');
        colorbar('FontSize', 20,'FontName','TimesNewRoman');
        
        cL = caxis;  
        caxis(ax1, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_x}(v_x, v_y)$');
        colormap(ax1, bluewhitered);

        
        ax2 = nexttile; % Plot tavg_cex_VxVz
        tmpf = squeeze(tavg_cex_VxVz);
        [~, h2] = contourf(VZ_ZX, VX_ZX, tmpf, 50);

        hold on;
        set(h2,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        cL = caxis;  
        caxis(ax2, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_z/v_{ti}$', '$v_x/v_{ti}$', '$C_{E_x}(v_z, v_x)$');
        
        colormap(ax2, bluewhitered);

        ax3 = nexttile; % Plot tavg_cex_VyVz
        tmpf = squeeze(tavg_cex_VyVz);
        [~, h3] = contourf(VZ_ZY, VY_ZY, tmpf, 50);

        hold on;
        set(h3,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        cL = caxis;  
        caxis(ax3, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_z/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_x}(v_z, v_y)$');
        
        colormap(ax3, bluewhitered);

        % ========== Cey =============================
        ax4 = nexttile;
        % Plot tavg_cey_VxVy
        tmpf = squeeze(tavg_cey_VxVy);
        [~, h4] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h4,'edgecolor','none');
        colorbar('FontSize', 20,'FontName','TimesNewRoman');
        
        cL = caxis;  
        caxis(ax4, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_y}(v_x, v_y)$');
        colormap(ax4, bluewhitered);

        
        ax5 = nexttile; % Plot tavg_cey_VxVz
        tmpf = squeeze(tavg_cey_VxVz);
        [~, h5] = contourf(VZ_ZX, VX_ZX, tmpf, 50);

        hold on;
        set(h5,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        cL = caxis;  
        caxis(ax5, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_z/v_{ti}$', '$v_x/v_{ti}$', '$C_{E_y}(v_z, v_x)$');
        
        colormap(ax5, bluewhitered);

        ax6 = nexttile; % Plot tavg_cey_VyVz
        tmpf = squeeze(tavg_cey_VyVz);
        [~, h6] = contourf(VZ_ZY, VY_ZY, tmpf, 50);

        hold on;
        set(h6,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        cL = caxis;  
        caxis(ax6, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_z/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_y}(v_z, v_y)$');
        
        colormap(ax6, bluewhitered);

        % ========== Ceperp =============================
        ax7 = nexttile;
        % Plot tavg_ceperp_VxVy
        tmpf = squeeze(tavg_ceperp_VxVy);
        [~, h7] = contourf(VX, VY, transpose(tmpf), 50);
        hold on;
        % xline(1.135, 'LineStyle','--', 'LineWidth', 3);
        % xline(-1.135, 'LineStyle','--', 'LineWidth', 3);
        
        set(h7,'edgecolor','none');
        colorbar('FontSize', 20,'FontName','TimesNewRoman');
        
        cL = caxis;  
        caxis(ax7, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        
        format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_\perp}(v_x, v_y)$');
        colormap(ax7, bluewhitered);

        
        ax8 = nexttile; % Plot tavg_ceperp_VxVz
        tmpf = squeeze(tavg_ceperp_VxVz);
        % tmpf = squeeze(tavg_cex_VxVz+tavg_cey_VxVz);
        [~, h8] = contourf(VZ_ZX, VX_ZX, tmpf, 50);

        hold on;
        set(h8,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        cL = caxis;  
        caxis(ax8, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_z/v_{ti}$', '$v_x/v_{ti}$', '$C_{E_\perp}(v_z, v_x)$');
        
        colormap(ax8, bluewhitered);

        ax9 = nexttile; % Plot tavg_ceperp_VyVz
        tmpf = squeeze(tavg_ceperp_VyVz);
        [~, h9] = contourf(VZ_ZY, VY_ZY, tmpf, 50);

        hold on;
        set(h9,'edgecolor','none');
        colorbar('FontSize',16,'FontName','TimesNewRoman');
        cL = caxis;  
        caxis(ax9, [-max(abs(cL)) max(abs(cL))]); 
        
        daspect([1 1 1]);
        format_subplot('$v_z/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_\perp}(v_z, v_y)$');
        
        colormap(ax9, bluewhitered);
    case -1
        disp_meth = 2;
        [VX_XYZ, VY_XYZ, VZ_XYZ] = meshgrid(vx_values, vy_values, vz_values);
        tmpf = squeeze(tavg_ceperp);
        flatten_tmpf = tmpf(:);
        figure;
    if disp_meth == 1
        VX_XYZ = VX_XYZ(:); % Convert to column vector
        VY_XYZ = VY_XYZ(:);
        VZ_XYZ = VZ_XYZ(:);
        scatter3(VX_XYZ, VY_XYZ, VZ_XYZ, 36, flatten_tmpf, 'filled');
    elseif disp_meth == 2
        if output_gif        
        % Define the output GIF filename
        gif_filename = '3D_rotation_525.gif';
        end
        slice(VX_XYZ, VY_XYZ, VZ_XYZ,tmpf, ...
            mean(vx_values), ...
            mean(vy_values), ...
            []);
        hold on;
        % Define the tilt angle (adjust as needed)
        angles = [pi/3, pi/6, -pi/6, -pi/3]; % 45-degree tilt
        
        for iangles = 1:length(angles)
        % Generate a grid for the slicing surface
        [xsurf, zsurf] = meshgrid(linspace(min(vx_values), max(vx_values), 56), ...
                                  linspace(min(vz_values), max(vz_values), 40));
                              
        % Compute z-values for the tilted plane using the plane equation
        ysurf = xsurf * tan(angles(iangles));
        slice(VX_XYZ, VY_XYZ, VZ_XYZ, tmpf, xsurf, ysurf, zsurf); % Add the tilted slice
        end
        % shading interp;
    elseif disp_meth == 3
        isosurface(VX_XYZ, VY_XYZ, VZ_XYZ,tmpf, mean(flatten_tmpf));
        hold on;
        % isosurface(VX_XYZ, VY_XYZ, VZ_XYZ,tmpf, mean(flatten_tmpf));
        isosurface(VX_XYZ, VY_XYZ, VZ_XYZ,tmpf, -mean(flatten_tmpf));
        % isosurface(VX_XYZ, VY_XYZ, VZ_XYZ,tmpf, 0);
        view(3);
        camlight;
        lighting phong;
    end
        colorbar;
        cL = caxis;
        caxis([-max(abs(cL)) max(abs(cL))]);
        daspect([1 1 1]);
        colormap("bluewhitered");
        xlim([-4, 4]);
        ylim([-4, 4]);
        zlim([-4, 4]);
        xlabel('VX-axis');
        ylabel('VY-axis');
        zlabel('VZ-axis');
        title('$C_{E_\perp}(v_x, v_y, v_z)$', "Interpreter","latex");
        grid on;
if output_gif
% Set up GIF writing
frame_delay = 0.1; % Time per frame (seconds)
num_frames = 36; % Number of frames for a smooth rotation
angles = linspace(0, 360, num_frames); % Define rotation angles

for i = 1:num_frames
    view(angles(i), 30); % Rotate the view (azimuth, elevation)
    drawnow; % Update the figure

    % Capture the frame
    frame = getframe(gcf);
    img = frame2im(frame);
    [A, map] = rgb2ind(img, 256); % Convert to indexed image

    % Write to GIF file
    if i == 1
        imwrite(A, map, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', frame_delay);
    else
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end

disp(['GIF saved as ', gif_filename]);
end

    case -2
        figure;
            tmpf = squeeze(tavg_ceperp_VperpVz);
            [~, h3] = contourf(VZ, VPERP, tmpf, 50);
            hold on;
            if plot_circs
                for iR = 1:nR
                    plot(x_cirs(iR, :), y_cirs(iR, :), 'LineStyle','-', 'LineWidth', 1, 'Color','#6aa84f');
                end
            end
            if field_choice == -653 || field_choice == -6531
                xline(n_pos1_mode1, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x1, 'LineStyle',':', 'LineWidth', 2);
                xline(n_pos1_mode2, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x2, 'LineStyle',':', 'LineWidth', 2);
                yline(1, 'LineStyle',':', 'LineWidth', 2);
            else        
                xline(n_pos1_mode, 'LineStyle','--', 'LineWidth', 2);
                xline(center_x, 'LineStyle',':', 'LineWidth', 2);
                xline(n_neg1_mode, 'LineStyle','--', 'LineWidth', 2);
                yline(1, 'LineStyle',':', 'LineWidth', 2);
            end
            
            set(h3,'edgecolor','none');
            colorbar('FontSize', 20,'FontName','TimesNewRoman');
            
            xlim([-4,4]);
            % ylim([min(vx_values) - dv / 2, max(vx_values) + dv / 2]);
            
            % Sets the colormap limits such that the center of the color bar is 0
            cL = caxis;  
            caxis([-max(abs(cL)) max(abs(cL))]); 
            
            daspect([1 1 1]);
            format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$', '$C_{E_\perp}(v_\parallel, v_\perp)$');
            colormap(bluewhitered);

end


% ================ functions ================
% formatting plots
function format_subplot(xlabel_text, ylabel_text, title_text)
    grid on;
    set(gca, 'Fontsize', 24, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold', 'LineWidth', 2);
    xlabel(xlabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel(ylabel_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 24, 'FontWeight', 'bold');
    title(title_text, 'Interpreter', 'latex', 'FontName', 'TimesNewRoman', 'FontSize', 24, 'FontWeight', 'bold');
end
