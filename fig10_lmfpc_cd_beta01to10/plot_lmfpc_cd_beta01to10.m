% This script reads iSHCDV21 data, re-calculate and plot FPC quantities
% The FPC re-calculation process requires electric field. 
% So it's important to cross check the local electric field and the
% electric field used in LM

close all
clear
clc

save_figure_local = true;
keep_anno = false;
plot_logf = true;
logdyn = 6;

figure_number = 5;
switch(figure_number)
    case 1 % Set field_choice to -6301 in set_params.m
        anno_labels = ["$(a)$", "$(b)$", "$(c)$", "$(d)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20250101_113049.mat");
    case 2 % Set field_choice to -6303 in set_params.m
        anno_labels = ["$(e)$", "$(f)$", "$(g)$", "$(h)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20250101_113157.mat");
    case 3 % Set field_choice to -63 in set_params.m
        anno_labels = ["$(i)$", "$(j)$", "$(k)$", "$(l)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20241117_185132.mat");
    case 4 % Set field_choice to -633 in set_params.m
        anno_labels = ["$(m)$", "$(n)$", "$(o)$", "$(p)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20250101_113245.mat");
    case 5 % Set field_choice to -6310 in set_params.m
        anno_labels = ["$(q)$", "$(r)$", "$(s)$", "$(t)$"];
        load("../lmfpc_matlab/data/iSHCDV21_20250101_113403.mat");
end

% Cross-checking parameters from mat data and local parameters
[q_local, m_local, mime_local, tite_local, ~, ~, vtic_local, field_choice_local, em_eps_local, ~, t_init_local, delta_phi_local, kvalue_local] = set_params;

params_from_mat = struct( ...
    'q', q, ...
    'm', m, ...
    'mime', mime, ...
    'tite', tite, ...
    'vtic', vtic, ...
    'field_choice', field_choice, ...
    'em_eps', em_eps, ...
    't_init', t_init, ...
    'delta_phi', delta_phi ...
);

params_local = struct( ...
    'q', q_local, ...
    'm', m_local, ...
    'mime', mime_local, ...
    'tite', tite_local, ...
    'vtic', vtic_local, ...
    'field_choice', field_choice_local, ...
    'em_eps', em_eps_local, ...
    't_init', t_init_local, ...
    'delta_phi', delta_phi_local ...
);

disp("Start cross checking parameters...");
if isequal(params_from_mat, params_local)
    disp("Correct local field.");
else
    disp("Wrong local field. Please press Ctrl+C.");
    disp("Differences:");
    show_params_differences(params_from_mat, params_local);
    pause;
end

if size(kvalue_local, 2) == 2
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
disp("Start calculation...");
for it_final = 1:nt_final
    for ix = 1:nx
        E0 = elecfield(t_final(it_final), [xval(ix), yval(ix), zval(ix)]);    
        % ======================================================================
        % Re-calculate vx derivatives of f at each x and t
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
        % Re-calculate vy derivatives of f at each x and t
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
        % Re-calculate vz derivatives of f at each x and t
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

                % transform ceperp
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
% Plot semi-circles 
if plot_circs
    th_cirs = linspace(0, pi, 100);   % Angular range for semicircles
    nth_cirs = length(th_cirs);
    % radii = linspace(0.5, 2, 5); % Define a range of radii for the family of circles
    radii = [0.5, 1, 1.5, 2, 2.5, 3, 3.5];
    nR = length(radii);
    center_x = 2. * pi / waveT / kvalue_local;           % x-coordinate of circle centers
    n_pos1_mode = (2. * pi / waveT - 1.) / kvalue_local;
    n_neg1_mode = (2. * pi / waveT + 1.) / kvalue_local;
    
    x_cirs = zeros(nR);
    y_cirs = zeros(nR);
    
    for iR = 1:nR
        for ith_cirs = 1:nth_cirs
            x_cirs(iR, ith_cirs) = radii(iR) .* cos(th_cirs(ith_cirs)) + center_x;
            y_cirs(iR, ith_cirs) = radii(iR) .* sin(th_cirs(ith_cirs));
        end
    end
else
    center_x1 = 2. * pi / waveT / kvalue_local(1);           % x-coordinate of circle centers
    n_pos1_mode1 = (2. * pi / waveT - 1.) / kvalue_local(1);
    n_neg1_mode1 = (2. * pi / waveT + 1.) / kvalue_local(1);
    center_x2 = 2. * pi / waveT / kvalue_local(2);           % x-coordinate of circle centers
    n_pos1_mode2 = (2. * pi / waveT - 1.) / kvalue_local(2);
    n_neg1_mode2 = (2. * pi / waveT + 1.) / kvalue_local(2);
end

disp("Calculation finished. Start plotting...");

scrsz = get(0, "ScreenSize");
hLF2 = figure('Position', [0 0 scrsz(3) 0.25*scrsz(3)], "Visible", "off");
t = tiledlayout(1, 4);

% =========================================
ax1 = nexttile; % Plot f_VperpVz at the last time slices
tmpf = squeeze(f_VperpVz);
if plot_logf
    tmpf = log10(tmpf);
    lgmax = max(tmpf,[],'all');
    tmpf(tmpf < lgmax - logdyn) = lgmax - logdyn;
end

[~, h1] = contourf(VZ, VPERP, tmpf(:, :, end-1), 50);
hold on;
annotation('textbox', [0.05, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(1), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);
annotation('textbox', [0.28, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(2), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);
annotation('textbox', [0.54, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(3), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);
annotation('textbox', [0.77, 0.44, 0.5, 0.5], "Interpreter", "latex", "String", anno_labels(4), 'FitBoxToText','on', "EdgeColor","none", "FontSize",28);

if size(kvalue_local,2) == 2
    param_line = sprintf("Field Choice: %d, $\\beta_i = %2.1f, \\mathbf{r} = (%3.2f, %3.2f, %3.2f), RSR = %3.2f, \\delta \\phi = %3.2f \\pi, k_{\\parallel, 1} \\rho_i = %4.3f, k_{\\parallel, 2} \\rho_i = %4.3f, k_{\\perp} \\rho_i = 0.01, t_i = %1.1d T, t_f = (%1.1d, %4.3f T; %4.3f T), (n_{v_x}, n_{v_y}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$", ...
        field_choice, bi, xval, yval, zval, em_eps, delta_phi/pi, kvalue_local(1), kvalue_local(2), t_init/waveT, t_final(1), t_final(end-1)/waveT, dt_final/waveT, nvx, nvy, nvz);
else
    param_line = sprintf("Field Choice: %d, $\\beta_i = %2.1f, \\mathbf{r} = (%3.2f, %3.2f, %3.2f), RSR = %3.2f, \\delta \\phi = %3.2f \\pi, k_{\\parallel} \\rho_i = %4.3f, k_{\\perp} \\rho_i = 0.01, t_i = %1.1d T, t_f = (%1.1d, %4.3f T; %4.3f T), (n_{v_x}, n_{v_y}, n_{v_z}) = (%1.1d, %1.1d, %1.1d)$", ...
        field_choice, bi, xval, yval, zval, em_eps, delta_phi/pi, kvalue_local, t_init/waveT, t_final(1), t_final(end-1)/waveT, dt_final/waveT, nvx, nvy, nvz);
end

if keep_anno
    annotation('textbox', [0.07, 0.48, 0.5, 0.5], "Interpreter", "latex", "String", param_line, 'FitBoxToText','on', "EdgeColor","none", "FontSize",14);
end

set(h1,'edgecolor','none');
colorbar('FontSize',20,'FontName','TimesNewRoman');
if plot_circs
    for iR = 1:nR
        plot(x_cirs(iR, :), y_cirs(iR, :), 'LineStyle','-', 'LineWidth', 1, 'Color','#6aa84f');
    end
end
if size(kvalue_local, 2) == 2
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

format_subplot('$v_\parallel/v_{ti}$', '$v_\perp/v_{ti}$', "$\log_{10}f(v_\parallel, v_\perp)$");
colormap(ax1, plasma);

% =========================================
ax2 = nexttile; % Plot tavg_ceperp_VperpVz
tmpf = squeeze(tavg_ceperp_VperpVz);
[~, h2] = contourf(VZ, VPERP, tmpf, 50);
hold on;

beta_label = sprintf("$\\beta_i = %2.1f$", bi);

text(-3.8, 1.1, beta_label,'Interpreter','latex', ...
    'FontName','TimesNewRoman','FontSize',24, ...
    'FontWeight','bold', 'HorizontalAlignment','left', ...
    'VerticalAlignment', 'top', ...
    'Rotation',90);

if plot_circs
    for iR = 1:nR
        plot(x_cirs(iR, :), y_cirs(iR, :), 'LineStyle','-', 'LineWidth', 1, 'Color','#6aa84f');
    end
end
if size(kvalue_local, 2) == 2
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
ax3 = nexttile; % Plot tavg_cex_VxVy
tmpf = squeeze(tavg_cex_VxVy);
[~, h3] = contourf(VX, VY, transpose(tmpf), 50);
hold on;

set(h3,'edgecolor','none');
colorbar('FontSize', 20,'FontName','TimesNewRoman');

% Sets the colormap limits such that the center of the color bar is 0
cL = caxis;  
caxis(ax3, [-max(abs(cL)) max(abs(cL))]); 

daspect([1 1 1]);

format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_x}(v_x, v_y)$');
colormap(ax3, bluewhitered);
    
% =========================================
ax4 = nexttile; % Plot tavg_cey_VxVy
tmpf = squeeze(tavg_cey_VxVy);
[~, h4] = contourf(VX, VY, transpose(tmpf), 50);
hold on;

set(h4,'edgecolor','none');
colorbar('FontSize', 20,'FontName','TimesNewRoman');

% Sets the colormap limits such that the center of the color bar is 0
cL = caxis;  
caxis(ax4, [-max(abs(cL)) max(abs(cL))]); 

daspect([1 1 1]);
format_subplot('$v_x/v_{ti}$', '$v_y/v_{ti}$', '$C_{E_y}(v_x, v_y)$');

colormap(ax4, bluewhitered);

disp("Plotting finished.");

if save_figure_local
    disp("Start saving figures...");
    if keep_anno
        pngname = sprintf("./png/iSHCDV21_%s_anno.png", time_suffix);
        epsname = sprintf("./eps/iSHCDV21_%s_anno.eps", time_suffix);
        pdfname = sprintf("./pdf/iSHCDV21_%s_anno.pdf", time_suffix);       
    else
        pngname = sprintf("./png/iSHCDV21_%s.png", time_suffix);
        epsname = sprintf("./eps/iSHCDV21_%s.eps", time_suffix);
        pdfname = sprintf("./pdf/iSHCDV21_%s.pdf", time_suffix);
    end
    % Save png and eps
    print(hLF2, pngname, '-dpng', '-r150');  % 150 dpi resolution
    print(hLF2, epsname, '-painters','-depsc','-r150');

    % Save pdf
    % Set paper units to inches (or points, cm, etc.)
    set(hLF2, 'PaperUnits', 'inches');
    % Match the figure size (convert pixels to inches, or set directly)
    fig_width = 18;   % in inches 18
    fig_height = 4.5;   % in inches 4.5
    set(hLF2, 'PaperSize', [fig_width fig_height], 'PaperPosition', [0 0 fig_width fig_height]);
    print(hLF2, pdfname, '-dpdf', '-vector');
    disp("Figures saved.");
else
    disp("Figures not saved.");
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

% Display differences between fieldnames of s1 and s2
function show_params_differences(s1, s2)
    params = fieldnames(s1);
    for i = 1:numel(params)
        param_name = params{i};
        param1 = s1.(param_name);
        param2 = s2.(param_name);
        if ~isequal(param1, param2)
            fprintf("Parameter '%s' differs: \n", param_name);
            disp([' From mat data: ', mat2str(param1)]);
            disp([' Local: ', mat2str(param2)]);
        end
    end
end