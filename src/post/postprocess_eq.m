% src/post/postprocess_eq.m
function postprocess_eq(result_or_file)
% Postprocess equilibrium fit results
%   - plot fits to data
%   - plot spline surfaces (univariate or surface_uv or surface_I1I2)
%   - plot parameter activation from Jacobian J
%
% Usage:
%   postprocess_eq('fit_result_eq.mat')
%   postprocess_eq(result_struct)

    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    if ischar(result_or_file) || isstring(result_or_file)
        tmp = load(result_or_file);
        if isfield(tmp,'result')
            result = tmp.result;
        else
            error('MAT-file does not contain variable "result".');
        end
    else
        result = result_or_file;
    end

    % load components
    J        = result.J;
    model    = result.model;
    prot     = result.prot;
    data     = result.data;
    penalty  = result.penalty;

    % create directory for saving postprocessing results
    model_dir = build_model_dir(model, penalty);
    out_dir   = fullfile('postprocessing', model_dir);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % ----------------------------
    % 1) Plot fits to data
    % ----------------------------
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        D = data.(mode);
        p = prot.(mode);

        Pfit = zeros(size(D.P11));
        for k = 1:numel(D.P11)
            F = p.Fbar(D.lambda(k));

            [a,b]  = model.eval_ab(model, F);
            P = eq_P_from_F(a,b,F);
            Pfit(k) = P(1,1);
        end

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];

        plot(D.lambda, D.P11, 'k.', 'MarkerSize', 11); % data
        plot(D.lambda, Pfit, 'LineWidth',1.8);

        % add annotation for MSE and r^2 value between data and fit
        mse    = mean((D.P11 - Pfit).^2) * 1e6; % convert to kPa^2 for more interpretable units
        r2      = 1 - sum((D.P11 - Pfit).^2) / sum((D.P11 - mean(D.P11)).^2);
        annotation_str = sprintf('MSE = %.2f kPa^2\nR^2 = %.4f', mse, r2);
        annotation('textbox', [0.15, 0.7, 0.2, 0.2], 'String', annotation_str, ...
            'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'none', ...
            'FontSize', 10, 'FontName', 'Times New Roman');

        xlabel('$\lambda$'); ylabel('$P_{11}$ [MPa]');
        grid on;
        fontsize(gca, 12, "Points");
        fontname(gca, 'Times New Roman');   
        exportgraphics(gcf, fullfile(out_dir, sprintf('eq_fit_data_%s.pdf', mode)), 'ContentType', 'vector');
    end

    % ----------------------------
    % 2) Plot splines (univariate or surface_uv or surface_I1I2)
    % ----------------------------
    if strcmp(model.type, 'univariate')
        plot_univariate(model, out_dir);
    elseif strcmp(model.type, 'surface_uv')
        plot_surface_uv(model, data, prot, out_dir);
    elseif strcmp(model.type, 'surface_I1I2')
        plot_surface_I1I2(model, data, prot, out_dir); % same plotting as surface_uv but with different mapping
    else
        warning('Unknown model.type="%s", skipping spline plots.', model.type);
    end

    
    % ---------------------------------
    % 3) Parameter activation from J
    % ---------------------------------

    % count data points
    Ndata = 0;
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        Ndata = Ndata + numel(data.(mode).P11);
    end
    J = J(1:Ndata, :); % only use data Jacobian for activation analysis (exclude penalty rows)
    fprintf('Computing parameter activation from Jacobian J of size %d x %d\n', size(J,1), size(J,2));
    act = parameter_activation_from_J(J, model);

    if strcmp(model.type, 'surface_uv')
        plot_surface_activation_uv(act, model, data, prot, out_dir);
    elseif strcmp(model.type, 'surface_I1I2')
        plot_surface_activation_I1I2(act, model, data, prot, out_dir);
    elseif strcmp(model.type, 'univariate')
        plot_univariate_activation(act, model, data, prot, out_dir);
    end
end


% function to generate postprocessing folder based on model type and number of parameters
function out_dir = build_model_dir(model, penalty)
    switch model.type
        case 'univariate'
            n1 = numel(model.I1_vals);
            if model.use_I2
                n2 = numel(model.I2_vals);
            else                
                n2 = 0;
            end
            out_dir = sprintf('univariate_I1_%d_I2_%d', n1, n2);
        case 'surface_uv'
            nu = numel(model.u_sites);
            nv = numel(model.v_sites);
            lambda_pen = penalty.lambda_pen;
            out_dir = sprintf('surface_uv_%dx%d_lambda_%.0e', nu, nv, lambda_pen);
        case 'surface_I1I2'
            nu = numel(model.I1_sites);
            nv = numel(model.I2_sites);
            lambda_pen = penalty.lambda_pen;
            out_dir = sprintf('surface_i1i2_%dx%d_lambda_%.0e', nu, nv, lambda_pen);
        otherwise
            error('Unknown model type: %s', model.type);
    end 
end

% ---------------------------- helpers ----------------------------
function plot_univariate(model, out_dir)

    lw = 2;
    ms_site     = 40; % scatter size
    col.spline  = [0.00 0.30 0.60];   % muted blue
    col.site    = [0.55 0.10 0.10];   % muted dark red
    figsize     = [8,6];
    fs          = 12;

    y_lim = [0, 15];

    I1_plot     = linspace(min(model.I1_knots), max(model.I1_knots), 100);
    W_I1_plot   = fnval(model.WI1_spline, I1_plot);
    W_I1_ctrl   = fnval(model.WI1_spline, model.I1_knots);

    fig = figure('Color','w'); hold on;
    fig.Units = 'centimeters';
    fig.Position(3:4) = figsize;


    plot(I1_plot, W_I1_plot, '-', 'LineWidth',lw, 'Color', col.spline);

    % interpolation values
    scatter(model.I1_knots, W_I1_ctrl, ms_site, ...
        'Marker', 'o', ...
        'MarkerFaceColor', col.site, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.95);

    xlabel('$\bar{I}_1$'); ylabel('$W(\bar{I}_1)$'); grid on;
    ylim(y_lim);
    fontsize(gca, fs, "Points");
    fontname(gca, 'Times New Roman');
    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_WI1.pdf'), 'ContentType', 'vector');

    if model.use_I2
        I2_plot     = linspace(min(model.I2_knots), max(model.I2_knots), 100);
        W_I2_plot   = fnval(model.WI2_spline, I2_plot);
        W_I2_ctrl   = fnval(model.WI2_spline, model.I2_knots);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = figsize;

        plot(I2_plot, W_I2_plot, '-', 'LineWidth',lw, 'Color', col.spline);
        scatter(model.I2_knots, W_I2_ctrl, ms_site, ...
            'Marker', 'o', ...
            'MarkerFaceColor', col.site, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.95);

        xlabel('$\bar{I}_2$'); ylabel('$W(\bar{I}_2)$'); grid on;
        ylim(y_lim);
        fontsize(gca, fs, "Points");
        fontname(gca, 'Times New Roman');
        exportgraphics(gcf, fullfile(out_dir, 'eq_fit_WI2.pdf'), 'ContentType', 'vector');
    end
end


function plot_surface_uv(model, data, prot, out_dir)
% Plots for a surface model W(u,v):
%   1) parameter domain (u,v): show optimized interpolation values
%   2) physical domain (I1,I2): show surface + data only

    % ---------- style ----------
    ms_data = 25;   % scatter size
    ms_site = 14;
    alpha_surf = 0.85;

    col.site  = [0.55 0.10 0.10];   % muted dark red
    col.data  = [1 1 1];            % white for contrast against surface
    col.bound = [0.20 0.20 0.20];   % dark gray

    lw_bound = 0.9;

    fs = 12;
    figsize = [8,6];

    % ---------- precompute site values in (u,v) ----------
    [UU, VV]    = ndgrid(model.u_sites, model.v_sites);
    W_sites_uv  = fnval(model.Wsurf, [UU(:).'; VV(:).']);

    %% =========================================================
    %% Panel 1: W(u,v)
    %% =========================================================
    fig = figure('Color','w','Renderer','painters');
    fig.Units = 'centimeters';
    fig.Position(3:4) = figsize;
    hold on;

    % --- surface ---
    nPlot = 40; % number of points in each direction for plotting surface
    u_grid     = linspace(min(model.u_sites), max(model.u_sites), nPlot);
    v_grid     = linspace(min(model.v_sites), max(model.v_sites), nPlot);
    [u_grid, v_grid]  = ndgrid(u_grid, v_grid);
    W_grid = fnval(model.Wsurf, [u_grid(:).'; v_grid(:).']);
    W_grid = reshape(W_grid, size(u_grid));
    surf(u_grid, v_grid, W_grid, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', alpha_surf);
    colormap(get_parula_colormap(256));

    % --- interpolation values ---
    scatter3(UU(:), VV(:), W_sites_uv(:), ms_site, ...
        'Marker', 'o', ...
        'MarkerFaceColor', col.site, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.95); 

    xlabel('$\xi$', 'Interpreter', 'latex');
    ylabel('$\eta$', 'Interpreter', 'latex');
    zlabel('$W(\xi, \eta)$', 'Interpreter', 'latex');

    grid off;
    box on;
    view(-35, 22);
    camproj('orthographic');
    axis tight;

    set(gca, ...
        'FontSize', 11, ...
        'FontName', 'Times New Roman', ...
        'LineWidth', 0.8, ...
        'TickDir', 'out');

    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_surface_uv.pdf'), 'ContentType', 'vector');
    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_surface_uv.png'), 'Resolution', 1200);


    %% =========================================================
    %% Panel 2: W(I1,I2)
    %% =========================================================
    fig = figure('Color','w','Renderer','painters');
    fig.Units = 'centimeters';
    fig.Position(3:4) = figsize;
    hold on;

    npoints = 40;
    u_plot = linspace(min(model.u_sites), max(model.u_sites), npoints);
    v_plot = linspace(min(model.v_sites), max(model.v_sites), npoints);
    [u_grid, v_grid] = ndgrid(u_plot, v_plot);

    I1_grid = u_grid * (model.I1_max - 3) + 3; % inverse mapping of u to I1

    flow = zeros(size(u_plot));
    fup  = zeros(size(u_plot));
    for j = 1:numel(u_plot)
        I1 = u_plot(j) * (model.I1_max - 3) + 3; % inverse mapping of u to I1
        [flow(j), fup(j)] = invariant_bounds_I2_from_I1(I1);
    end

    if model.use_polyconvex_I2
        flow = flow.^(3/2) - 3.0 * sqrt(3);
        fup  = fup.^(3/2)  - 3.0 * sqrt(3);
    end

    FLOWg = repmat(flow(:), 1, numel(v_plot));
    FUPg  = repmat(fup(:),  1, numel(v_plot));
    I2_grid = FLOWg + v_grid .* (FUPg - FLOWg);

    W_grid = fnval(model.Wsurf, [u_grid(:).'; v_grid(:).']);
    W_grid = reshape(W_grid, size(u_grid));

    surf(I1_grid, I2_grid, W_grid, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', alpha_surf);
    colormap(get_parula_colormap(256));

 
    % admissible domain bounds
    I1_bound = u_plot * (model.I1_max - 3) + 3; % inverse mapping of u to I1
    flow_bound = flow;
    fup_bound  = fup;

    % floor height for drawing admissible region
    z0 = min(W_grid(:));

     % shaded admissible domain
    fill3([I1_bound, fliplr(I1_bound)], ...
          [flow_bound, fliplr(fup_bound)], ...
          z0 * ones(1,2*numel(I1_bound)), ...
          [0.88 0.88 0.88], ...
          'EdgeColor','none', ...
          'FaceAlpha',0.55);

    % admissible boundary curves
    plot3(I1_bound, flow_bound, z0 * ones(size(I1_bound)), ...
        '--', 'Color',[0.15 0.15 0.15], 'LineWidth',lw_bound);

    plot3(I1_bound, fup_bound,  z0 * ones(size(I1_bound)), ...
        '--', 'Color',[0.15 0.15 0.15], 'LineWidth',lw_bound);  

    xlabel('$\bar{I}_1$', 'Interpreter', 'latex');
    ylabel('$\bar{I}_2$', 'Interpreter', 'latex');
    zlabel('$W(\bar{I}_1, \bar{I}_2)$', 'Interpreter', 'latex');

    grid off;
    box on;
    view(-35, 22);
    camproj('orthographic');
    axis tight;

    set(gca, ...
        'FontSize', 11, ...
        'FontName', 'Times New Roman', ...
        'LineWidth', 0.8, ...
        'TickDir', 'out');

    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_surface_uv_physical.pdf'), 'ContentType', 'vector');
    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_surface_uv_physical.png'), 'Resolution', 1200);
end



function cmp = get_parula_colormap(n)
    % Get parula colormap with n colors
    cmp = parula(n);
end



function plot_surface_I1I2(model, data, prot, out_dir)
% Plots for a surface model W(I1,I2):
%     1) parameter domain (I1,I2): show optimized interpolation values

    % ---------- style ----------
    ms_data = 25;   % scatter size, not radius
    ms_site = 14;
    alpha_surf = 0.85;

    col.site  = [0.55 0.10 0.10];   % muted dark red
    col.data  = [1 1 1];   % white for contrast against surface

    % ---------- figure ----------
    fig = figure('Color','w', 'Renderer','painters');
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8,6];
    hold on;

    % ---------- surface ----------
    nPlot = 40; % number of points in each direction for plotting surface
    I1_grid     = linspace(min(model.I1_sites), max(model.I1_sites), nPlot);
    I2_grid     = linspace(min(model.I2_sites), max(model.I2_sites), nPlot);
    [I1g, I2g]  = ndgrid(I1_grid, I2_grid);

    W_grid = fnval(model.Wsurf, [I1g(:).'; I2g(:).']);
    W_grid = reshape(W_grid, size(I1g));
    surf(I1g, I2g, W_grid, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', alpha_surf);
    colormap(get_parula_colormap(256));

    % --- interpolation values ---
    [I1g, I2g] = ndgrid(model.I1_sites, model.I2_sites);
    W_sites = fnval(model.Wsurf, [I1g(:).'; I2g(:).']);
    scatter3(I1g(:), I2g(:), W_sites(:), ms_site, ...
        'Marker', 'o', ...
        'MarkerFaceColor', col.site, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.95); 

    % ---------- axes ----------
    xlabel('$\bar{I}_1$', 'Interpreter', 'latex');
    ylabel('$\bar{I}_2$', 'Interpreter', 'latex');
    zlabel('$W(\bar{I}_1, \bar{I}_2)$', 'Interpreter', 'latex');

    grid off;
    box on;
    view(-35, 22);
    camproj('orthographic');
    axis tight;

    set(gca, ...
        'FontSize', 11, ...
        'FontName', 'Times New Roman', ...
        'LineWidth', 0.8, ...
        'TickDir', 'out');

    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_surface_I1I2.pdf'), 'ContentType', 'vector');
    exportgraphics(gcf, fullfile(out_dir, 'eq_fit_surface_I1I2.png'), 'Resolution', 1200);
end



function plot_sampled_invariants(model, data, prot, ms_data, color, map_to_uv_space, z_offset)

    mk.UT = 'o';
    mk.ET = 'square';
    mk.PS = 'diamond';
    
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        D = data.(mode);
        p = prot.(mode);

        switch mode
            case 'UT', m = mk.UT;
            case 'ET', m = mk.ET;
            case 'PS', m = mk.PS;
            otherwise, m = 'o';
        end

        I1_data = zeros(numel(D.lambda),1);
        I2_data = zeros(numel(D.lambda),1);
        W_data  = zeros(numel(D.lambda),1);

        for k = 1:numel(D.lambda)
            F = p.Fbar(D.lambda(k));
            C = F.' * F;

            I1 = trace(C);
            I2 = 0.5 * (I1^2 - trace(C * C));

            if map_to_uv_space
                [u, v] = map_I1I2_to_uv(I1, I2, model.delta, model.I1_max, model.use_polyconvex_I2);
                u = min(max(u, model.u_sites(1)), model.u_sites(end));
                v = min(max(v, model.v_sites(1)), model.v_sites(end));

                I1_data(k) = u;
                I2_data(k) = v;

                W_data(k)  = fnval(model.Wsurf, [u; v]);
            else
                I1_data(k) = I1;
                I2_data(k) = I2;
                if model.use_polyconvex_I2
                    I2_data(k) = I2_data(k)^(3/2) - 3.0 * sqrt(3);
                end

                W_data(k)  = fnval(model.Wsurf, [I1; I2]);
            end

        end

        if isempty(z_offset)
            height = 0;
        else
            height = z_offset + W_data;
        end

        scatter3(I1_data, I2_data, height, ms_data, ...
            'Marker', m, ...
            'MarkerFaceColor', color, ...
            'MarkerEdgeColor', [0.15 0.15 0.15], ...
            'MarkerFaceAlpha', 0.90);
    end
end



function act = parameter_activation_from_J(J, model)
    % Compute parameter activation from Jacobian J of size (Ndata x nParams)


    raw = vecnorm(J, 2, 1).'; % nParams x 1 vector of raw activations (L2 norm of each column)
    rel = raw / max(raw);

    act = struct();
    act.raw = raw;
    act.rel = rel;
    act.logrel = log10(rel + 1e-16);

    if strcmp(model.type, 'surface_uv')
        nu = numel(model.u_sites);
        nv = numel(model.v_sites);
        act.grid_raw = reshape(raw, [nu, nv]);
        act.grid_rel = reshape(rel, [nu, nv]);
        act.grid_logrel = reshape(act.logrel, [nu, nv]);
    elseif strcmp(model.type, 'surface_I1I2')
        nu = numel(model.I1_sites);
        nv = numel(model.I2_sites);
        act.grid_raw = reshape(raw, [nu, nv]);
        act.grid_rel = reshape(rel, [nu, nv]);
        act.grid_logrel = reshape(act.logrel, [nu, nv]);
    elseif strcmp(model.type, 'univariate')
        n1 = numel(model.I1_vals);
        act.I1_raw = raw(1:n1);
        act.I1_rel = rel(1:n1);
        act.I1_logrel = act.logrel(1:n1);

        if model.use_I2
            n2 = numel(model.I2_vals);
            act.I2_raw = raw(n1+1:n1+n2);
            act.I2_rel = rel(n1+1:n1+n2);
            act.I2_logrel = act.logrel(n1+1:n1+n2);
        end
    end
end


function plot_surface_activation_uv(act, model, data, prot, out_dir)

    fig = figure('Color','w');
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 6];
    hold on;

    imagesc(model.u_sites, model.v_sites, act.grid_logrel.');
    set(gca, 'YDir', 'normal');
    colormap(get_parula_colormap(256));

    [UU, VV] = meshgrid(model.u_sites, model.v_sites);
    plot(UU(:), VV(:), 'o', ...
        'MarkerSize', 2.5, ...
        'MarkerFaceColor', [0.1 0.1 0.1], ...
        'MarkerEdgeColor', 'none');
    colorbar;

    plot_sampled_invariants(model, data, prot, 25, [1 1 1], true, []);

    xlabel('$\xi$','Interpreter','latex');
    ylabel('$\eta$','Interpreter','latex');
    fontsize(gca, 12, "Points");
    fontname(gca, 'Times New Roman');
    
    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_uv.pdf'), 'ContentType', 'vector');
    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_uv.png'), 'Resolution', 1200);
end


function plot_surface_activation_physical(act, model, out_dir)
    fig = figure('Color','w');
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 6];

    [UU, VV] = ndgrid(model.u_sites, model.v_sites);
    I2_sites = zeros(size(UU));

    for i = 1:numel(UU)
        u = UU(i) * (model.I1_max - 3) + 3; % inverse mapping of u to I1
        v = VV(i);
        [flow, fup] = invariant_bounds_I2_from_I1(u);
        I2_sites(i) = flow + v * (fup - flow);
    end
    if model.use_polyconvex_I2
        I2_sites = I2_sites.^(3/2) - 3.0 * sqrt(3);
    end

    fig = figure('Color','w');
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 6];

    scatter(UU(:), I2_sites(:), 120, act.grid_logrel(:), 'filled');
    colorbar;
    xlabel('$\bar{I}_1$','Interpreter','latex');
    ylabel('$\bar{I}_2$','Interpreter','latex');
    grid on;

    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_physical_I1I2.pdf'), 'ContentType', 'vector');
    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_physical_I1I2.png'), 'Resolution', 1200);
end



function plot_surface_activation_I1I2(act, model, data, prot, out_dir)
    fig = figure('Color','w');
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 6];
    hold on;

    imagesc(model.I1_sites, model.I2_sites, act.grid_logrel.');
    set(gca, 'YDir', 'normal');
    colormap(get_parula_colormap(256));

    [I1g, I2g] = meshgrid(model.I1_sites, model.I2_sites);
    plot(I1g(:), I2g(:), 'o', ...
        'MarkerSize', 2.5, ...
        'MarkerFaceColor', [0.1 0.1 0.1], ...
        'MarkerEdgeColor', 'none');
    colorbar;

    plot_sampled_invariants(model, data, prot, 25, [1 1 1], false, []);  

    xlabel('$\bar{I}_1$','Interpreter','latex');
    ylabel('$\bar{I}_2$','Interpreter','latex');
    fontsize(gca, 12, "Points");
    fontname(gca, 'Times New Roman');

    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_I1I2.pdf'), 'ContentType', 'vector');
    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_I1I2.png'), 'Resolution', 1200);
end


function plot_univariate_activation(act, model, data, prot, out_dir)

    % Compute invariants 
    I1 = [];
    I2 = [];
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        D = data.(mode);
        p = prot.(mode);
        for k = 1:numel(D.lambda)
            F = p.Fbar(D.lambda(k));
            C = F.' * F;

            I1(end+1) = trace(C);
            I2(end+1) = 0.5 * (I1(end)^2 - trace(C * C));
            if model.use_polyconvex_I2
                I2(end) = I2(end)^(3/2) - 3.0 * sqrt(3);
            end
        end
    end


    fs = 12;
    col.spline = [0.00 0.30 0.60];
    col.rug    = [0.2 0.2 0.2];

    fig = figure('Color','w');
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 6];

    ax = axes(fig);
    hold(ax, 'on');

    % Main activation plot
    plot(ax, model.I1_knots, act.I1_logrel, '-o', ...
        'LineWidth', 1.5, ...
        'Color', col.spline, ...
        'MarkerFaceColor', col.spline);

    % Rug plot of sampled I1 values (all modes combined)
    rugplot(ax, I1, ...
        'Height', 0.06, ...
        'Color', [0.3 0.3 0.3], ...
        'LineWidth', 0.5);

    xlabel('$\bar{I}_1$','Interpreter','latex');
    grid on;
    fontsize(gca, fs, "Points");
    fontname(gca, 'Times New Roman');
    exportgraphics(gcf, fullfile(out_dir, 'eq_activation_I1.pdf'), 'ContentType', 'vector');

    if model.use_I2
        fig = figure('Color','w');
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];

        ax = axes(fig);
        hold(ax, 'on');

        plot(ax, model.I2_knots, act.I2_logrel, '-o', ...
            'LineWidth', 1.5, ...
            'Color', col.spline, ...
            'MarkerFaceColor', col.spline);

        y_lim = ax.YLim;
        y_lim(2) = 0; y_lim(1) = -1.5;
        ax.YLim = y_lim;

        rugplot(ax, I2, ...
            'Height', 0.06, ...
            'Color', [0.3 0.3 0.3], ...
            'LineWidth', 0.5);

        
        xlabel('$\bar{I}_2$','Interpreter','latex');
        grid on;
        fontsize(gca, fs, "Points");
        fontname(gca, 'Times New Roman');
        exportgraphics(gcf, fullfile(out_dir, 'eq_activation_I2.pdf'), 'ContentType', 'vector');
    end
end




function rugplot(ax, x, varargin)
%RUGPLOT  Draw a rug plot on axis ax for data x

    % --- defaults ---
    opts.Height    = 0.05;         % relative height (fraction of axis)
    opts.Color     = [0.2 0.2 0.2];
    opts.LineWidth = 0.5;
    opts.Location  = 'bottom';     % 'bottom' or 'top'

    % --- parse optional args ---
    for i = 1:2:numel(varargin)
        opts.(varargin{i}) = varargin{i+1};
    end

    x = x(:);  % ensure column

    yl = ylim(ax);
    dy = yl(2) - yl(1);
    h  = opts.Height * dy;

    switch lower(opts.Location)
        case 'bottom'
            y0 = yl(1);
            y1 = y0 + h;
        case 'top'
            y1 = yl(2);
            y0 = y1 - h;
        otherwise
            error('Unknown Location option.');
    end

    hold(ax, 'on');

    for i = 1:numel(x)
        line(ax, [x(i), x(i)], [y0, y1], ...
            'Color', opts.Color, ...
            'LineWidth', opts.LineWidth);
    end

    % restore original y-limits
    ylim(ax, yl);
end