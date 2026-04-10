% main/run_fit_eq.m
function run_fit_eq(opts)
    % opts.type = 'univariate' | 'surface_uv' | 'surface_i1i2'

    % default model: univariate spline with W(I1) + W(I2)
    if nargin < 1 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'type') || isempty(opts.type), opts.type = 'univariate'; end
    
    % add src/ to path 
    addpath(genpath('..'));

    % ----------------------------
    % 1) Load data
    % ----------------------------
    type = 'treloar';
    modes = {'UT', 'ET', 'PS'}; % fit all modes jointly

    switch type
        case 'treloar'
            data = get_data(type, [], modes);
        otherwise
            error('Unknown data type');
    end

    fprintf('Loaded data type: %s\n', type);
    fprintf('Modes: '); fprintf('%s ', modes{:}); fprintf('\n');
    for mode = modes
        mode = mode{1};
        D = data.(mode);
        fprintf(' Mode: %s, number of data points: %d\n', mode, numel(D.P11));
    end

    % ----------------------------
    % 2) Model + protocol handles
    % ----------------------------
    prot = protocol_map();

    I1_max = -inf; I2_max = -inf;
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        p = prot.(mode);
        lam_max = max(data.(mode).lambda);
        I1_max = max(I1_max, p.I1_max(lam_max));
        I2_max = max(I2_max, p.I2_max(lam_max));
    end

    c = 1.001; % small safety factor for interpolation range of splines
    model = eq_params_init_wrapper(I1_max * c, I2_max * c, opts);

    % ----------------------------
    % 3) Pack parameters for solver
    % ----------------------------
    x0 = model.pack(model);
    lb = model.lb();
    ub = model.ub();

    fprintf('Number of parameters to fit: %d\n', numel(x0));
    fprintf('  lb: '); fprintf('%.3f ', lb); fprintf('\n');
    fprintf('  ub: '); fprintf('%.3f ', ub); fprintf('\n');
    fprintf('  x0: '); fprintf('%.3f ', x0); fprintf('\n');

    [Aineq, bineq] = constraints_hyper(model);
    fprintf('Number of linear constraints: %d\n', size(Aineq,1));

    % ----------------------------
    % 4) Build constant operators once
    % ----------------------------
    t_start = tic();
    [A_data, y_data] = build_linear_system(model, prot, data);
    fprintf('Linear data system built in %.2f seconds.\n', toc(t_start));

    penalty = [];
    if strcmpi(model.type,'surface_uv') || strcmpi(model.type,'surface_i1i2')
        penalty = build_second_derivative_penalty_operator(model, 1.0); % unscaled operator
        fprintf('Penalty operator built with %d rows.\n', size(penalty.J,1));
    end

    % ----------------------------
    % 5) Choose algorithm
    % ----------------------------
    algorithm = 'lsqlin';   % 'lsqlin' | 'lsqnonlin' | 'pareto'

    switch algorithm
        case 'lsqlin'
            lambda_pen = choose_default_lambda(model);

            opts = optimoptions('lsqlin', ...
                'Algorithm','interior-point', ...
                'Display','iter', ...
                'OptimalityTolerance',1e-12, ...
                'FunctionTolerance',1e-12, ...
                'StepTolerance',1e-12);

            [A_aug, y_aug] = assemble_augmented_system(A_data, y_data, penalty, lambda_pen);

            [rvec0, ~]  = objective_multi(x0, model, prot, data, scaled_penalty(penalty, lambda_pen));
            r_lin       = A_aug*x0 - y_aug;
            fprintf('Checking consistency of linear system at initial guess:\n');
            fprintf('  ||(A*x0 - y) - rvec0|| = %.6e\n', norm(r_lin - rvec0));

            [x_opt, resnorm, residual, exitflag, output, lambda] = ...
                lsqlin(A_aug, y_aug, Aineq, bineq, [], [], lb, ub, x0, opts); %#ok<ASGLU>

            [r,J] = objective_multi(x_opt, model, prot, data, scaled_penalty(penalty, lambda_pen));

            fprintf('Final resnorm: %.6e\n', resnorm);
            fprintf('Final sum of squares: %.6e\n', sum(residual.^2));
            fprintf('Sum of squares from objective_multi: %.6e\n', sum(r.^2));

            % unpack final model for later analysis
            model = model.unpack(x_opt, model);

            % save results to a .mat file for later analysis
            result = struct();
            result.x_opt        = x_opt;
            result.residual     = residual;
            result.J            = J;
            result.model        = model;
            result.prot         = prot;
            result.data         = data;
            result.penalty      = scaled_penalty(penalty, lambda_pen);
            result.lambda_pen   = lambda_pen;
            save('fit_result_eq.mat', 'result');

        case 'lsqnonlin'

            % optimization options
            tol = 1e-10;
            opts = optimoptions('lsqnonlin', ...
                'Display','iter', ...
                'OptimalityTolerance',tol, ...
                'FunctionTolerance',tol, ...
                'StepTolerance',tol, ...
                'ConstraintTolerance',tol, ...
                'SpecifyObjectiveGradient',true);

            % objective function handle for lsqnonlin
            obj = @(x) objective_multi(x, model, prot, data, scaled_penalty(penalty, choose_default_lambda(model)));

            % check the provided gradients against finite differences at the initial guess
            fd_opts         = optimoptions('lsqnonlin', 'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-6);
            [valid, error]  = checkGradients(obj, x0, fd_opts, Display='on');

            % Run optimization
            t_start = tic();
            [x_opt,resnorm,residual,~,~,~,J] = lsqnonlin(obj, x0, lb, ub, Aineq, bineq, [], [], [], opts);
            t_end = toc(t_start);
            fprintf('Optimization completed in %.2f seconds.\n', t_end);

            % print resnorm and sum of squares
            fprintf('Final resnorm: %.6e\n', resnorm);
            fprintf('Final sum of squares: %.6e\n', sum(residual.^2)); 

            model = model.unpack(x_opt, model);

            % save results to a .mat file for later analysis
            result = struct();
            result.x_opt        = x_opt;
            result.residual     = residual;
            result.J            = J;
            result.model        = model;
            result.prot         = prot;
            result.data         = data;
            result.penalty      = penalty;
            save('fit_result_eq.mat', 'result');

        case 'pareto'
            
            lambda_list = logspace(-10, 2, 80); %  pareto sweep over these values of lambda_pen

            % optimization options for lsqlin in the pareto sweep
            opts = optimoptions('lsqlin', ...
                'Algorithm','interior-point', ...
                'Display','off', ...
                'OptimalityTolerance',1e-12, ...
                'FunctionTolerance',1e-12, ...
                'StepTolerance',1e-12);

            sweep = run_pareto_sweep(A_data, y_data, penalty, ...
                                     Aineq, bineq, lb, ub, x0, ...
                                     lambda_list, opts, model, prot, data);
                        
            save('lambda_sweep_eq.mat', 'sweep');

            fprintf('\nPareto sweep completed.\n');
            fprintf('Selected lambda_pen = %.6e\n', sweep.lambda_star);
            fprintf('Selected phi_data   = %.6e\n', sweep.phi_data(sweep.idx_corner));
            fprintf('Selected phi_pen    = %.6e\n', sweep.phi_pen(sweep.idx_corner));

            % unpack the selected solution and evaluate the objective for later analysis
            x_opt = sweep.x_all(:, sweep.idx_corner);
            lambda_pen = sweep.lambda_star;

            model = model.unpack(x_opt, model);
            [r,J] = objective_multi(x_opt, model, prot, data, scaled_penalty(penalty, lambda_pen));

            % save results to a .mat file for later analysis
            result = struct();
            result.x_opt        = x_opt;
            result.residual     = r;
            result.J            = J;
            result.model        = model;
            result.prot         = prot;
            result.data         = data;
            result.penalty      = scaled_penalty(penalty, lambda_pen);
            save('fit_result_eq.mat', 'result');

        otherwise
            error('Unknown algorithm.');
    end
end

% -----------------------------------------------------------------
% ------------------------ local functions ------------------------
% -----------------------------------------------------------------

function [rvec, J] = objective_multi(x, model, prot, data, penalty)
    % Objective function for nonlinear least squares optimization (lsqnonlin):
    %
    %   minimizes  ||rvec||^2 = sum_over_modes w_mode*(1/Nm)*sum_k (P11_pred - P11_data)^2 + lambda_pen * ||penalty.J * x||^2
    %   the jacbian J is the sensitivity of the residuals to the parameters x, which is needed for gradient-based optimization.

    model = model.unpack(x, model); % unpack the current parameter vector "x" into the model strucutre

    % Count total residuals across all modes
    Ntot = 0;
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        Ntot = Ntot + numel(data.(mode).P11);
    end

    % Initialize residual vector and Jacobian matrix
    rvec = zeros(Ntot, 1);
    if nargout > 1
        J = zeros(Ntot, numel(x));
    end

    % Loop over modes 
    row = 0;
    for i = 1:numel(data.modes) % UT, ET, PS
        mode = data.modes{i};
        D = data.(mode);
        p = prot.(mode);

        w_mode = 1.0;
        if isfield(D,'w') && ~isempty(D.w)
            w_mode = D.w;
        end

        % Scale residuals by sqrt(w_mode / Nm) 
        Nm      = numel(D.P11); % Number of data points for this mode
        scale   = sqrt(w_mode / Nm);

        % Loop over data points in this mode
        for k = 1:Nm 
            row = row + 1;

            % kinematics for this data point
            lam = D.lambda(k); 
            F   = p.Fbar(lam);

            % Piola prediction for this data point
            [a,b]  = model.eval_ab(model, F);
            P_pred = eq_P_from_F(a,b,F); % 3x3 matrix

            rvec(row) = scale * (P_pred(1,1) - D.P11(k)); % residual for P11

            if nargout > 1
                [da_dx, db_dx] = model.eval_dabdx(model, F);
                dP_dx = eq_dPdx_from_F(da_dx, db_dx, F);

                % residual is r = scale * (P11_pred - P11_data) --> dr/dx = scale * dP11_dx 
                J(row, :) = scale * reshape(dP_dx(1,:), 1, []);
            end
        end
    end

    % concatenate penalty residuals and jacobian if provided
    if isfield(penalty,'J') && ~isempty(penalty.J)
        r_pen = penalty.J * x; % nPenalty x 1
        rvec  = [rvec; r_pen];
        if nargout > 1
            J = [J; penalty.J]; % (Ntot+nPenalty) x nParams
        end
    end
end



function lambda_pen = choose_default_lambda(model)
    % default values for lambda_pen: 
    % 
    %       penalty: lambda_pen * ||J*x||^2 where J is the second derivative penalty operator

    if strcmpi(model.type,'surface_uv')
        lambda_pen = 5.917686e-05;
    elseif strcmpi(model.type,'surface_i1i2')
        lambda_pen = 9.712740e-04;
    else
        lambda_pen = 0.0; % no penalty for univariate splines
    end
end



function penalty = build_second_derivative_penalty_operator(model, lambda_pen)
% Penalty on directional second derivatives of surface spline:
%   lambda * \int (W_uu^2 + W_vv^2) dOmega
%
% Returned as residual vector for lsqnonlin so that ||rvec||^2
% approximates the integral.

    if nargin < 2
        lambda_pen = 1.0;
    end

    if lambda_pen <= 0
        penalty.J = [];
        penalty.lambda = 0.0;
        return;
    end

    % derivative sensitivity splines
    Buu = fnder(model.Wsurf_dtheta, [2 0]);
    Bvv = fnder(model.Wsurf_dtheta, [0 2]);

    % unique knot lines define spline elements
    ku = unique(model.Wsurf.knots{1});
    kv = unique(model.Wsurf.knots{2});

    % remove zero-width knot spans
    ku = ku(:).';
    kv = kv(:).';

    % domain bounds for quadrature
    umin = model.Wsurf.knots{1}(1);
    umax = model.Wsurf.knots{1}(end);
    vmin = model.Wsurf.knots{2}(1);
    vmax = model.Wsurf.knots{2}(end);

    Lu = umax - umin;
    Lv = vmax - vmin;

    % 4-point Gauss on [-1,1]
    xi = [-0.861136311594053, -0.339981043584856, ...
           0.339981043584856,  0.861136311594053];
    wi = [ 0.347854845137454,  0.652145154862546, ...
           0.652145154862546,  0.347854845137454];


    J = {}; 

    for i = 1:numel(ku)-1
        ua = ku(i); ub = ku(i+1);
        if ub <= ua, continue; end

        for j = 1:numel(kv)-1
            va = kv(j); vb = kv(j+1);
            if vb <= va, continue; end

            % rectangle mapping from [-1,1]^2
            Ju = 0.5*(ub - ua);
            Jv = 0.5*(vb - va);
            detJ = Ju * Jv;

            for a = 1:numel(xi)
                uq = Ju*xi(a) + 0.5*(ua + ub);
                wu = wi(a);

                for b = 1:numel(xi)
                    vq = Jv*xi(b) + 0.5*(va + vb);
                    wq = wi(b);

                    % sensitivities of second derivatives
                    phi_uu = fnval(Buu, [uq; vq]);   % nParams x 1
                    phi_vv = fnval(Bvv, [uq; vq]);   % nParams x 1

                    % residual for Wuu and Wvv
                    wscale_uu = sqrt(lambda_pen * wu * wq * detJ * (Lu^3 / Lv));
                    wscale_vv = sqrt(lambda_pen * wu * wq * detJ * (Lv^3 / Lu));

                    J{end+1} = wscale_uu * phi_uu(:).';
                    J{end+1} = wscale_vv * phi_vv(:).';
                    
                end
            end
        end
    end

    penalty.J = vertcat(J{:}); % nPenalty x nParams
    penalty.lambda_pen = lambda_pen;
end



function [A,y] = build_linear_system(model, prot, data)
    % Build the linear system A*x = y for the data residuals, where x is the parameter vector.
    % This is used for the 'lsqlin' algorithm where we solve a linear least squares problem with linear constraints. 


    % Count total residuals across all modes
    Ntot = 0;
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        Ntot = Ntot + numel(data.(mode).P11);
    end
    nParams = numel(model.pack(model));

    A = zeros(Ntot, nParams);
    y = zeros(Ntot, 1);

    row = 0;

    % Loop over modes
    for ii = 1:numel(data.modes) % UT, ET, PS
        mode = data.modes{ii};
        D = data.(mode);
        p = prot.(mode);

        w_mode = 1.0;
        if isfield(D,'w') && ~isempty(D.w)
            w_mode = D.w;
        end

        % Scale residuals by sqrt(w_mode / Nm) to ensure balanced contributions from each mode in the least squares objective
        Nm      = numel(D.P11); % Number of data points for this mode
        scale   = sqrt(w_mode / Nm);

        % Loop over data points in this mode
        for k = 1:Nm 

            % kinematics for this data point
            lam = D.lambda(k); 
            F   = p.Fbar(lam);
            C   = F.' * F;
            I1  = trace(C);
            I2  = 0.5*(I1^2 - trace(C*C));

            row = row + 1;
            rowA = zeros(1, nParams); % this will be the row of A corresponding to this data point
            col = 0;

            if model.type == "univariate"
                % univariate splines: W(I1) + W(I2)
    
                [phiI1,~] = model.eval_dabdx(model, F);
                for j = 1:numel(model.I1_vals)
                    col = col + 1;
                    Pj = eq_P_from_F(phiI1(j), 0.0, F); % sensitivity of P11 to this parameter
                    rowA(col) = scale * Pj(1,1);
                end
                
                if model.use_I2
                    sIdx  = numel(model.I1_vals);
                    [~, phiI2] = model.eval_dabdx(model, F);
                    for j = 1:numel(model.I2_vals)
                        col = col + 1;
                        Pj = eq_P_from_F(0.0, phiI2(sIdx+j), F); % sensitivity of P11 to this parameter
                        rowA(col) = scale * Pj(1,1);
                    end
                end
            % surface spline: W(I1,I2) with coupled parameters
            elseif model.type == "surface_I1I2"
                [phiI1, phiI2] = model.eval_dabdx(model, F);
                for j = 1:numel(model.Wvals(:))
                    col = col + 1;
                    Pj = eq_P_from_F(phiI1(j), phiI2(j), F); % sensitivity of P11 to this parameter
                    rowA(col) = scale * Pj(1,1);
                end
            elseif model.type == "surface_uv"
                [da_dx, db_dx] = model.eval_dabdx(model, F);
                for j = 1:numel(model.Wvals(:))
                    col = col + 1;
                    Pj = eq_P_from_F(da_dx(j), db_dx(j), F); % sensitivity of P11 to this parameter
                    rowA(col) = scale * Pj(1,1);
                end
            else
                error('Unknown model type');
            end


            A(row, :) = rowA;
            y(row) = scale * D.P11(k); % this is the target value for P11, scaled by sqrt(w_mode / Nm)
        end
    end
end



function [A_aug, y_aug] = assemble_augmented_system(A_data, y_data, penalty, lambda_pen)

    if isempty(penalty) || ~isfield(penalty,'J') || isempty(penalty.J) || lambda_pen <= 0
        A_aug = A_data;
        y_aug = y_data;
    else
        A_aug = [A_data; sqrt(lambda_pen) * penalty.J];
        y_aug = [y_data; zeros(size(penalty.J,1),1)];
    end
end



function pen_scaled = scaled_penalty(penalty, lambda_pen)

    pen_scaled = penalty;

    if isempty(penalty) || ~isfield(penalty,'J') || isempty(penalty.J)
        return;
    end

    pen_scaled.J = sqrt(lambda_pen) * penalty.J;
    pen_scaled.lambda_pen = lambda_pen;
end



function sweep = run_pareto_sweep(A_data, y_data, penalty, ...
                                  Aineq, bineq, lb, ub, x0, ...
                                  lambda_list, opts, model, prot, data)

    nLam = numel(lambda_list);
    nPar = numel(x0);

    phi_data = zeros(nLam,1);
    phi_pen  = zeros(nLam,1);
    phi_tot  = zeros(nLam,1);
    x_all    = zeros(nPar, nLam);
    exitflag = zeros(nLam,1);

    x_init = x0;   % warm start

    for k = 1:nLam
        lam = lambda_list(k);
        fprintf('Pareto sweep: lambda = %.3e (%d/%d)\n', lam, k, nLam);

        [A_aug, y_aug] = assemble_augmented_system(A_data, y_data, penalty, lam);

        [x_opt, ~, ~, ef] = lsqlin(A_aug, y_aug, Aineq, bineq, [], [], lb, ub, x_init, opts);

        x_all(:,k) = x_opt;
        exitflag(k) = ef;
        x_init = x_opt;

        r_data = A_data * x_opt - y_data;
        phi_data(k) = norm(r_data)^2;

        if isempty(penalty) || ~isfield(penalty,'J') || isempty(penalty.J)
            phi_pen(k) = 0.0;
        else
            r_pen = penalty.J * x_opt;
            phi_pen(k) = norm(r_pen)^2;
        end

        phi_tot(k) = phi_data(k) + lam * phi_pen(k);
    end


    idx_corner = estimate_lcurve_corner(phi_data, phi_pen);

    lambda_corner = lambda_list(idx_corner);
    lambda_scaled = lambda_corner * 0.1;
    [~, idx_corner] = min(abs(lambda_list - lambda_scaled));
   
    lambda_star = lambda_list(idx_corner);

    sweep = struct();
    sweep.lambda_list = lambda_list;
    sweep.phi_data = phi_data;
    sweep.phi_pen = phi_pen;
    sweep.phi_tot = phi_tot;
    sweep.x_all = x_all;
    sweep.exitflag = exitflag;
    sweep.idx_corner = idx_corner;
    sweep.lambda_star = lambda_star;
    sweep.model = model;
    sweep.prot = prot;
    sweep.data = data;
end


function idx = estimate_lcurve_corner(phi_data, phi_pen)

    xd = log(phi_data(:));
    yd = log(phi_pen(:));

    n = numel(xd);
    kappa = nan(n,1);

    for i = 2:n-1
        p1 = [xd(i-1), yd(i-1)];
        p2 = [xd(i),   yd(i)];
        p3 = [xd(i+1), yd(i+1)];

        a = norm(p2 - p1);
        b = norm(p3 - p2);
        c = norm(p3 - p1);

        area2 = abs(det([p2-p1; p3-p1]));

        if a > 0 && b > 0 && c > 0
            kappa(i) = 2 * area2 / (a*b*c);
        end
    end

    [~, idx] = max(kappa);
end