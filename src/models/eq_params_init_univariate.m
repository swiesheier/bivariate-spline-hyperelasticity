function model = eq_params_init_univariate(I1_max, I2_max, opts)

    model.type = 'univariate';
    if ~isfield(opts,'use_I2') || isempty(opts.use_I2)
        model.use_I2 = true;
        model.use_polyconvex_I2 = true;
    else
        model.use_I2 = logical(opts.use_I2);
    end

    % print maximum invariants (and polyconvex I2 if applicable) for reference
    fprintf('eq_params_init_univariate: I1_max = %.3f, I2_max = %.3f', I1_max, I2_max);
    if model.use_I2 && model.use_polyconvex_I2
        I2_max_poly = I2_max^(3/2) - 3*sqrt(3);
        fprintf(', I2_max_poly = %.3f', I2_max_poly);
    end
    fprintf('\n');
    

    model.k  = 4;  % cubic splines
    model.n1 = 20; % number of interpolation points for W(I1)
    model.n2 = 5;  % number of interpolation points for W(I2)

    % --- W(I1) ---
    model.I1_knots = linspace(3, I1_max, model.n1); % sites
    model.I1_vals  = linspace(0, 20, model.n1);     % initial guess
    model.WI1_spline         = spapi(model.k, model.I1_knots, model.I1_vals); % spline representation of W(I1)
    model.WI1_spline_dtheta1 = spapi(model.k, model.I1_knots, eye(model.n1)); % parameter sensitivity of W(I1) w.r.t. I1_vals

    % --- W(I2) ---
    if model.use_I2
        if model.use_polyconvex_I2
            I2_max_poly = I2_max^(3/2) - 3*sqrt(3);
            model.I2_knots = linspace(0, I2_max_poly, model.n2); % sites
        else
            model.I2_knots = linspace(3, I2_max, model.n2); % sites
        end
        model.I2_vals  = linspace(0, 40, model.n2); % initial guess
        model.WI2_spline         = spapi(model.k, model.I2_knots, model.I2_vals); % spline representation of W(I2)
        model.WI2_spline_dtheta2 = spapi(model.k, model.I2_knots, eye(model.n2)); % parameter sensitivity of W(I2) w.r.t. I2_vals
    else
        model.I2_knots = [];
        model.I2_vals  = [];
        model.WI2_spline = [];
        model.WI2_spline_dtheta2 = [];
    end

    % --- Packing / Unpacking ---
    model.pack   = @(m) pack_model(m);
    model.unpack = @(x,m) unpack_model(x,m);

    % --- Bounds ---
    model.lb = @() build_lb(model);
    model.ub = @() build_ub(model);

    % --- Formulation-agnostic wrappers ---
    model.eval_ab    = @(m,F) eval_ab_univariate(m,F);
    model.eval_dabdx = @(m,F) eval_dabdx_univariate(m,F);
end

function x = pack_model(model)
    x = model.I1_vals(:);
    if model.use_I2
        x = [x; model.I2_vals(:)];
    end
end

function model = unpack_model(x, model)
    % unpack x into model.I1_vals and model.I2_vals, then rebuild splines
    n1 = numel(model.I1_vals);
    model.I1_vals = reshape(x(1:n1), size(model.I1_vals));

    if model.use_I2
        n2 = numel(model.I2_vals);
        model.I2_vals = reshape(x(n1+1:n1+n2), size(model.I2_vals));
    end

    % rebuild splines with updated values
    model.WI1_spline = spapi(model.k, model.I1_knots, model.I1_vals);
    if model.use_I2
        model.WI2_spline = spapi(model.k, model.I2_knots, model.I2_vals);
    end
end

function lb = build_lb(model)
    % enforce W(I1=3)=0 and W(I2=0)=0 by setting first value to 0 and rest >= 0
    lb = [0.0; zeros(numel(model.I1_vals)-1,1)];
    if model.use_I2
        lb = [lb; 0.0; zeros(numel(model.I2_vals)-1,1)];
    end
end

function ub = build_ub(model)
    % enforce W(I1=3)=0 and W(I2=0)=0 by setting first value to 0 and rest <= Inf
    ub = [0.0; Inf(numel(model.I1_vals)-1,1)];
    if model.use_I2
        ub = [ub; 0.0; Inf(numel(model.I2_vals)-1,1)];
    end
end

function [a,b] = eval_ab_univariate(model, F)
    % Given deformation gradient F, compute invariants I1, I2 and then evaluate a = dPsi/dI1, b=dPsi/dI2

    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));
    J = sqrt(det(C));

    a = fnval(fnder(model.WI1_spline, 1), I1); % dW(I1)/dI1
    if model.use_I2
        if model.use_polyconvex_I2
            I2_eff = I2^(3/2) - 3*sqrt(3);
            db_dI2eff = fnval(fnder(model.WI2_spline, 1), I2_eff);
            dI2eff_dI2 = (3/2)*sqrt(I2);
            b = db_dI2eff * dI2eff_dI2;
        else
            b = fnval(fnder(model.WI2_spline, 1), I2); % dW(I2)/dI2
        end
    else
        b = 0.0;
    end
end

function [da_dx, db_dx] = eval_dabdx_univariate(model, F)
    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));

    n1 = model.WI1_spline_dtheta1.number;
    da_dx1 = fnval(fnder(model.WI1_spline_dtheta1, 1), I1); % dW(I1)/dI1 sensitivity w.r.t. I1_vals

    if model.use_I2
        n2 = model.WI2_spline_dtheta2.number;
        if model.use_polyconvex_I2
            I2_eff = I2^(3/2) - 3*sqrt(3);
            db_dI2eff = fnval(fnder(model.WI2_spline_dtheta2, 1), I2_eff);
            dI2eff_dI2 = (3/2)*sqrt(I2);
            db_dx2 = db_dI2eff * dI2eff_dI2;
        else
            db_dx2 = fnval(fnder(model.WI2_spline_dtheta2, 1), I2); % dW(I2)/dI2 sensitivity w.r.t. I2_vals
        end
    else
        n2 = 0;
        db_dx2 = zeros(0,1);
    end

    nParams = n1 + n2;
    da_dx = zeros(nParams,1);
    db_dx = zeros(nParams,1);

    da_dx(1:n1) = da_dx1;
    if n2 > 0
        db_dx(n1+1:end) = db_dx2;
    end
end
