function model = eq_params_init_surface_I1I2(I1_max, I2_max, opts)

    model.type = 'surface_I1I2';
    model.use_polyconvex_I2 = true;

    % spline degree in each direction 
    model.k = {4,4}; 

    % number interpolation sites in each direction 
    if ~isfield(opts,'nu') || isempty(opts.nu), opts.nu = 20; end
    if ~isfield(opts,'nv') || isempty(opts.nv), opts.nv = 5; end
    nu = opts.nu;
    nv = opts.nv;

    model.I1_max = I1_max; % maximum I1 from training data (used to set interpolation sites)
    model.I2_max = I2_max; % maximum I2 from training data (used to set interpolation sites)

    % print maximum invariants (and polyconvex I2 if applicable) for reference
    fprintf('eq_params_init_surface_I1I2: I1_max = %.3f, I2_max = %.3f', I1_max, I2_max);
    if model.use_polyconvex_I2
        I2_max_poly = I2_max^(3/2) - 3*sqrt(3);
        fprintf(', I2_max_poly = %.3f', I2_max_poly);
    end
    fprintf('\n');

    % interpolation sites in (I1,I2) space
    model.I1_sites = linspace(3, I1_max, nu);
    if model.use_polyconvex_I2
        I2_max_poly = I2_max^(3/2) - 3*sqrt(3);
        model.I2_sites = linspace(0, I2_max_poly, nv);
    else
        model.I2_sites = linspace(3, I2_max, nv);
    end

    % initial guess for W(I1,I2) as sum of univariate splines in I1 and I2 directions, evaluated at (I1_sites, I2_sites)
    W1              = spapi(4, model.I1_sites, linspace(0, 20, nu));
    W2              = spapi(4, model.I2_sites, linspace(0, 40, nv));

    model.Wvals = zeros(nu, nv);
    [uu, vv]    = ndgrid(model.I1_sites, model.I2_sites);
    for j = 1:nv
        v = vv(1,j);
        W2_vals = fnval(W2, v);
        model.Wvals(:,j) = fnval(W1, uu(:,1)) + W2_vals;
    end

    % bivariate spline representation of W(I1,I2) surface
    model.Wsurf = spapi(model.k, {model.I1_sites, model.I2_sites}, model.Wvals);
  
    % bivariate spline representation of sensitivity of W(I1,I2) surface w.r.t. all parameters
    nParams = nu*nv;
    Eval    = reshape(eye(nParams), [nParams, nu, nv]);
    model.Wsurf_dtheta = spapi(model.k, {model.I1_sites, model.I2_sites}, Eval);

    % packing/unpacking
    model.pack   = @(m) m.Wvals(:);
    model.unpack = @(x,m) unpack_surface(x,m);

    % bounds: anchor at (I1=3,I2=0) with W=0
    lb = zeros(nParams,1);
    ub = Inf(nParams,1);
    lb(1) = 0; ub(1) = 0;

    model.lb = @() lb;
    model.ub = @() ub;

    % wrappers (formulation-agnostic API)
    model.eval_ab    = @(m,F) eval_ab_surface_I1I2(m,F);
    model.eval_dabdx = @(m,F) eval_dabdx_surface_I1I2(m,F);
end

function model = unpack_surface(x, model)
    % unpack x into model.Wvals, then rebuild Wsurf and Wsurf_dtheta

    nu = numel(model.I1_sites);
    nv = numel(model.I2_sites);
    model.Wvals = reshape(x, [nu, nv]);
    model.Wsurf = spapi(model.k, {model.I1_sites, model.I2_sites}, model.Wvals);
end

function [a,b] = eval_ab_surface_I1I2(model, F)
    % Given deformation gradient F, compute invariants I1, I2 and then evaluate a = dPsi/dI1, b=dPsi/dI2 using the bivariate spline surface W(I1,I2)

    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));

    if model.use_polyconvex_I2
        I2_eff = I2^(3/2) - 3*sqrt(3);
        db_dI2eff = fnval(fnder(model.Wsurf, [0 1]), [I1; I2_eff]);
        dI2eff_dI2 = (3/2)*sqrt(I2);

        a = fnval(fnder(model.Wsurf, [1 0]), [I1; I2_eff]);
        b = db_dI2eff * dI2eff_dI2;
    else
        a = fnval(fnder(model.Wsurf, [1 0]), [I1; I2]);
        b = fnval(fnder(model.Wsurf, [0 1]), [I1; I2]);
    end
end

function [da_dx, db_dx] = eval_dabdx_surface_I1I2(model, F)
    % Given deformation gradient F, compute invariants I1, I2 and then evaluate da_dx = d(dPsi/dI1)/dx and db_dx = d(dPsi/dI2)/dx where x are the parameters defining W(I1,I2) surface

    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));

    if model.use_polyconvex_I2
        I2_eff = I2^(3/2) - 3*sqrt(3);
        db_dI2eff = fnval(fnder(model.Wsurf_dtheta, [0 1]), [I1; I2_eff]);
        dI2eff_dI2 = (3/2)*sqrt(I2);

        da_dx = fnval(fnder(model.Wsurf_dtheta, [1 0]), [I1; I2_eff]);  % nParams×1
        db_dx = db_dI2eff * dI2eff_dI2;
    else
        da_dx = fnval(fnder(model.Wsurf_dtheta, [1 0]), [I1; I2]);  % nParams×1
        db_dx = fnval(fnder(model.Wsurf_dtheta, [0 1]), [I1; I2]);
    end
end
