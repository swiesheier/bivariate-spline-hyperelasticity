function model = eq_params_init_surface_uv(I1_max, I2_max, opts)

    model.type = 'surface_uv';
    model.use_polyconvex_I2 = true;

    % spline degree in each direction 
    model.k = {4,4}; 

    % number interpolation sites in each direction 
    if ~isfield(opts,'nu') || isempty(opts.nu), opts.nu = 20; end % 20
    if ~isfield(opts,'nv') || isempty(opts.nv), opts.nv = 5; end % 5
    nu = opts.nu;
    nv = opts.nv;

    model.I1_max = I1_max;  % maximum I1 from training data, used for affine mapping to u in [0,1]
    model.delta = 3e-5;     % regularization of admissible band at origin in (I1,I2) space, used for mapping to v in [0,1]

    % print maximum invariants (and polyconvex I2 if applicable) for reference
    fprintf('eq_params_init_surface_uv: I1_max = %.3f', I1_max);
    if model.use_polyconvex_I2
        I2_max_poly = I2_max^(3/2) - 3*sqrt(3);
        fprintf(', I2_max = %.3f, I2_max_poly = %.3f', I2_max, I2_max_poly);
    end
    fprintf('\n');

    % interpolation sites in (u,v) space
    model.u_sites = linspace(0, 1, nu); % u in [0,1] maps to I1 in [3, I1max]
    model.v_sites = linspace(0, 1, nv); % v in [0,1] maps to I2 in [flow(I1), fup(I1)] via regularized mapping

    % initial guess
    W1              = spapi(4, model.u_sites, linspace(0, 20, nu));
    [~, fup_max]    = invariant_bounds_I2_from_I1(I1_max);
    I2_sites        = linspace(3, fup_max, nu);
    W2              = spapi(4, I2_sites, linspace(0, 40, nu));

    model.Wvals = zeros(nu, nv);
    [uu, vv]    = ndgrid(model.u_sites, model.v_sites);
    for j = 1:nv
        v = vv(1,j);
        I2_mapped = 3 + v * (fup_max - 3);
        W2_vals = fnval(W2, I2_mapped);
        model.Wvals(:,j) = fnval(W1, uu(:,1)) + W2_vals;
    end

    % bivariate spline representation of W(u,v) surface
    model.Wsurf = spapi(model.k, {model.u_sites, model.v_sites}, model.Wvals);
  
    % bivariate spline representation of sensitivity of W(u,v) surface w.r.t. all parameters
    nParams = nu*nv;
    Eval    = reshape(eye(nParams), [nParams, nu, nv]);
    model.Wsurf_dtheta = spapi(model.k, {model.u_sites, model.v_sites}, Eval);

    % packing/unpacking
    model.pack   = @(m) m.Wvals(:);
    model.unpack = @(x,m) unpack_surface(x,m);

    % bounds: anchor at (u=3,v=0)
    lb = zeros(nParams,1);
    ub = Inf(nParams,1);
    lb(1) = 0; ub(1) = 0;

    model.lb = @() lb;
    model.ub = @() ub;

    % wrappers (formulation-agnostic API)
    model.eval_ab    = @(m,F) eval_ab_surface_uv(m,F,model.delta);
    model.eval_dabdx = @(m,F) eval_dabdx_surface_uv(m,F,model.delta);
end

function model = unpack_surface(x, model)
    % unpack x into model.Wvals, then rebuild Wsurf (fixed sites) --- Wsurf_dtheta unchanged since fixed sites

    nu = numel(model.u_sites);
    nv = numel(model.v_sites);
    model.Wvals = reshape(x, [nu, nv]);
    model.Wsurf = spapi(model.k, {model.u_sites, model.v_sites}, model.Wvals);
end

function [a,b] = eval_ab_surface_uv(model, F, delta)
    % Given deformation gradient F, compute invariants I1, I2 and then evaluate a = dPsi/dI1, b=dPsi/dI2 

    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));

    % I1max/min 
    I1min = 3.0;
    I1max = model.I1_max;

    [u, v_eff, dvEff_dI1, dvEff_dI2] = map_I1I2_to_uv(I1, I2, delta, I1max, model.use_polyconvex_I2);

    % clamp (u,v_eff) to be within interpolation range to avoid NaNs from "fnval"
    u     = min(max(u,     model.u_sites(1)), model.u_sites(end));
    v_eff = min(max(v_eff, model.v_sites(1)), model.v_sites(end));

    Wu = fnval(fnder(model.Wsurf, [1 0]), [u; v_eff]);
    Wv = fnval(fnder(model.Wsurf, [0 1]), [u; v_eff]);

    a = Wu ./ (I1max - I1min) + Wv * dvEff_dI1;
    b = Wv * dvEff_dI2;
end

function [da_dx, db_dx] = eval_dabdx_surface_uv(model, F, delta)
    % Given deformation gradient F, compute invariants I1, I2 and then evaluate da_dx = d(dPsi/dI1)/dx and db_dx = d(dPsi/dI2)/dx where x are the parameters defining W(u,v) surface
    
    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));

    I1min = 3.0;
    I1max = model.I1_max;
    [u, v_eff, dvEff_dI1, dvEff_dI2] = map_I1I2_to_uv(I1, I2, delta, I1max, model.use_polyconvex_I2);

    u     = min(max(u,     model.u_sites(1)), model.u_sites(end));
    v_eff = min(max(v_eff, model.v_sites(1)), model.v_sites(end));

    Phi_u = fnval(fnder(model.Wsurf_dtheta, [1 0]), [u; v_eff]);  % nParams×1
    Phi_v = fnval(fnder(model.Wsurf_dtheta, [0 1]), [u; v_eff]);  % nParams×1

    da_dx = Phi_u ./ (I1max - I1min) + Phi_v * dvEff_dI1;
    db_dx = Phi_v * dvEff_dI2;
end
