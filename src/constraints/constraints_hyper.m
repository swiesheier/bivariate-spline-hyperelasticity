function [Aineq, bineq] = constraints_hyper(model)
% Linear inequality constraints for:
%   - univariate model: enforce W''(I1) >= 0, W'(I1) >= 0 and (optionally) W''(I2) >= 0 and W'(I2) >= 0 (convexity and monotonicity)
%   - surface_uv model: enforce What_uu >= 0, What_u >= 0 and What_vv >= 0, What_v >= 0 (directional convexity and monotonicity)

    switch lower(model.type)
        case 'univariate'
            [Aineq, bineq] = constraints_univariate(model);

        case {'surface_uv', 'surface_i1i2'} 
            [Aineq, bineq] = constraints_surface(model);

        otherwise
            Aineq = [];
            bineq = [];
    end

    % normalize each row (constraints are defined only up to scaling)
    if ~isempty(Aineq)
        s = vecnorm(Aineq, 2, 2);
        s(s==0) = 1;
        Aineq = Aineq ./ s;
        bineq = bineq ./ s;
    end
end

% =========================== UNIVARIATE ===========================
function [Aineq, bineq] = constraints_univariate(model)
    % W(I1) convexity
    x       = model.I1_knots;
    B1      = spapi(4, x, eye(numel(x)));
    B1dd    = fnder(B1, 2);
    A1      = -B1dd.coefs';
    b1      = zeros(size(A1,1), 1);

    % W(I1) monotonicity
    B1d     = fnder(B1, 1);
    A1m     = -B1d.coefs';
    b1m     = zeros(size(A1m,1), 1);

    Aineq = [A1; A1m];
    bineq = [b1; b1m];

    if model.use_I2
        % W(I2) convexity
        B2      = spapi(4, model.I2_knots, eye(numel(model.I2_knots)));
        B2dd    = fnder(B2, 2);
        A2      = -B2dd.coefs';
        b2      = zeros(size(A2,1), 1);

        % W(I2) monotonicity
        B2d     = fnder(B2, 1);
        A2m     = -B2d.coefs';
        b2m     = zeros(size(A2m,1), 1);

        A2 = [A2; A2m];
        b2 = [b2; b2m];

        Aineq = blkdiag(Aineq, A2);
        bineq = [bineq; b2];
    end
end

% ============================ SURFACE =============================
function [Aineq, bineq] = constraints_surface(model)
    % Requires model.Wsurf_dtheta: vector-valued spapi surface where component k is dWhat/dx_k
    assert(isfield(model,'Wsurf_dtheta') && ~isempty(model.Wsurf_dtheta), ...
        'constraints_surface: model.Wsurf_dtheta missing.');

    % Directional convexity: What_uu >= 0 and What_vv >= 0
    Buu = fnder(model.Wsurf_dtheta, [2 0]);
    Bvv = fnder(model.Wsurf_dtheta, [0 2]);

    % Convert coefficient tensor to linear constraints A*x <= 0
    eps_uu = 0.0; % small tolerance to help with conditioning
    eps_vv = 0.0;
    Auu = -reshape(Buu.coefs, size(Buu.coefs,1), []).';
    Avv = -reshape(Bvv.coefs, size(Bvv.coefs,1), []).';
    buu = eps_uu * ones(size(Auu,1), 1);
    bvv = eps_vv * ones(size(Avv,1), 1);


    % Directional monotonicity: What_u >= 0 and What_v >= 0 
    Bu = fnder(model.Wsurf_dtheta, [1 0]);
    Bv = fnder(model.Wsurf_dtheta, [0 1]);
    Au = -reshape(Bu.coefs, size(Bu.coefs,1), []).';
    Av = -reshape(Bv.coefs, size(Bv.coefs,1), []).';
    bu = zeros(size(Au,1), 1);
    bv = zeros(size(Av,1), 1);

    Aineq = [Auu; Avv; Au; Av];
    bineq = [buu; bvv; bu; bv];
end