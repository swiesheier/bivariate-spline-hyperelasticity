function model = eq_params_init(I1_max, I2_max, opts)
% Create equilibrium hyperelastic energy model either as:
%   - separable univariate splines:                 W(I1)+W(I2)
%   - coupled bivariate spline in (I1,I2) space:    W(I1,I2)
%   - coupled bivariate spline in (u,v) space:      W(u,v) where u = f(I1), v= f(I1,I2) based on admissible region of (I1,I2) space
%
% Usage:
%   model = eq_params_init(I1_max, I2_max);  % default = univariate
%   model = eq_params_init(I1_max, I2_max, struct('type','surface'));
%
% Fields returned (common API):
%   model.pack(model) -> x
%   model.unpack(x,model) -> model
%   model.lb() -> lb
%   model.ub() -> ub
%   model.eval_ab(model,F) -> [a,b]
%   model.eval_dabdx(model,F) -> [da_dx, db_dx]
%   model.nParams() -> number of params

    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'type') || isempty(opts.type), opts.type = 'univariate'; end

    switch lower(opts.type)
        case 'univariate'
            model = eq_params_init_univariate(I1_max, I2_max, opts);
        case 'surface_uv'
            model = eq_params_init_surface_uv(I1_max, I2_max, opts);
        case 'surface_i1i2'
            model = eq_params_init_surface_I1I2(I1_max, I2_max, opts);
        otherwise
            error('eq_params_init: unknown opts.type = %s', opts.type);
    end

    % Common helper
    model.nParams = @() numel(model.pack(model));
end
