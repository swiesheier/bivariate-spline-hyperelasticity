function [u, v, dv_dI1, dv_dI2, flow, fup] = map_I1I2_to_uv(I1, I2, delta, I1max, use_poly_I2)
% Map (I1,I2) to (u,v) with regularized admissible width.
%
% u = (I1 - 3)/(I1max - 3)
%
% If use_poly_I2 = false:
%   v = (I2  - flow ) / Delta_eff
%
% If use_poly_I2 = true:
%   I2p = I2^(3/2) - sqrt(3)
%   v   = (I2p - flow_p) / Delta_eff
%
% Returned dv_dI1 and dv_dI2 are always derivatives w.r.t. classical I1, I2.

    if nargin < 5
        use_poly_I2 = false;
    end

    u = (I1 - 3) / (I1max - 3);

    [flow, fup, dflow, dfup] = invariant_bounds_I2_from_I1(I1);

    if ~use_poly_I2
        % ----- classical I2 -----
        target = I2;
        flow_t = flow;
        fup_t  = fup;
        dflow_t = dflow;
        dfup_t  = dfup;
        dtarget_dI2 = 1.0;

    else
        % ----- polyconvex invariant I2p = I2^(3/2) - 3*sqrt(3) -----
        target = I2^(3/2) - 3.0 * sqrt(3);

        flow_t = flow^(3/2) - 3.0 * sqrt(3);
        fup_t  = fup^(3/2)  - 3.0 * sqrt(3);

        dflow_t = 1.5 * sqrt(flow) * dflow;
        dfup_t  = 1.5 * sqrt(fup)  * dfup;

        dtarget_dI2 = 1.5 * sqrt(I2);
    end

    Delta   = fup_t - flow_t;
    Dp      = dfup_t - dflow_t;
    Delta_eff = sqrt(Delta.^2 + delta.^2);

    num = target - flow_t;
    v   = num ./ Delta_eff;

    % dv/dI2
    dv_dI2 = dtarget_dI2 ./ Delta_eff;

    % dv/dI1
    dnum  = -dflow_t;
    dDeff = (Delta ./ Delta_eff) .* Dp;

    dv_dI1 = (dnum .* Delta_eff - num .* dDeff) ./ (Delta_eff.^2);
end