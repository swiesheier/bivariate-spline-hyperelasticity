function dP_dx = eq_dPdx_from_F(da_dx, db_dx, F)
% Compute Jacobian d vec(P) / d x for incompressible isotropic hyperelasticity
% given parameter sensitivities of a=dW/dI1 and b=dW/dI2:
%   da_dx: nParams x 1
%   db_dx: nParams x 1
%
% Output:
%   dP_dx : 9 x nParams, column k is dP(:)/dx_k

    assert(all(size(F) == [3 3]), 'F must be 3x3');
    assert(abs(det(F) - 1) < 1e-6, 'F must be volume-preserving (det(F)=1)');

    % invariants / kinematics
    C     = F.' * F;
    I1    = trace(C);
    I     = eye(3);
    K     = I1 * I - C;
    Cinv  = inv(C);
    FinvT = inv(F).';

    nParams = numel(da_dx);
    assert(numel(db_dx) == nParams, 'da_dx and db_dx must match length.');

    dP_dx = zeros(9, nParams);

    dSbar_da = 2 * I;
    dSbar_db = 2 * K;

    for k = 1:nParams
        da = da_dx(k);
        db = db_dx(k);

        dS_bar = dSbar_da * da + dSbar_db * db;

        dtrace = sum(sum(C .* dS_bar));     % C : dS_bar
        dS_iso = dS_bar - (1/3) * dtrace * Cinv;

        dP_iso = F * dS_iso;

        dp = F(3,3) * dP_iso(3,3);

        dP = dP_iso - dp * FinvT;
        dP_dx(:, k) = dP(:);
    end
end
