function P = eq_P_from_F(a, b, F)
% Compute P for incompressible isotropic hyperelasticity given
% a = dW/dI1 and b = dW/dI2 at current invariants.
%
% This is formulation-agnostic: caller provides a,b.

    assert(all(size(F) == [3 3]), 'F must be 3x3');
    assert(abs(det(F) - 1) < 1e-6, 'F must be volume-preserving (det(F)=1)');


    C = F.' * F;
    I1 = trace(C);
    I  = eye(3);

    % K = (I1*I - C)
    K = I1 * I - C;

    % 2nd PK (isochoric part)
    S_bar = 2 * (a * I + b * K);

    Cinv = inv(C);
    trace_term = sum(sum(C .* S_bar));        % C : S_bar
    S_iso = S_bar - (1/3) * trace_term * Cinv;

    % 1st PK
    P_iso = F * S_iso;

    % pressure elimination with P33=0 --> p = F33*Piso33
    p = F(3,3) * P_iso(3,3);
    FinvT = inv(F).';

    P = P_iso - p * FinvT; % total stress
end