function [flow, fup, dflow_dI1, dfup_dI1] = invariant_bounds_I2_from_I1(I1)
% Compute sharp bounds I2 = flow(I1), fup(I1) for admissible incompressible invariants
% using D(I1,I2)=0 (discriminant) and implicit differentiation.
%
% Paper: Dammaß et al. (2025), Theorem 2 and Eq. (32)

    if any(I1(:) < 3)
        error('I1 must be >= 3 for incompressible admissible invariants.');
    end

    % Vectorized over I1
    flow      = zeros(size(I1));
    fup       = zeros(size(I1));
    dflow_dI1 = zeros(size(I1));
    dfup_dI1  = zeros(size(I1));

    for idx = 1:numel(I1)
        a = I1(idx);

        % Discriminant polynomial P(a, I2) = 0 (scaled by 27)
        % From Eq. (32): I1^3 + I2^3 - 1/4 I1^2 I2^2 - 9/2 I1 I2 + 27/4 = 0
        % => I2^3 - (1/4)a^2 I2^2 - (9/2)a I2 + (a^3 + 27/4) = 0
        coeff = [ ...
            1, ...
            -0.25*a^2, ...
            -4.5*a, ...
            (a^3 + 27/4) ...
        ];

        r = roots(coeff);

        % keep real roots (numerically)
        r = r(abs(imag(r)) < 1e-6);
        r = real(r);

        % keep those in the admissible "B" estimate region roughly (>=3)
        r = r(r >= 3 - 1e-10);

        if numel(r) < 2
            error('Could not find two admissible real roots for I1=%.6g', a);
        end

        r = sort(r(:));
        lo = r(1);
        up = r(end);

        flow(idx) = lo;
        fup(idx)  = up;

        % Implicit derivative f'(I1) = -P_I1 / P_I2 at (I1, I2=f(I1))
        % P(I1,I2) = I2^3 - 1/4 I1^2 I2^2 - 9/2 I1 I2 + I1^3 + 27/4
        %
        % P_I2 = 3 I2^2 - 1/2 I1^2 I2 - 9/2 I1
        % P_I1 = -(1/2) I1 I2^2 - 9/2 I2 + 3 I1^2

        dflow_dI1(idx) = - P_I1(a, lo) / P_I2(a, lo);
        dfup_dI1(idx)  = - P_I1(a, up) / P_I2(a, up);
    end
end

function val = P_I2(I1, I2)
    val = 3*I2.^2 - 0.5*I1.^2.*I2 - 4.5*I1;
end

function val = P_I1(I1, I2)
    val = -0.5*I1.*I2.^2 - 4.5*I2 + 3*I1.^2;
end
