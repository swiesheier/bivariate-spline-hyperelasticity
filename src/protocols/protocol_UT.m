% src/protocols/protocol_UT.m
function protocol = protocol_UT()
    protocol.name = "uniax";

    % deformation gradient for uniaxial stretch
    protocol.Fbar = @(lambda) diag([lambda, 1 / sqrt(lambda), 1 / sqrt(lambda)]);

    % I1/I2 as a function of lambda_max
    protocol.I1_max = @(lambda_max) lambda_max^2 + 2 * (1 / lambda_max);
    protocol.I2_max = @(lambda_max) 2 * lambda_max + (1 / lambda_max^2);
end
