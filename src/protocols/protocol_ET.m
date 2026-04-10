% src/protocols/protocol_ET.m
function protocol = protocol_ET()
    protocol.name = "biax";

    % deformation gradient for equi-biaxial stretch
    protocol.Fbar = @(lambda) diag([lambda, lambda, 1 / (lambda^2)]);

    % I1/I2 as a function of lambda_max
    protocol.I1_max = @(lambda_max) 2 * lambda_max^2 + (1 / (lambda_max^4));
    protocol.I2_max = @(lambda_max) lambda_max^4 + 2 * (1 / (lambda_max^2));
end
