% src/protocols/protocol_SS.m
function protocol = protocol_SS()
    protocol.name = "simple_shear";

    % deformation gradient for simple shear
    protocol.Fbar = @(lambda) [1, lambda, 0; 0, 1, 0; 0, 0, 1];

    % I1/I2 as a function of lambda_max
    protocol.I1_max = @(lambda_max) 3 + lambda_max^2;
    protocol.I2_max = @(lambda_max) 3 + lambda_max^2;
end
