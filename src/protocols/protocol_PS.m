% src/protocols/protocol_PS.m
function protocol = protocol_PS()
    protocol.name = "planar_shear";

    % deformation gradient for planar shear
    protocol.Fbar = @(lambda) [lambda, 0, 0; 0, 1, 0; 0, 0, 1 / lambda];

    % I1/I2 as a function of lambda_max
    protocol.I1_max = @(lambda_max) 1 + lambda_max^2 + (1 / (lambda_max^2));
    protocol.I2_max = @(lambda_max) 1 + lambda_max^2 + (1 / (lambda_max^2));
end
