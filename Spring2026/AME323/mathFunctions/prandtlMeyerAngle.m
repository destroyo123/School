% Prantl-Meyer Angle "nu"
function prandtlMeyerAngle = nu(M, gamma)
% MY_STATS Calculates permutations and combinations.
% prandtlMeyerAngle = nu(M, gamma) Calculates Prantl-Meyer Angle "nu"
    prandtlMeyerAngle = sqrt((gamma+1)/(gamma-1)) * atand(sqrt(((gamma-1)/(gamma+1))*(M.^2 - 1))) - atand(sqrt(M.^2 - 1));
    % OUTPUTS DEGREES
end