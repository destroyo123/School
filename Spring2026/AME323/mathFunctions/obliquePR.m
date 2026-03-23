% Oblique Shock Pressure Ratio
function pressureRatio = obliquePR(M, gamma, beta)
    pressureRatio = 1 + ( (2*gamma)/( gamma+1 ) * ( (M^2) * (sind(beta)^2) - 1 ) );
end