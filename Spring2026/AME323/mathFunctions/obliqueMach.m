% Oblique shock resultant Mach
function mach = obliqueMach(M, beta, delta, gamma)
    mach = ((( (gamma-1)*(M^2)*((sind(beta))^2)+2 ) / ( 2*gamma * (M^2) * ((sind(beta))^2) - (gamma-1) ) ) ^0.5 ) / ( sind(beta-delta) ) ;
end