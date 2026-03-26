%%
% *Prantl-Meyer Angle "nu"* 
% 
% *Inputs:* 
%
% * Local flow mach number 'M'
% * Local flow ratio of specific heats 'gamma' (usually 1.4)
% 
% 
% *Outputs:* 
%
%
% * Prandtl-Meyer angle "nu", in degrees.
%

function prandtlMeyerAngle = nu(M, gamma)
    prandtlMeyerAngle = sqrt((gamma+1)/(gamma-1)) * atand(sqrt(((gamma-1)/(gamma+1))*(M.^2 - 1))) - atand(sqrt(M.^2 - 1));
end