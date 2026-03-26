% *Mach from prandtl-meyer angle (IN DEGREES)*
% 
% *Inputs:*
%
%
% * Ratio of specific heats "gamma" of flow, usually 1.4
% * Prandtl-Meyer angle "nu" (degrees)
% 
% *Outputs:*
%
% * Local mach number "mach"
%

function mach = meyerMach(gamma, nu)
    mach = flowprandtlmeyer(gamma, nu, 'nu');
end