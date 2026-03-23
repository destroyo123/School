% Mach from prandtl-meyer angle (IN DEGREES)
function mach = meyerMach(gamma, nu)
    mach = flowprandtlmeyer(gamma, nu, 'nu');
end