%% AME323 Discussion 306

% Housekeeping:
close all;
clear;
clc;

%% The Easy Part
% Applying the given data to find the fixed parts of the problem

% Given data:
gamma = 1.4;
M1 = 2.8;
beta12 = 35; % degrees
beta13 = 40; % degrees


% First, find the flow angles at 2 and 3
delta12 = deltaFinder(M1, beta12, gamma);
delta13 = deltaFinder(M1, beta13, gamma);

% Then use those to find the mach at 2 and 3
M2 = obliqueMach(M1, beta12, delta12, gamma);
M3 = obliqueMach(M1, beta13, delta13, gamma);

% Also the pressure at 2 and 3
p2p1 = pr(M1, gamma, beta12); % "p2 / p1"
p3p1 = pr(M1, gamma, beta13); % "p3 / p1"

% CORRECT UP TO THIS POINT SO FAR!!!

%% The Hard Part

% Define local pressure ratios as functions of local beta
p5p3 = @(b35) pr(M3, gamma, b35);
p4p2 = @(b24) pr(M2, gamma, b24);

% Define local deflections as functions of local beta
d35 = @(b35) deltaFinder(M3, b35, gamma);
d24 = @(b24) deltaFinder(M2, b24, gamma);

% The System of Equations
% v(1) = local beta at state 2 (beta24)
% v(2) = local beta at state 3 (beta35)
fun = @(v) [
    p2p1 * p4p2(v(1)) - p3p1 * p5p3(v(2)); ...          % p4 = p5
    (delta12 - d24(v(1))) - (-delta13 + d35(v(2))) ...  % theta4 = theta5
];

% Initial Guess: Local shock angles are probably similar to the first ones
v0 = [beta12, beta13]; 

[v_result, ~] = fsolve(fun, v0);

% Results Extraction
beta24_local = v_result(1);
beta35_local = v_result(2);

% Calculate final shared properties
p_final = p2p1 * p4p2(beta24_local);
theta_final = delta12 - d24(beta24_local); % Should match -delta13 + d35...

fprintf("\n");
fprintf("------------------------------------------------------------\n");
fprintf("Results:\n")
fprintf("Final Pressure Ratio (p4/p1 or p5/p1): %.3f\n", p_final);
fprintf("Final Flow Angle δ' (degrees): %.3f\n", theta_final);
% fprintf("Local Beta 2-4: %.3f\n", beta24_local);
% fprintf("Local Beta 3-5: %.3f\n", beta35_local);
fprintf("------------------------------------------------------------\n");


%% Functions
% Equations and Relations I mostly made for other homeworks & discussion
% problems.

% Prantl-Meyer Angle "nu"
function prandtlMeyerAngle = nu(M, gamma)
    prandtlMeyerAngle = sqrt((gamma+1)/(gamma-1)) * atand(sqrt(((gamma-1)/(gamma+1))*(M.^2 - 1))) - atand(sqrt(M.^2 - 1));
    % OUTPUTS DEGREES
end

% Oblique Shock Pressure Ratio
function pressureRatio = pr(M, gamma, beta)
    pressureRatio = 1 + ( (2*gamma)/( gamma+1 ) * ( (M^2) * (sind(beta)^2) - 1 ) );
end

% Oblique shock resultant Mach
function mach = obliqueMach(M, beta, delta, gamma)
    mach = ((( (gamma-1)*(M^2)*((sind(beta))^2)+2 ) / ( 2*gamma * (M^2) * ((sind(beta))^2) - (gamma-1) ) ) ^0.5 ) / ( sind(beta-delta) ) ;
end

% Beta finder
function Beta=beta(M,theta,gamma,n)
    theta=theta*pi/180;             % convert to radians
    mu=asin(1/M);                   % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta);
    d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
    Beta=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)))*180/pi;
    % OUTPUTS DEGREES
end

% Mach from prandtl-meyer angle (IN DEGREES)
function mach = meyerMach(gamma, nu)
    mach = flowprandtlmeyer(gamma, nu, 'nu');
end

% Find delta when given M and beta (degrees)
function delta = deltaFinder(M, beta, gamma)
    delta = ( ( M^2 * (sind(beta)^2) ) - 1) / ( M^2 * (gamma + cosd(2*beta)) + 2 );
    delta = 2*cotd(beta)*delta;
    delta = atand(delta);
end
