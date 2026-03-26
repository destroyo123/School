%% AME 323 - Gas Dynamics, Project 1: Solving the Taylor-Maccoll equations
% *Authors:* Maren Kalberer, Etan Grant
%%
% *Date:* 3-27-26
%% Description:
% So this program is gonna do XYZ 
% by doing XYZ
% 
% 
% and you can reference it like TaylorMaccol(parameter1, parameter2, etc...)
% No idea if this is how they want it but probably.
% 
% Test Inline LaTeX format: $x^2+e^{\pi i}$ 
% Test block of LaTeX: 
% 
% $$e^{\pi i} + 1 = 0$$
% 
% 
% 
% Bottom text lol <3

%% Housekeeping:
% Requires having the Aerospace Toolbox installed for certain functions

% Clear old variable values, outputs, and close old figures
close all;
clear;
clc;

%% Given Values:
% Known relationships, and input parameters.

% Atmospheric Assumptions:
R_universal =  8.314; % J/(mol K ) Universal gas constant "R"
R_air = 287; % J/(kg-degreeK) Specific gas constant for air "R_air"
gamma = 1.4; % Ratio of specific heats for standard air
rho_amb = 1.225; % Air density kg/m^3
p_amb = 1; % 1 atm? Dunno enough to know what unit we want yet I guess
T_amb = 288.15; % Std sea level temeprature is 273.15 Kelvin?

% The given mach values for which we need to generate curves:
M_inputs = [1.5, 2.0, 5.0];

%% Setting Up the System of Equations
% Sets up the functions of θ, β, M, other crap (is this actually what we wanna do?)

%% Solving the System of Equations
% Solves the system and logs the results.

%% Plotting
% This is where the required plots are generated.

%% Functions Used
% Below is a set of functions used, based on lectures, NASA, and other
% sources


%%
% *Oblique shock resultant Mach*
% 
%
% *Inputs:*
%
% * Unperturbed flow mach number 'M' (often notated as 'M1' in diagrams)
% * Oblique shock angle 'beta' (in degrees)
% * Wall/flow deflection angle 'delta' (in degrees)
% * Flow ratio of specific heats 'gamma,' usually 1.4
%
% 
%
% *Outputs:*
%
% * Mach number after oblique shock & deflection, often notated as 'M2'
%
%
%
function mach = obliqueMach(M, beta, delta, gamma)
    mach = ((( (gamma-1)*(M^2)*((sind(beta))^2)+2 ) / ( 2*gamma * (M^2) * ((sind(beta))^2) - (gamma-1) ) ) ^0.5 ) / ( sind(beta-delta) ) ;
end

%%
% *Find oblique shock angle "beta" instead of using the table.*
% 
%
% *Inputs:*
%
% * Unperturbed mach "M" (AKA 'M1' in most notation),
% * wall/flow deflection angle "theta" (often notated as "delta"),
% * ratio of specific heats "gamma" of the flow,
% * "n" specifies weak (n=0) or strong (n=1) shock case. (Usually should be weak shock case.)
%
% 
% *Outputs:*
%
% * Oblique shock angle "beta" in degrees.
%
%
%
function Beta=beta(M,theta,gamma,n)
    theta=theta*pi/180;             % convert to radians
    mu=asin(1/M);                   % Mach wave angle
    c=tan(mu)^2;
    a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta);
    b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta);
    d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
    Beta=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)))*180/pi;
    %{
    Source:
    Based on an analytical solution to the theta-beta-Mach relation given in
    the following reference:  Rudd, L., and Lewis, M. J., "Comparison of
    Shock Calculation Methods", AIAA Journal of Aircraft, Vol. 35, No. 4,
    July-August, 1998, pp. 647-649.
    
    Original script by Chris Plumley, undergraduate, University of Maryland.
    %}
end

%%
% *Find delta when given M and beta (degrees)*
% 
% *Inputs:*
%
%
% * Unperturbed Mach (AKA M1)
% * oblique shock angle beta (degrees)
% * Ratio of specific heats, 'gamma' of flow, usually 1.4
% 
% *Outputs:*
%
% * Wall deflection angle 'delta' (in degrees)
%

function delta = deltaFinder(M, beta, gamma)
    delta = ( ( M^2 * (sind(beta)^2) ) - 1) / ( M^2 * (gamma + cosd(2*beta)) + 2 );
    delta = 2*cotd(beta)*delta;
    delta = atand(delta);
end

%% Broad Steps Needed%%
% 1. Assume is a wedge and get M2 and shock angle Beta (have functionsfor
% them)
% 2. Find V immediately after shock angle and break into r and theta components 
% 3. Want to find V_theta =0 for BC



%% INPUTS %%
M1 = [1.5, 2, 5];
beta_i = deg2rad(40);         %initial guess for shock
%already have gamma defined

%% FOCUS ON M1(1) FIRST %%

tcone = deltaFinder(M1(1), beta_i, gamma);       % Now have cone half angle assumed as wedge
M2 = obliqueMach(M1(1), beta_i, tcone, gamma);   % Now have Mach after shock assumed as wedge

dtcone = @(beta) ((deltaFinder(M1(1), beta+(1*10^-8), gamma)) - (deltaFinder(M1(1), beta-(1*10^-8), gamma)) ...
    / (2*10^-8));                                 % Now have how theta cone changes with guessed beta
beta_max = fsolve(dtcone, beta_i);                % Have maximum beta possible from given mach
tcone_max = deltaFinder(M1(1), beta_max, gamma);  % Have maximum wedge half angle from max beta



%% MIGHT NOT NEED BC GEOMETRIC APPROACH %%

v_after = ((2/((gamma-1)*M2^2))+1)^-0.5;         % Now have velocity immediately after shock assumed as wedge

vr = v_after * cos(beta_i-tcone);                % Radial component
vt = v_after * sin(beta_i-tcone);                % Theta component
%%%%%%%%%%%%%%%%%%

%% TESTING %%
Results = zeros(1,3);

if vt < 0
    beta_i = beta_i + (beta_i/2);
elseif vt > 0
    beta_i = beta_i - (beta_i/2);
else
    zeros(1) = tcone;
end




%% ODE45 Function %%

function dydt = taylormaccoll(theta, y, gamma)

%convert into 1st order ODE
%radial velocity
Vr = y(1);
%equivalent to dV_r/d(theta)
dVr = y(2); 

%Derivative of y(2) equation as shown in notes
numerator = ((gamma-1)/2)*(1-Vr^2-Vt^2)*(2*Vr+Vt*cot(theta))-(Vr*Vt^2);
denominator = Vt^2-((gamma-1)/2)*(1-Vr^2-Vt^2);

dy2dt = numerator / denominator;

%Set up y-variable for ode45 as collumn
dydt = [dVr, dy2dt]';

end

%%%%%%%%%%%%%%%%





