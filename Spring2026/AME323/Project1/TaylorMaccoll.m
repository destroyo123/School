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

%% Setting up the arrays and variables
% Sets up system of θ, β, M, other crap.

% number of curves to plot (same as number of mach inputs)
num_curves = length(M_inputs);

% Range of theta_c (cone half-angle θ) input valves to plot
theta_c_min = 0.5; % Left bound of the curve (degrees)
theta_c_max = 15; % Right bound of the curve (degrees)
num_points = 4;

theta_c = linspace(theta_c_min, theta_c_max, num_points); % 'x' coord values

% Set up an empty array for the output β values
% Each row: One M input
% Each column: The values of theta_c (see above definitions)
beta_outputs = zeros(num_curves, num_points);

% Finds the cone shock angles given the above-defined M, θ inputs.
for i = 1:num_curves
    
    % For each row (M value), run the beta, point-by point.
    % HORRIBLY inefficient, but we just need it to work first.
    for j = 1:num_points

        % Grab the input
        freestream_mach = M_inputs(i);
        cone_angle = theta_c(j);

        % Store the Output
        beta_outputs(i,j) = coneBeta(freestream_mach, cone_angle, gamma);

    end
end


%% Solving the System of Equations
% Solves the system and logs the results.

%% Plotting
% This is where all plots are generated.

% Set interpereter to Latex. Lets us use subscripts and greek characters.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

plot1 = true; % Toggle this plot on/off
if plot1
    fig1 = figure;
    ax1 = gca;
    
    % Text format
    lineWidth = 1.25;
    axisFontSize = 18;
    legendFontSize = 14;
    hold on;
    
    colors = ['r', 'g', 'b']; % colors to cycle through
    for i = 1:num_curves
        % Make a custon name "M = X.XX" for each curve and store for the legend
        name = "$M = " + sprintf("%.2f", M_inputs(i))+"$";
    
        % Grab the x 'θ' and y 'β' values
        x = theta_c; % Row vector of cone angles
        y = beta_outputs(i,:); % Row vector of shock angles
    
        % Cycle through the colors (however many colors are defined)
        color = colors(mod(i-1,length(colors))+1);
        
        % Plot the given curve, including the name (referenced in legend)
        plot(x, y, "LineWidth", lineWidth, "Color", color, DisplayName=name);
    end
    
    grid on;
    
    % Put the legend and axis labels
    legend(FontSize=legendFontSize,Location="northwest");
    xlabel("$\theta_{cone}$", fontsize = axisFontSize);
    ylabel("$\beta_{cone}$", FontSize = axisFontSize);
    
    % Sets axes at origin.
    ax1.XAxisLocation = 'origin';
    ax1.YAxisLocation = 'origin';
    
    hold off; % Done with plot :)
end

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
% 3. Want to find V_theta = 0 for BC

function outputbeta = coneBeta(Mach, theta, gamma)
    %% 1. Setup Inputs %%
    M1 = Mach;
    theta_cone_input = theta; % degrees

    %% 2. Initial Guess %%
    % Get an initial guess for the shock angle beta (assuming it's a wedge)
    % This provides a starting point for fzero
   

    %% 3. Shooting Method (The Fix) %%
    % fzero tries different beta values until the error is zero
    options = optimset('TolX', 1e-4); % Loosen the tolerance for performance

    %beta_low = asind(1/M1)+0.1; % low-end guess
    beta_low = asind(1/M1)+0.1; % low-end guess
    beta_high = 89; % high-end guess (tee-hee)

    outputbeta = fzero(@(b) solve_for_theta(b, M1, gamma, theta_cone_input), [beta_low beta_high]);

end % End of main function

%% --- HELPER FUNCTIONS --- %%

function error = solve_for_theta(beta_guess, M1, gamma, theta_target)

    % 1. Initial conditions at shock boundary (Dimensionless)
    beta_rad = deg2rad(beta_guess);
    
    delta = deltaFinder(M1, beta_guess, gamma);
    M2 = obliqueMach(M1, beta_guess, delta, gamma);

    % V_prime is the velocity magnitude normalized by V_max
    V_prime = ( (2 / ((gamma - 1) * M2^2)) + 1 )^-0.5;

    % Density ratio across the shock
    M1n2 = (M1 * sin(beta_rad))^2;
    eps = ((gamma - 1) * M1n2 + 2) / ((gamma + 1) * M1n2);
    
    theta_rad = deg2rad(theta_target); % Bookkeeping teehee
    beta_rad = deg2rad(beta_guess);

    % Initial velocity components (Vr and Vt)
    vr0 = V_prime * cos(beta_rad-theta_rad); 
    vt0 = V_prime  * sin(beta_rad-theta_rad); 

    % 2. Integrate ODE from shock (beta) down toward axis (0)
    y0 = [vr0; vt0];
    tspan = [beta_rad, deg2rad(0.1)];
    
    % Solve the Taylor-Maccoll ODE using ODE45
    % but set the tolerances to be looser
    ode45options = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);
    [theta_out, y_out] = ode45(@(t, y) taylormaccoll(t, y, gamma), tspan, y0, ode45options);
    
    
    % 3. Find where Vt crosses zero (the physical cone surface)
    Vt_vals = y_out(:,2);
    
    % Find zero crossing
    idx = find(Vt_vals(1:end-1).*Vt_vals(2:end) < 0, 1);
    
    if isempty(idx)
        theta_found = theta_out(end);
    else
        theta_found = theta_out(idx);
    end
    
    % 4. Return the difference between solved theta and target theta
    % except note that theta_target is input in degrees.
    error = theta_found - deg2rad(theta_target); % Outputs radians
end

function dydt = taylormaccoll(theta, y, gamma)
    Vr = y(1);
    Vt = y(2); 

    % The Taylor-Maccoll Equation
    numerator = ((gamma-1)/2)*(1-Vr^2-Vt^2)*(2*Vr+Vt*cot(theta))-(Vr*Vt^2);
    denominator = Vt^2-((gamma-1)/2)*(1-Vr^2-Vt^2);

    dVr = Vt;
    dVt = numerator / denominator;

    dydt = [dVr; dVt];

end