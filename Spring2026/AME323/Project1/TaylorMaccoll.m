%% AME 323 - Gas Dynamics, Project 1: Solving the Taylor-Maccoll equations
% (Text below the table of contents in the PDF generated when you hit
% "publish")

%% Description:
%{
So this program is gonna do XYZ
by doing XYZ

and you can reference it like TaylorMaccol(parameter1, parameter2, etc...)
No idea if this is how they want it but probably.

Bottom text lol <3
%}
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
% sources.