%% Aero 405 Project 1
% Gerard Boberg, Ivan Cheng, and Arseniy Kotov
% 26 May 2015
%
% California Polytechnic State University, San Luis Obispo
% Aerospace engineering undergraduate program
%
% METHOD OF CHARACTERISTICS 
% Models supersonic flow thru a nozzle geometry.
% Assumes Callorically perfect gas
clc
close all
clear all
format compact


% imports
addpath( 'moc_solver' )

%% Setup Global variables.
% CHANGE THINGS HERE for different gasses.

n = 34; % number of characteristic lines


R     = 287;  % J / kg K
T0    = 2500; % K
P0    = 5e6;  % Pa    5MPa = 5x10^6
gamma = 1.4;


%% Initialize data for the algorithim

% Setup nozzle geometry
y_throat = 0.1307;     % meters, throat radius
theta_max_nozzle = 35; % degrees, maximum wall angle
n_nozzle = 100;         % number of points to render of the wall geometry
[ x_nozzle, y_nozzle ] = nozzle_geo( y_throat, theta_max_nozzle, n_nozzle );


thermo.gamma = gamma;
thermo.R     = R;
thermo.T0    = T0;

%% Run the Method of Characteristics
[ x, y, u, v, a, M ] = moc_iterative_solver( x_nozzle, y_nozzle, n,...
                                               thermo, y_throat );

x = real( x );
y = real( y );
M = real( M );
% solve out the thermodynamic properties based on Mach
[ P_nozzlethroat, P_static_wall ] = thermo_relation(...
                                   gamma, M, M(end,:), T0, P0, R );


                               
%% Plot the results
figure
surfc(x,y,M, 'EdgeAlpha', 0.1)
%contour( x, y, M )
view(2)
colorbar; 
hold on
plot(x_nozzle,y_nozzle,'LineWidth',4,'Color','k')
title( 'Mach number variation' )
xlabel('meters')
ylabel('meters')

figure
surfc(x,y,P_nozzlethroat, 'EdgeAlpha', 0.1)
%contour( x, y, M )
view(2)
colorbar; 
hold on
plot(x_nozzle,y_nozzle,'LineWidth',4,'Color','k')
title( 'Pressure variation' )
xlabel('meters')
ylabel('meters')