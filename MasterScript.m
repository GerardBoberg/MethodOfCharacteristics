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

n = 5; % number of characteristic lines


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
[ x, y, u, v, a, M ] = moc_iterative_solver( x_nozzle, y_nozzle, n, thermo );

%figure();
%plot( x_nozzle, y_nozzle, 'b', x(:,1), y(:,1), 'r' );
%axis equal;