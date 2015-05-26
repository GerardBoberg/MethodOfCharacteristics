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

n = 125; % number of characteristic lines


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


% Smooth out the points along the symmetry line
L = size( x, 2 );
x_prime = [ zeros( 1, L ); x ];
y_prime = [ zeros( 1, L ); y ];
M_prime = [ zeros( 1, L ); M ];

for ii = 1:size( x_prime, 1 )
    for jj = 2:size( x_prime, 2 )
        if( x_prime( ii, jj ) == 0 )
            x_prime( ii, jj ) = NaN;
        end
        if( M_prime( ii, jj ) == 0 )
            M_prime( ii, jj ) = NaN;
        end
    end
end
for kk = 1:L
    if( (y_prime(2,kk) > 0 ) &&...
        (y_prime(2, kk-1) == 0) && (y_prime(2,kk+1) == 0) )
        x_prime( 1, kk ) = (x_prime(2,kk-1) + x_prime(2,kk+1))/2;
        y_prime( 1, kk ) = 0;
        M_prime( 1, kk ) = (M_prime(2,kk-1) + M_prime(2,kk+1))/2;
        
        if( ~(x_prime( 1, kk-1 )  > 0) )
            x_prime( 1, kk-1 ) = x_prime( 1, kk );
            M_prime( 1, kk-1 ) = M_prime( 1, kk );
            y_prime( 1, kk-1 ) = y_prime( 1, kk );
        end
    end
end


% solve out the thermodynamic properties based on Mach
[ P_nozzlethroat, P_static_wall ] = thermo_relation(...
                                   gamma, M_prime, M(end,:), T0, P0, R );

                               
%% Plot the results

% 1 -- Mach number variation
figure
C = min( M_prime, 5 );
surfc(x_prime,y_prime,M_prime, C, 'EdgeAlpha', 0.1)
axis( [0, 0.5, 0, 0.14] );
caxis( [1,6] );
colormap( hsv );
%contour( x, y, M )
view(2)
colorbar; 
hold on
plot(x_nozzle,y_nozzle,'LineWidth',4,'Color','k')
title( 'Mach number variation in the nozzle' )
xlabel('location in nozzle, meters')
ylabel('location in nozzle, meters')

% 2 -- Pressure variation
figure
surfc(x_prime,y_prime,P_nozzlethroat, 'EdgeAlpha', 0.1)
%contour( x, y, M )
axis( [0, 0.5, 0, 0.14] );
view(2)
colorbar; 
colormap( hsv );
caxis( [0, 2.5e6] )
hold on
plot(x_nozzle,y_nozzle,'LineWidth',4,'Color','k')
title( 'Pressure variation in the nozzle' )
xlabel('nozzle location, meters')
ylabel('nozzle location, meters')

% 3 -- P_static along the wall vs. x-location
figure();
plot( x(end,1:length(P_static_wall)), P_static_wall, 'r-x');
title( 'Pressure variation along the wall' );
xlabel( 'wall x location, meters' );
ylabel( 'Pressure, Pa' );

% 4 -- Characteristic Line intersections
figure();
hold on;
plot( x_prime, y_prime, 'c.' );
axis( [0, 0.5, 0, 0.14] );
plot(x_nozzle,y_nozzle,'LineWidth',2,'Color','k')
title( 'Characteristic Line intersections' );
xlabel( 'nozzle x location, meters' );
ylabel( 'nozzle y location, meters' );

% 5 -- Final Nozzle geometry w/ points where char lines intersected wall
figure();
hold on;
plot(x_nozzle,y_nozzle,'LineWidth',1,'Color','k')
plot( x(end,:), y(end,:), 'rx' );
axis( [0, 0.04, 0.125, 0.145] );
title( 'Characteristic line intersection with nozzle wall' );
xlabel( 'nozzle x location, meters' );
ylabel( 'nozzle y location, meters' );

% 6 -- output of the values of the symmetry intersections
sym_loc( :, 1 ) = x( 1, 1:2:ceil(2.2*n) );
sym_loc( :, 2 ) = y( 1, 1:2:ceil(2.2*n) );
sym_loc( :, 3 ) = M( 1, 1:2:ceil(2.2*n) );
display( '---- symmetry intersections -----' )
display( '  x            y         Mach' );
display( num2str( sym_loc, 3 ) );

% 7 -- quiver for debug purposes
figure();
quiver( x, y, u, v );
title( 'Quiver plot of velocities' );
xlabel( 'Nozzle x location, meters' );
ylabel( 'Nozzle y location, meters' );
