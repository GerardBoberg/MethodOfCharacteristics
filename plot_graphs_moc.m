function [] = plot_graphs_moc(M, theta, x, y, x_nozzle, y_nozzle, ...
                             P_nozzlethroat, P_static_wall, n)

%need x_prime, y_prime, u, v, x, y
u = M.*cosd(theta);        %left in Mach x-y components, not veloc
v = M.*sind(theta);

%x_prime and y_prime will be x, y 
%which will be spit out by iterative_solver

% 1 -- Mach number variation
figure
C = min( M, 5 );
surfc(x,y,M, C, 'EdgeAlpha', 0.1)
axis( [0, 0.5, 0, 0.14] );
caxis( [1,6] );
%colormap( hsv );
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
surfc( x,y,P_nozzlethroat, 'EdgeAlpha', 0.1)
%contour( x, y, M )
axis( [0, 0.5, 0, 0.14] );
view(2)
colorbar; 
%colormap( hsv );
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
plot( x, y, 'c.' );
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

end