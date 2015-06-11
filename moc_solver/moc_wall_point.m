function [ x3, y3, slope3, Mach3 ] = moc_wall_point( data_1,...
                                                f_wall, f_wall_der, x_star )
%MOC_WALL_POINT Summary of this function goes here
%   Detailed explanation goes here
global gamma;
% Assume -- Data_1 is below
%              -------
%        ---o--   3
%   -----  /
%---      /
%    1   o

%% Extract data from the inputs
x1     = data_1(1);
y1     = data_1(2);
slope1 = data_1(3); % slope of streamline
Mach1  = data_1(4);

% Find prandtl-meyer angle (nu) and mach angle (mu)
[ ~, nu1, mu1 ]    = flowprandtlmeyer( gamma, Mach1, 'mach' );

%% Use fsolve to determine where x3 lives
f = @(x)( x - find_expected_x( x, x1, y1, slope1, nu1, mu1,...
                                                  f_wall, f_wall_der ) );

dy = f_wall( x1 ) - y1; % dist point 1 <--> wall same order of magnitude as 


options = optimset('Display', 'off');
x3 = fsolve( f, x1+dy, options );

if( x3 >= x_star )
    error( 'ERROR:MOC:MISSED_WALL',...
           'wall boundary has eneded, exceeded maximum x value' );
end

%% Now pull out all of the information about point 3
y3     = f_wall( x3 );
slope3 = atand( f_wall_der( x3 ) );

nu_13  = ( slope3 - slope1 ) + nu1;
[ Mach3, ~, ~ ] = flowprandtlmeyer( gamma, nu_13, 'nu' );


end

%% Private function used in fsolve in the main solver
function xr = find_expected_x( xg, x1, y1, slope1, nu1, mu1, f_wall, f_wall_der )
global gamma;
slope3 = atand( f_wall_der( xg ) );
nu_13  = ( slope3 - slope1 ) + nu1;
[ ~, ~, mu3 ] = flowprandtlmeyer( gamma, nu_13, 'nu' );
char_13 = tand( ( (slope1 + mu1) + (slope3 + mu3) ) /2 );

F  = @(x)( (y1 + char_13*(x-x1)) - f_wall(x) );


options = optimset('Display', 'off');
xr = fsolve( F, xg, options );
end

