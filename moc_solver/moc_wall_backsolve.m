function [ x3, y3, slope3, Mach3 ] = moc_wall_backsolve( data_1, data_2,...
                                               f_wall, f_wall_der,...
                                               x_star, y_star )
%MOC_WALL_POINT Summary of this function goes here
%   Detailed explanation goes here
global gamma;

% Assume -- Data_1 is above, Data_2 is below, A is between 1 and 2
%              -------
%        ---o--   3
% 1 ---o-  /  
%---    \ /  
%   A    o 
%         \
%   2       o

%% Extract data from the inputs
x1     = data_1(1);
y1     = data_1(2);
slope1 = data_1(3); % slope of streamline, in degrees
Mach1  = data_1(4);

x2     = data_2(1);
y2     = data_2(2);
slope2 = data_2(3);
Mach2  = data_2(4);

% We know where point 3 is, and the slope of the streamline.
x3     = x_star;
y3     = y_star;
slope3 = atand( f_wall_der( x3 ) );
% only need to find the Mach, but relies on point a between 1 and 2

%% Use fsolve to determine where point a lives. 
f = @(ag)( x3 - find_expected_x( ag,...
                                     x1, y1, slope1, Mach1,...
                                     x2, y2, slope2, Mach2,...
                                     f_wall, x3, slope3 ) );
ag = fsolve( f, 0.5 );

%% Determine the mach number based off our located point A
% assume lin variation
slopea = slope1 + ( slope2 - slope1 ) * ag; 
Macha  = Mach1 + ( Mach2 - Mach1 ) * ag;
[ ~, nua, ~ ] = flowprandtlmeyer( gamma, Macha, 'mach' );

% Use to solve for Mach at point 3
nu_a3  = ( slope3 - slopea ) + nua;
[ Mach3, ~, ~ ] = flowprandtlmeyer( gamma, nu_a3, 'nu' );
end


%% Private function used in fsolve in the main solver
function xr = find_expected_x( ag,...
                                     x1, y1, slope1, Mach1,...
                                     x2, y2, slope2, Mach2,...
                                     f_wall, x3, slope3 )
global gamma;

% ag is a percentage between 0 and 1. Force inside those bounds.
if     ( ag <= 0 )
    ag = 0;
elseif ( ag >= 1 )
    ag = 1;
end

% derive the information of point a. Linear variation between 1 and 2
xa = x1 + ( x2 - x1 ) * ag;
ya = y1 + ( y2 - y1 ) * ag;
slopea = slope1 + ( slope2 - slope1 ) * ag;
Ma = Mach1 + ( Mach2 - Mach1 ) * ag;

[ ~, nua, mua ] = flowprandtlmeyer( gamma, Ma, 'mach' );

% Use the slope information to find the Mach angle at 3
nu_a3  = ( slope3 - slopea ) + nua;
[ ~, ~, mu3 ] = flowprandtlmeyer( gamma, nu_a3, 'nu' );
char_a3 = tand( ( (slopea + mua) + (slope3 + mu3) ) /2 );

F  = @(x)( (ya + char_a3*(x-xa)) - f_wall(x) );
xr = fsolve( F, x3 );
end
