function [ x3, y3, slope3, Mach3 ] = moc_interior_point( data_1, data_2 )
%MOC_INTERIOR_POINT Summary of this function goes here
%   Detailed explanation goes here
global gamma;

% Assume -- Data_1 is above, Data_2 is below
%
%    1   o
%         \
%          o  3
%         /
%        /
%    2  o

%% Extract data from the inputs
x1     = data_1(1);
y1     = data_1(2);
slope1 = data_1(3); % slope of streamline, in degrees
Mach1  = data_1(4);

x2     = data_2(1);
y2     = data_2(2);
slope2 = data_2(3);
Mach2  = data_2(4);


% Find prandtl-meyer angle (nu) and mach angle (mu)
[ ~, nu1, mu1 ]    = flowprandtlmeyer( gamma, Mach1, 'mach' );
[ ~, nu2, mu2 ]    = flowprandtlmeyer( gamma, Mach2, 'mach' );

%% Slove for angles of point 3
% First, solve for the theta value of point 3, using prandtl-meyer angles
slope3 = ( (slope1 + nu1) + (slope2 - nu2) )/ 2;

% Now, do prandtl-meyer expansion for point 3
nu_13 = ( slope1 - slope3 ) + nu1;
nu_23 = ( slope2 - slope3 ) + nu2;
nu_avg = ( nu_13 + nu_23) / 2; % 1 and 2 should be = ; average to be safe.
[ Mach3, ~, mu3] = flowprandtlmeyer( gamma, nu_avg, 'nu' );

%% Slove for the location of point 3
% Find the slope of the characteristic lines 
char_1_dn = ( slope1 - mu1 );
char_2_up = ( slope2 + mu2 );

char_3_dn = ( slope3 - mu3 );
char_3_up = ( slope3 + mu3 );

char_13 = tand( ( char_1_dn + char_3_dn ) / 2 );
char_23 = tand( ( char_2_up + char_3_up ) / 2 );

% Use the slopes to solve for x and y locations
epsilon = x1 - x2;
dx = ( (y1 - y2) - (char_23 * epsilon) ) / ( char_23 - char_13 );

x3 = x1 + dx;
y3 = y1 + ( dx * char_13);
end
