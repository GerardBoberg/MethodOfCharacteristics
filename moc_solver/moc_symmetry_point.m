function [ x3, y3, slope3, Mach3 ] = moc_symmetry_point( data_1 )
%MOC_INTERIOR_POINT Summary of this function goes here
%   Detailed explanation goes here
global gamma;

% Assume -- Data_1 is above
%
%    1   o
%         \
%          o  3

%% Extract data from the inputs
x1     = data_1(1);
y1     = data_1(2);
slope1 = data_1(3); % slope of streamline
Mach1  = data_1(4);

% Find prandtl-meyer angle (nu) and mach angle (mu)
[ ~, nu1, mu1 ]    = flowprandtlmeyer( gamma, Mach1, 'mach' );

%% Slove for angles of point 3
% Because point 3 is on the symmetry line, the symmetry line must be a
% streamline in the flow. Theta is the slope of the streamline. Therefore
slope3 = 0;

% Now, do prandtl-meyer expansion for point 3
nu_13 = ( slope1 - slope3 ) + nu1;
[ Mach3, ~, mu3] = flowprandtlmeyer( gamma, nu_13, 'nu' );

%% Slove for the location of point 3
% Find the slope of the characteristic lines 
slope_1_dn = slope1 - mu1;
slope_3_dn = slope3 - mu3;

char_13 = tand( ( slope_1_dn + slope_3_dn ) / 2 );

% Point 3 is along the symmetry line, which occurs at:
y3 = 0;
x3 = x1 + ( (1/char_13) * (y3 - y1) );
end
