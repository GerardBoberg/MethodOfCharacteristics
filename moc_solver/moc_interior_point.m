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
slope1 = data_1(3); % slope of streamline
Mach1  = data_1(4);


x2     = data_2(1);
y2     = data_2(2);
slope2 = data_2(3); % slope of streamline
Mach2  = data_2(4);

% Find prandtl-meyer angle (nu) and mach angle (mu)
[ ~, ~, mu1 ]    = flowprandtlmeyer( gamma, Mach1, 'mach' );
[ ~, ~, mu2 ]    = flowprandtlmeyer( gamma, Mach2, 'mach' );

%% Use fsolve to determine where x3 lives
f = @(x) ( find_expected_x( x, ...
                               x1, y1, slope1, Mach1, mu1,...
                               x2, y2, slope2, Mach2, mu2 ) );

dy = abs( y2 - y1 ); % dy order of magnitude as dx
theta_12 = (slope1 + slope2) / 2;
M_12     = (Mach1 + Mach2) / 2;

x0 = [  x1 + dy;  (y1+y2) / 2; M_12 * cosd(theta_12); M_12 * sind(theta_12)]; 

options = optimset('Display', 'off');
data3 = fsolve( f, x0, options );

%% Now pull out all of the information about point 3
x3 = data3( 1 );
y3 = data3( 2 );

M_cos_t = data3( 3 );
M_sin_t = data3( 4 );

slope3  = atand( M_sin_t / M_cos_t );
Mach3      = M_cos_t / cosd( slope3 );
end

%% Private function used in fsolve in the main solver
function er = find_expected_x( data3, ...
                                      x1, y1, slope1, M1, mu1,...
                                      x2, y2, slope2, M2, mu2 )
global gamma;

% Extract data from input
%x3 = data3( 1 );
y3 = data3( 2 );
M_cos_t = data3( 3 );
M_sin_t = data3( 4 );

theta3  = atand( M_sin_t / M_cos_t );
M3      = M_cos_t / cosd( theta3 );
[ ~, ~, mu3 ] = flowprandtlmeyer( gamma, M3, 'mach' );

% Calculate the characteristic slopes
char_1_dn = tand( slope1 - mu1 );
char_2_up = tand( slope2 + mu2 );

char_3_dn = tand( theta3 - mu3 );
char_3_up = tand( theta3 + mu3 );

char_13 = (char_1_dn + char_3_dn ) / 2;
char_23 = (char_2_up + char_3_up ) / 2;

% calculate the parameters s, q, r
theta_13 = (slope1 + theta3) / 2;
theta_23 = (slope2 + theta3) / 2;

M13 = (M1+M3) / 2;
M23 = (M2+M3) / 2;

a13_over_a3 = (1 + ((gamma-1)/2)*M3^2) / (1+((gamma-1)/2)*M13^2);
a23_over_a3 = (1 + ((gamma-1)/2)*M3^2) / (1+((gamma-1)/2)*M23^2);

s_13 = (a13_over_a3) * M13 * sind( theta_13 ) / ((y1+y3)/2);
s_23 = (a23_over_a3) * M23 * sind( theta_23 ) / ((y2+y3)/2);

q_13 = M13^2 * cosd( theta_13 ) - 1;
q_23 = M23^2 * cosd( theta_23 ) - 1;

r_13 = 2 * M13^2 * sind( theta_13 ) * cosd( theta_13 ) - q_13 * char_13;
r_23 = 2 * M23^2 * sind( theta_23 ) * cosd( theta_23 ) - q_23 * char_23;

% calculate 
% setup the matrix
A = [
        char_13,   -1, 0,  0;...
        char_23,   -1, 0,  0;...
        -s_13,  0,  q_13,   r_13;...
        -s_23,  0,  q_23,   r_23;
    ];

b = [
        (char_13 * x1) - y1;...
        (char_23 * x2) - y2;...
        -(s_13 * x1) + (q_13 * M1 * cosd( slope1 )) + (r_13 * M1 * sind( slope1 ));...
        -(s_23 * x2) + (q_23 * M2 * cosd( slope2 )) + (r_23 * M2 * sind( slope2 ));
    ];


% Calculate the error value
er = ( A * data3 ) - b;
end