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
[ ~, ~, mu1 ]    = flowprandtlmeyer( gamma, Mach1, 'mach' );

%% Use fsolve to determine where x3 lives
f = @(x)( find_expected_x( x,...
                            x1, y1, slope1, Mach1, mu1,...
                            f_wall, f_wall_der ) );

dy = f_wall( x1 ) - y1; % dy same order of magnitude as dx
x0 = [  x1 + dy;  f_wall( x1+dy ); Mach1 * cosd( slope1 ); Mach1 * sind( slope1 )]; 


options = optimset('Display', 'off');
data3 = fsolve( f, x0, options );

%% Extract information from data3
x3 = data3( 1 );

if( x3 >= x_star )
    error( 'ERROR:MOC:MISSED_WALL',...
           'wall boundary has eneded, exceeded maximum x value' );
end % error gets caught in moc_iterative_solver

y3 = f_wall( x3 );

M_cos_t = data3( 3 );
M_sin_t = data3( 4 );

slope3  = atand( M_sin_t / M_cos_t );
Mach3   = M_cos_t / cosd( slope3 );
end

%% Private function used in fsolve in the main solver
function er = find_expected_x( data3,...
                                    x1, y1, slope1, M1, mu1,...
                                    f_wall, f_wall_der )
global gamma;

% Extract data from input
x3 = data3( 1 );
y3 = data3( 2 );
M_cos_t = data3( 3 );
M_sin_t = data3( 4 );
alpha = f_wall_der( x3 );

theta3  = atand( M_sin_t / M_cos_t );
M3      = M_cos_t / cosd( theta3 );
[ ~, ~, mu3 ] = flowprandtlmeyer( gamma, M3, 'mach' );

% Calculate the characteristic slopes
char_1_up = tand( slope1 + mu1 );
char_3_up = tand( theta3 + mu3 );

char_13 = (char_1_up + char_3_up ) / 2;

% calculate the parameters s, q, r
theta_13 = (slope1 + theta3) / 2;

M13 = (M1+M3) / 2;

a13_over_a3 = (1 + ((gamma-1)/2)*M3^2) / (1+((gamma-1)/2)*M13^2);

s_13 = (a13_over_a3) * M13 * sind( theta_13 ) / ((y1+y3)/2);

q_13 = M13^2 * cosd( theta_13 ) - 1;

r_13 = 2 * M13^2 * sind( theta_13 ) * cosd( theta_13 ) - q_13 * char_13;

% calculate 
% setup the matrix
A = [
        char_13,    -1, 0,  0;...
            0,       1, 0,  0;...
        -s_13,       0,  q_13,   r_13;...
            0,       0,  alpha,   -1;
    ];

b = [
        (char_13 * x1) - y1;...
        f_wall( x3 );...
        -(s_13 * x1) + (q_13 * M1 * cosd( slope1 )) + (r_13 * M1 * sind( slope1 ));...
        0
    ];


% Calculate the error value
er = ( A * data3 ) - b;
end

