function [ x3, y3, u3, v3, a3 ] = moc_wall_backsolve( data_1, data_2,...
                                       ~, f_wall_der, x_star, y_star )
%MOC_WALL_POINT Summary of this function goes here
%   Detailed explanation goes here
global a0;
global gamma;
%% Setup initial itteration
tol      = 1e-4; % Tollerance to stop itterating after
max_runs = 5;  % If it takes this long, something's wrong
% extract data from data_1 and data_2
x_1 = data_1( 1 );
y_1 = data_1( 2 );
u_1 = data_1( 3 );
v_1 = data_1( 4 );
a_1 = data_1( 5 );

x_2 = data_2( 1 );
y_2 = data_2( 2 );
u_2 = data_2( 3 );
v_2 = data_2( 4 );
a_2 = data_2( 5 );

% to make sure we don't stop on the first round
x_prev = 100*ones( 6, 1 );
x_next = x_prev * 0;

% initial conditions
x3 = x_star;
y3 = y_star;
u3 = u_2;
v3 = v_2;

y_a = (y_1 + y_2) / 2;
u_a = (u_1 + u_2) / 2;
v_a = (v_1 + v_2) / 2;

u_12 = (u_1 + u_2) / 2;
v_12 = (v_1 + v_2) / 2;

y_12 = (y_1 + y_2) / 2;
a_12 = sqrt( a0^2 - (gamma-1)/2 * ( u_12^2 + v_12^2 ) );

lambda_2_1 = lambda( u_1, v_1, a_1, -1, u_12, a_12 );
lambda_2_2 = lambda( u_2, v_2, a_2, -1, u_12, a_12 );

lambda_12 = ( lambda_2_1 + lambda_2_2 ) / 2;

S_12 = ( a_12^2 * v_12 ) / y_12;
Q_12 = u_12^2 - a_12^2;
R_12 = 2 * u_12 * v_12 - Q_12 * lambda_12;



%% Main loop
not_conv = true;
counter = 1;
while( not_conv )
    
    %% Setup variables that rely on initial conditions
    u_a3 = ( u_a + u3 ) / 2;
    v_a3 = ( v_a + v3 ) / 2;
    y_a3 = ( y_a + y3 ) / 2;
    a_a = sqrt( a0^2 - (gamma-1)/2 * ( u_a^2 + v_a^2 ) );
    a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );
    a_a3 = sqrt( a0^2 - (gamma-1)/2 * (u_a3^2 + v_a3^2) );
    
    lambda_1_a = lambda( u_a, v_a, a_a, +1, u_a3, a_a3 );
    lambda_1_3 = lambda( u3, v3, a3, +1, u_a3, a_a3 );
    lambda_a3 = ( lambda_1_a + lambda_1_3 ) / 2;


    S_a3 = ( a_a3^2 * v_a3 ) / y_a3;
    Q_a3 = u_a3^2 - a_a3^2;
    R_a3 = 2 * u_a3 * v_a3 - Q_a3 * lambda_a3;
    
    theta = atan( v3/u3 );
    

    %% Solve for the next itteration 
    A = [...
            lambda_12, -1,     0,           0,     0,     0;...
            lambda_a3, -1,     0,           0,     0,     0;...
                -S_12,  0,  Q_12,        R_12,     0,     0;...
                -S_a3,  0,  Q_a3,        R_a3, -Q_a3, -R_a3;...
            v_2 - v_1,  0,     0, -(x_2 - x_1),     0,     0;...
                    0,  0,     0,   0,  sin( theta ), -cos( theta ) ...
        ];

    B = [...
            lambda_12 * x_1 - y_1;...
            lambda_a3 * x3 - y3;...
            -S_12 * x_1 + Q_12 * u_1 + R_12 * v_1;...
            -S_a3 * x3;...
            (v_2-v_1)*x_1 - (x_2-x_1)*v_1;...
            0 ...
        ];
        
    x_next = A \ B;
    
    %% Setup variables that don't rely on initial condition.
    y_a = x_next( 2 );
    u_a = x_next( 3 );
    v_a = x_next( 4 );
    u3 = x_next( 5 );
    v3 = x_next( 6 );
      
      
    % Check for convergance
    if( max( abs( x_next - x_prev ) ) < tol )
        not_conv = false;
    end
    x_prev = x_next
    
    
    % Check for infinite loop
    counter = counter + 1;
    if( counter > max_runs )
        error( 'ERROR:MOC:FAILED_TO_CONVERGE', 'counter exceeded max_runs' );
    end
    
end % End loop



%% Return the values we converged to
x3 = x_star;
y3 = y_star;
u3 = x_next( 5 );
v3 = x_next( 6 );
a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );


end

