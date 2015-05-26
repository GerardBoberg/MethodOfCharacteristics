function [ x3, y3, u3, v3, a3 ] = moc_interior_point( data_1, data_2 )
%MOC_INTERIOR_POINT Summary of this function goes here
%   Detailed explanation goes here
global a0;
global gamma;
%% Setup initial itteration
tol      = 1e-4; % Tollerance to stop itterating after
max_runs = 500;  % If it takes this long, something's wrong

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
x_prev = 100*ones( 4, 1 );
x_next = x_prev * 0;


% initial conditions
u_13 = u_1;
v_13 = v_1;
u_12 = (u_1 + u_2)/2;
a_12 = (a_1 + a_2) / 2;

u_23 = u_2;
v_23 = v_2;

if( abs( y_1 ) < tol*10 )
    y_13 = 1;
else
    y_13 = y_1;
end
y_23 = y_2;

lambda_1_1 = lambda( u_1, v_1, a_1, +1, u_12, a_12 );
lambda_2_2 = lambda( u_1, v_1, a_1, -1, u_12, a_12 );

lambda_13 = lambda_1_1;
lambda_23 = lambda_2_2;



%% Main loop
not_conv = true;
counter = 1;
while( not_conv )
    
    %% Setup variables that rely on initial conditions
    a_13 = sqrt( a0^2 - (gamma-1)/2 * (u_13^2 + v_13^2) );
    a_23 = sqrt( a0^2 - (gamma-1)/2 * (u_23^2 + v_23^2) );

    Q_13 = u_13^2 - a_13^2;
    Q_23 = u_23^2 - a_23^2;

    R_13 = (2 * u_13 * v_13) - (Q_13 * lambda_13);
    R_23 = (2 * u_23 * v_23) - (Q_23 * lambda_23);

    if( abs( y_13 ) < tol*10 )
        S_13 = 0;
    else
        S_13 = (a_13^2 * v_13) / y_13;
    end
    S_23 = (a_23^2 * v_23) / y_23;

    %% Solve for the next itteration 
    A = [...
            lambda_13, -1,    0,    0;...
            lambda_23, -1,    0,    0;...
                -S_13,  0, Q_13, R_13;...
                -S_23,  0, Q_23, R_23 ...
        ];

    B = [...
            lambda_13 * x_1 - y_1;...
            lambda_23 * x_2 - y_2;...
            -S_13 * x_1 + Q_13 * u_1 + R_13 * v_1;...
            -S_23 * x_2 + Q_23 * u_2 + R_23 * v_2 ...
        ];
    x_next = A \ B;
    
    %% Setup variables that don't rely on initial condition.
    y3 = x_next( 2 );
    u3 = x_next( 3 );
    v3 = x_next( 4 );
    a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );
    
    u_13 = ( u_1 + u3 ) / 2;
    v_13 = ( v_1 + v3 ) / 2;

    u_23 = ( u_2 + u3 ) / 2;
    v_23 = ( v_2 + v3 ) / 2;

    y_13 = ( y_1 + y3 ) / 2;
    y_23 = ( y_2 + y3 ) / 2;
    
    
    a_13 = sqrt( a0^2 - (gamma-1)/2 * (u_13^2 + v_13^2) );
    a_23 = sqrt( a0^2 - (gamma-1)/2 * (u_23^2 + v_23^2) );
    lambda_1_3 = lambda( u3, v3, a3, +1, u_13, a_13 );
    lambda_2_3 = lambda( u3, v3, a3, -1, u_23, a_23 );

    lambda_13 = ( lambda_1_1 + lambda_1_3 ) / 2;
    lambda_23 = ( lambda_2_2 + lambda_2_3 ) / 2;
      
      
    % Check for convergance
    if( max( abs( x_next - x_prev ) ) < tol )
        not_conv = false;
    end
    x_prev = x_next;
    
    
    % Check for infinite loop
    counter = counter + 1;
    if( counter > max_runs )
        error( 'ERROR:MOC:FAILED_TO_CONVERGE',...
            'counter exceeded max_runs in moc_interior, A' );
    end
end % End loop



%% Return the values we converged to
x3 = x_next( 1 );
y3 = x_next( 2 );
u3 = x_next( 3 );
v3 = x_next( 4 );
a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );

end
