function [ x3, y3, u3, v3, a3 ] = moc_symmetry_point( data_2 )
%MOC_SYMETRY_POINT Summary of this function goes here
%   Detailed explanation goes here

%% Setup initial itteration
tol      = 1e-4; % Tollerance to stop itterating after
max_runs = 500;  % If it takes this long, something's wrong

% extract data from data_2
x_2 = data_2( 1 );
y_2 = data_2( 2 );
u_2 = data_2( 3 );
v_2 = data_2( 4 );
a_2 = data_2( 5 );

% to make sure we don't stop on the first round
x_prev = 100*ones( 4, 1 );
x_next = x_prev * 0;

u_23 = u_2;
v_23 = v_2;

y_23 = y_2;

lambda_2_2 = ( u_2 * v_2 - a_2 * sqrt( u_2^2 + v_2^2 - a_2^2 ) ) ...
          / ( u_2^2 - a_2^2 );

lambda_23 = lambda_2_2;



%% Main loop
not_conv = true;
counter = 1;
while( not_conv )
    
    %% Setup variables that rely on initial conditions
    a_23 = sqrt( a0^2 - (gamma-1)/2 * (u_23^2 + v_23^2) );

    Q_23 = u_23^2 - a_23;

    R_23 = 2 * u_23 * v_23 - Q_23 * lambda_23;

    S_23 = (a_23^2 * v_23) / y_23;

    %% Solve for the next itteration 
    A = [...
            lambda_23, -1,    0,    0;...
                    0,  1,    0,    0;...
                -S_23,  0, Q_23, R_23;...
                    0,  0,    0,    1 ...
        ];

    B = [...
            lambda_23 * x_2 - y_2;...
                                0;...
            -S_23 * x_2 + Q_23 * u_2 + R_23 * v_2;...
                                0 ...
        ];
        
    x_next = A \ B;
    
    %% Setup variables that don't rely on initial condition.
    u3 = x_next( 3 );
    v3 = 0;
    a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );

    u_23 = ( u_2 + u3 ) / 2;
    v_23 = ( v_2 + v3 ) / 2;

    y_23 = ( y_2 + v3 ) / 2;
    
    
    lambda_2_3 = (u3*v3 - a3 * sqrt(u3^2 + v3^2 - a3^2) ) / (u3^2-a3^2);

    lambda_23 = ( lambda_2_2 + lambda_2_3 ) / 2;
      
      
    % Check for convergance
    if( max( abs( x_next - x_prev ) ) < tol )
        not_conv = false;
    end
    x_prev = x_next;
    
    
    % Check for infinite loop
    counter = counter + 1;
    if( counter > max_runs )
        error( 'ERROR:MOC:FAILED_TO_CONVERGE', 'counter exceeded max_runs' );
    end
end % End loop



%% Return the values we converged to
x3 = x_next( 1 );
y3 = 0;
u3 = x_next( 3 );
v3 = 0;
a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );


end

