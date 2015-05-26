function [ x3, y3, u3, v3, a3 ] = moc_wall_point( data_1,...
                                                f_wall, f_wall_der, x_star )
%MOC_WALL_POINT Summary of this function goes here
%   Detailed explanation goes here
global a0;
global gamma;
%% Setup initial itteration
tol      = 1e-4; % Tollerance to stop itterating after
max_runs = 50;  % If it takes this long, something's wrong
begin_checking_out_of_bounds = 5; % don't start immediatly

% extract data from data_1 and data_2
x_1 = data_1( 1 );
y_1 = data_1( 2 );
u_1 = data_1( 3 );
v_1 = data_1( 4 );
a_1 = data_1( 5 );

% to make sure we don't stop on the first round
x_prev = 100*ones( 4, 1 );
x_next = x_prev * 0;

% initial conditions
x3 = x_1;
y3 = f_wall( x3 );
alpha = f_wall_der( x3 );

u_13 = u_1;
v_13 = v_1;

y_13 = y_1;

lambda_1_1 = lambda( u_1, v_1, a_1, +1 );
lambda_13 = lambda_1_1;

%% Main loop
not_conv = true;
counter = 1;
while( not_conv )
    
    %% Setup variables that rely on initial conditions
    a_13 = sqrt( a0^2 - (gamma-1)/2 * (u_13^2 + v_13^2) );

    Q_13 = u_13^2 - a_13^2;

    R_13 = 2 * u_13 * v_13 - Q_13 * lambda_13;

    S_13 = (a_13^2 * v_13) / y_13;
    

    %% Solve for the next itteration 
    A = [...
            lambda_13, -1,     0,    0;...
                alpha, -1,     0,    0;...
                -S_13,  0,  Q_13, R_13;...
                    0,  0, alpha,   -1 ...
        ];

    B = [...
            lambda_13 * x_1 - y_1;...
            alpha * x3 - y3;...
            -S_13 * x_1 + Q_13 * u_1 + R_13 * v_1;...
            0 ...
        ];
        
    x_next = A \ B;
    
    %% Setup variables that don't rely on initial condition.
    x3 = x_next( 1 );
    y3 = f_wall( x3 );
    u3 = x_next( 3 );
    v3 = x_next( 4 );
    a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );
    
    u_13 = ( u_1 + u3 ) / 2;
    v_13 = ( v_1 + v3 ) / 2;

    y_13 = ( y_1 + y3 ) / 2;
    
    
    lambda_1_3 = lambda( u3, v3, a3, +1, u_13, a_13 );
    lambda_13 = ( lambda_1_1 + lambda_1_3 ) / 2;
      
      
    % Check for convergance
    if( max( abs( x_next - x_prev ) ) < tol )
        not_conv = false;
    end
    x_prev = x_next;
    
    
    % Check for infinite loop
    counter = counter + 1;
    if( counter > max_runs )
        error( 'ERROR:MOC:FAILED_TO_CONVERGE',...
            'counter exceeded max_runs in moc_wall, C' );
    end
    
    % Check for fail to converge / converging outside bounds
    if( (counter > begin_checking_out_of_bounds) && (x3 > x_star) )
        error( 'ERROR:MOC:MISSED_WALL', 'wall boundary has ended' );
    end
end % End loop



%% Return the values we converged to
x3 = x_next( 1 );
y3 = x_next( 2 );
u3 = x_next( 3 );
v3 = x_next( 4 );
a3 = sqrt( a0^2 - (gamma-1)/2 * ( u3^2 + v3^2 ) );


end

