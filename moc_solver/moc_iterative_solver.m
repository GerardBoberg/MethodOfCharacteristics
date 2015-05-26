function [ x, y, u, v, a, M ] = moc_iterative_solver( ...
                                        x_nozzle, y_nozzle, n, thermo )
%MOC_ITERATIVE_SOLVER Summary of this function goes here
%   Detailed explanation goes here

N_POLY = 6;

%% Setup variables
% Setup matricies
dimensions = [ n, 10 * n ];

x = zeros( dimensions );
y = zeros( dimensions );
u = zeros( dimensions );
v = zeros( dimensions );
a = zeros( dimensions );

% Setup initial data line
global gamma;
gamma = thermo.gamma;
R     = thermo.R;
T0    = thermo.T0;
[ x(:,1), y(:,1), u(:,1), v(:,1), a(:,1) ] = initial_data_line...
                                     ( y_throat, y_throat/2, gamma, R, T0, n );
global a0;
a0 = sqrt( gamma * R * T0 );
% wall gemoetry and function handles
wall_poly     = polyfit( x_nozzle, y_nozzle, N_POLY );
wall_poly_der = polyder( wall_poly );

f_wall     = @(x)( polyval( wall_poly, x ) );
f_wall_der = @(x)( polyval( wall_poly_der, x ) ); 

x_star = x_nozzle(end);
y_star = y_nozzle(end);

%% Main Loop
% indexing:
% i3 [ 3  ~  ~  ~  ~ ]
% i2 [ 2  3  3  ~  ~ ]
% i1 [ 1  2  2  3  3 ]
%     j1 j2 j3 j4 j5

not_done  = true;
% n = number of given initial characteristic lines
char_line = 0; 
while( not_done )    % For each characteristic line, we don't know how many
    char_line = char_line + 1;
    
    %% SETUP INDEXING
    % Setup the indexing for the starting point of the char line
    %   the top-left most point that we will cast down to the syemmtry
    if( char_line <= n ) % char_line inside 1:n
        % We're still on the initial data line. Just setup index.
        ii = char_line;
        jj = 1;
    else % if char line > n
        % We have left the initial data line. The intial point needs to be
        %   cast onto the wall in order to continue with index setup.
        ii = n;
        jj = 1 + 2 * ( char_line - n );
        
        % data 1 is below, data 2 is above ( and a step behind )
        data_1 = [ x(ii-1, jj-1), y(ii-1, jj-1),...
                   u(ii-1, jj-1), v(ii-1, jj-1), a(ii-1 ,jj-1) ];
        data_2 = [ x(ii  , jj-2), y(ii  , jj-2),...
                   u(ii  , jj-2), v(ii  , jj-2), a(ii   ,jj-2) ];
        
        try
            % Cast C
            [ x(ii,jj), y(ii,jj), u(ii,jj), v(ii,jj), a(ii,jj) ] = ...
                moc_wall_point( data_1, f_wall, f_wall_der, x_star );
        catch e
            if( strcmp( e.identifier, 'ERROR:MOC:MISSED_WALL' ) )
                % Cast D
                [ x(ii,jj), y(ii,jj), u(ii,jj), v(ii,jj), a(ii,jj) ] = ...
                    moc_wall_backsolve( ...
                        data_1, data_2, f_wall, f_wall_der, x_star, y_star );
                not_done = false;
% exit flag is right here ^
%   occurs when the cast for a char line onto the wall misses the wall, for
%   char line > n, the number of initial data points
            else
                throw( e );
            end
        end
    end % END SETUP INDEX
    
    %% CAST THE LINE DOWN TO THE SYMMETRY PLANE
    % ii and jj are setup for the initial data point. 
    % Follow a char_line down until it hits the symmetry plane.
    % 
    % The indexing works as follows:
    % i3 [ 3  ~  ~  ~  ~ ]
    % i2 [ 2  3  3  ~  ~ ]
    % i1 [ 1  2  2  3  3 ]
    %     j1 j2 j3 j4 j5
    i_values = floor( ii:-0.5:1 ); % traces down the line
    j_values = jj:jj+(2*(ii-1));
    final_value = length( i_values );
    
    for index = 2:final_value % for every point on the line until symmetry
        
        % setup indexing
        ic = i_values( index ); % ic, jc are the current index
        jc = j_values( index );
        
        i2 = i_values( index-1 ); % i2, j2 are the previous, above index
        j2 = j_values( index-1 );
        
        i1 = i2 - 1; % i1, j1 are the previous, below index
        j1 = j2;
        
        % data_2 is above, data_1 is below
        data_2 = [ x(i2, j2), y(i2,j2), u(i2,j2), v(i2,j2), a(i2,j2) ];
        
        if( index == final_value ) % the last point on the line
            % Cast B
            [ x(ic,jc), y(ic,jc), u(ic,jc), v(ic,jc), a(ic,jc) ] = ...
                moc_symmetry_point( data_2 );
            
        else % every other interrior point
            data_1 = [ x(i1, j1), y(i1,j1), u(i1,j1), v(i1,j1), a(i1,j1) ];
            
            % Cast A
            [ x(ic,jc), y(ic,jc), u(ic,jc), v(ic,jc), a(ic,jc) ] = ...
                moc_interior_point( data_1, data_2 );
        end
    end
    
end


% cleanup at the end
M = (u.^2 + v.^2) ./ a.^2;
end
