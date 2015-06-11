function [ x, y, slope, M ] = moc_iterative_solver( ...
                                  x_nozzle, y_nozzle, n, thermo, y_throat )
%MOC_ITERATIVE_SOLVER Summary of this function goes here
%   Detailed explanation goes here

N_POLY = 4; % degree of polynomial for fitting to nozzle x & y

%% Setup variables
% Setup matricies
dimensions = [ n, 3 * n ];

x = zeros( dimensions );
y = zeros( dimensions );
slope = zeros( dimensions ); % slope of the streamline in degrees (i,j)
Mach = zeros( dimensions );

% Setup initial data line
global gamma;
gamma = thermo.gamma;
[ x(:,1), y(:,1), slope(:,1), Mach(:,1) ] = initial_data_line...
                                     ( y_throat, y_throat/2, gamma, n );
figure();
plot( x(:,1), y(:,1), 'bx' );
                                 
% create wall gemoetry function handles
wall_poly     = polyfit( x_nozzle, y_nozzle, N_POLY );
wall_poly_der = polyder( wall_poly );

f_wall     = @(x)( polyval( wall_poly, x ) );
f_wall_der = @(x)( polyval( wall_poly_der, x ) ); 

% define the final points along the fixed wall
x_star = x_nozzle(end);
y_star = y_nozzle(end);

%% Main Loop
% indexing:
% i3 [ 3  ~  ~  ~  ~ ]
% i2 [ 2  3  3  ~  ~ ]
% i1 [ 1  2  2  3  3 ]
%     j1 j2 j3 j4 j5

figure();
hold on
plot(x_nozzle,y_nozzle,'LineWidth',4,'Color','k')


not_done  = true;
% n = number of given initial characteristic lines
char_line = 0; 
while( not_done )    % For each characteristic line, we don't know how many
    char_line = char_line + 1;
    
    % ----- This is a two-part loop -----
    % part 1: setup indexing. Handle upper wall boundary
    % part 2: cast down to the symmetry plane
    
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
        ii = n; % ii ~ height
        jj = 1 + 2 * ( char_line - n );
        
        % data 1 is below, data 2 is above ( and a step behind )
        data_1 = [ x(ii-1, jj-1), y(ii-1, jj-1),...
                   slope(ii-1, jj-1), Mach(ii-1, jj-1) ];
        data_2 = [ x(ii  , jj-2), y(ii  , jj-2),...
                   slope(ii  , jj-2), Mach(ii  , jj-2) ];
        
        try
            % Cast C -- upper wall boundary
            [ x(ii,jj), y(ii,jj), slope(ii,jj), Mach(ii,jj) ] = ...
                moc_wall_point( data_1, f_wall, f_wall_der, x_star );
            %plot( [], [], 'r.'
        catch e % if C goes past x_star, it throws an exception
            if( strcmp( e.identifier, 'ERROR:MOC:MISSED_WALL' ) )
                % Cast D -- backsolve wall boundary
                [ x(ii,jj), y(ii,jj), slope(ii,jj), Mach(ii,jj) ] = ...
                    moc_wall_backsolve( data_2, data_1,... % flip intended
                                        f_wall, f_wall_der,...
                                        x_star, y_star );
                not_done = false;
% exit flag is right here ^
%   occurs when the cast for a char line onto the wall misses the wall, for
%   char line > n, the number of initial data points
            else
                throw( e );
            end
        end
        %plot( real(x(ii,jj)), real(y(ii,jj)), 'rx' );
    end % END SETUP INDEX
    
    %% CAST THE LINE DOWN TO THE SYMMETRY PLANE
    % ii and jj are setup for the initial data point. 
    % Follow a char_line down until it hits the symmetry plane.
    % 
    % The indexing works as follows:
    % i3 [ 3  ~  4  ~  5  ~  ~  ~ ]
    % i2 [ 2  3  3  4  4  5  5  ~ ]
    % i1 [ 1  2  2  3  3  4  4  5 ]
    %     j1 j2 j3 j4 j5
    i_values = floor( ii:-0.5:1 ); % traces down the line
    j_values = jj:jj+(2*(ii-1));
    final_value = length( i_values );
    
    for index = 2:final_value % for every point on the line until symmetry
        % setup indexing
        ic = i_values( index ); % ic, jc are the current index
        jc = j_values( index );
        
        i1 = i_values( index-1 ); % i1, j1 are the previous, above index
        j1 = j_values( index-1 );
        
        i2 = i1 - 1; % i2, j2 are the previous, below index
        j2 = j1;
        
        % data_1 is above, data_2 is below
        %fprintf( ['---char_line = ', num2str( char_line ), ' ---\n' ] );
        data_1 = [ x(i1,j1), y(i1,j1), slope(i1,j1), Mach(i1,j1) ];
        
        if( index == final_value ) % the last point on the line
            % Cast B -- symmetry line
            [ x(ic,jc), y(ic,jc), slope(ic,jc), Mach(ic,jc) ] = ...
                moc_symmetry_point( data_1 );
            plot( [x(ii,jj), x(ic,jc)], [y(ii,jj), y(ic,jc)], 'b' );
            
        else % every other interior point
            data_2 = [ x(i2, j2), y(i2,j2), slope(i2,j2), Mach(i2,j2) ];
            
            % Cast A -- interior point
            [ x(ic,jc), y(ic,jc), slope(ic,jc), Mach(ic,jc) ] = ...
                moc_interior_point( data_1, data_2 );
            plot( x(ic,jc), y(ic,jc), 'r*' );
        end
        
    end
    
end

%% Cleanup at the end

% Lerp mach numbers to fill in the gaps on the top surface. 
% Makes everything look nicer.
for kk = 2:size( Mach, 2 )
    % if we're inbetween two points that have defined mach values
    if( ~(Mach(end,kk) > 0 ) && (Mach(end, kk-1) > 0) && (Mach(end,kk+1) > 0) )
        display( [ 'interp filling top at kk = ', num2str( kk ) ] );
        
        % assume linear variation, and fill in the values
        x( end, kk ) = (x(end,kk-1) + x(end,kk+1))/2;
        y( end, kk ) = f_wall( x(end,kk) );
        Mach( end, kk ) = (Mach(end,kk-1) + Mach(end,kk+1))/2;
        slope( end,kk ) = atand( f_wall_der( x(end,kk) ) );
    end
end

%
% Smooth out the points along the symmetry line
L = size( x, 2 );
x = [ zeros( 1, L ); x ];
y = [ zeros( 1, L ); y ];
slope = [ zeros( 1, L ); slope ];
M = [ zeros( 1, L ); Mach ];

% NaN doesn't graph, prevents crazy behaviour around zero
for ii = 1:size( x, 1 ) 
    for jj = 2:size( x, 2 )
        if( (x( ii, jj ) == 0) || ( M( ii, jj ) == 0 ) )
            x( ii, jj ) = NaN;
            y( ii, jj ) = NaN;
            slope( ii, jj ) = NaN;
            M( ii, jj ) = NaN;
        end
    end
end

% Fill in spaces along the symmetry line
for kk = 2:L
    if( (y(2,kk) > 0 ) &&...
        (y(2, kk-1) == 0) && (y(2,kk+1) == 0) )
        x( 1, kk ) = (x(2,kk-1) + x(2,kk+1))/2;
        y( 1, kk ) = 0;
        M( 1, kk ) = (M(2,kk-1) + M(2,kk+1))/2;
        
        if( ~(x( 1, kk-1 )  > 0) )
            x( 1, kk-1 ) = x( 1, kk );
            M( 1, kk-1 ) = M( 1, kk );
            y( 1, kk-1 ) = y( 1, kk );
        end
    end
end
end
