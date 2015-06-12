function [ x_p, y_p, t_p, M_p ] = moc_clean_top_row( x, y, t, M )
%MOC_CLEAN_TOP_ROW Summary of this function goes here
%   Detailed explanation goes here

L = size( x, 2 );

% pre-allocate
x_p = [ x(1,:); x; x(end,:) ];
y_p = [ y(1,:); y; y(end,:) ];
t_p = [ t(1,:); t; t(end,:) ];
M_p = [ M(1,:); M; M(end,:) ];

for ii = 2:L % too lazy to make this a MATLAB matrix operation
    x_p( 1, ii ) = ( x_p( 2, ii-1 ) + x_p( 2, ii ) ) / 2;
    y_p( 1, ii ) = ( y_p( 2, ii-1 ) + y_p( 2, ii ) ) / 2;
    t_p( 1, ii ) = ( t_p( 2, ii-1 ) + t_p( 2, ii ) ) / 2;
    M_p( 1, ii ) = ( M_p( 2, ii-1 ) + M_p( 2, ii ) ) / 2;
    
    
    x_p( end, ii ) = ( x_p( end-1, ii-1 ) + x_p( end-1, ii ) ) / 2;
    y_p( end, ii ) = ( y_p( end-1, ii-1 ) + y_p( end-1, ii ) ) / 2;
    t_p( end, ii ) = ( t_p( end-1, ii-1 ) + t_p( end-1, ii ) ) / 2;
    M_p( end, ii ) = ( M_p( end-1, ii-1 ) + M_p( end-1, ii ) ) / 2;
end

figure();
plot( x_p(end,:), y_p(end,: ) );

end

