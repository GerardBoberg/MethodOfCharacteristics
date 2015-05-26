function [ x, y, u, v, a ] = initial_data_line( yt, rt, g, R, T0, n )
%INITIALDATALINE Summary of this function goes here
%   Detailed explanation goes here

alpha   = sqrt( 2 / (rt * yt * (g+1) ) );
epsilon = yt / 8 * sqrt( (g+1) * 2 / (rt / yt) );

y = linspace( 0, yt, n );
x_bar = -alpha * y.^2 * (g+1) / 8;
x = x_bar + epsilon;

a_star = sqrt( g * R * 0.8333*T0 );
a0     = sqrt( g * R * T0 );

u_prime = (alpha * x_bar) + (alpha^2 * y.^2 * (g + 1) / 4);
u = a_star * ( 1 + u_prime );
a = sqrt( a0^2 - (g-1)/2*(u.^2) );

v = zeros( size( y ) );
end

