function [ x, y, slope, Mach ] = initial_data_line( yt, rt, g, n )
%INITIALDATALINE Summary of this function goes here
%   Detailed explanation goes here

alpha   = sqrt( 2 / (rt * yt * (g+1) ) );
epsilon = yt / 8 * sqrt( (g+1) * 2 / (rt / yt) );

y = linspace( yt/1000, yt, n );
y = yt * ( y./yt ).^(1/2); % transforms the vector space to be more top-heavy
x_bar = -alpha * y.^2 * (g+1) / 8;
x = x_bar + epsilon;

u_prime = (alpha * x_bar) + (alpha^2 * y.^2 * (g + 1) / 4);
Mach  = ( 1.0001 + u_prime );
slope = zeros( size( Mach ) );

end

