function [ L ] = lambda( u, v, a, sign, u13, a13  )
%LAMBDA Summary of this function goes here
%   Detailed explanation goes here

if( nargin <= 4 )
    L = ((u*v) + sign*a*sqrt( u^2 + v^2 - a^2 ) ) / ( u^2 - a^2);
else
    L = ((u*v) + sign*a13*sqrt( u^2 + v^2 - a13^2 ) ) / ( u13^2 - a13^2 ) ;
end

