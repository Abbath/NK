function [ res ] = Fresnel( N )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
F = 0;
Num = 0;
Den = 0;
Num = N - 1;
Den = N + 1;
res = Num / Den;
end

