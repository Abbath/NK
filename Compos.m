function [ res ] = Compos( A, B )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
X = A + B;
Y = A * B;
Y = complex(real(Y) + 1, imag(Y));
res = X / Y;
end

