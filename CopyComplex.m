function [ beta ] = CopyComplex( MaxLay, n )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
beta = zeros(1, MaxLay);
for i = 1:MaxLay
   beta(i) = n(i); 
end

end

