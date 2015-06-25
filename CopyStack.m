function [ d_to ] = CopyStack( MaxLay, d_from )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
d_to = zeros(1, MaxLay);
for i = 1:MaxLay
   d_to(i) = d_from(i); 
end

end

