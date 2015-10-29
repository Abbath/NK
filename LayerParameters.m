function [ F, alfa ] = LayerParameters( MaxLay, Ndx )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F = zeros(1, MaxLay+2);
alfa = zeros(1, MaxLay+2);
for lay = 1:MaxLay+2
    F(lay) = Fresnel(Ndx(lay));
end
alfa(2) = F(2);
for lay = 3:MaxLay+2
    negF = -F(lay-1);
    alfa(lay) = Compos(F(lay), negF);
end
end

