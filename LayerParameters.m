function [ F, alfa ] = LayerParameters( MaxLay, Ndx )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
lay = 0;
negF = 0;
for lay = 1:MaxLay
    F(lay) = Fresnel(Ndx(lay));
end
alfa(1) = F(1);
for lay = 2:MaxLay
    negF = -F(lay-1);
    alfa(lay) = Compos(F(lay), negF);
end
end

