function [Gauss_x] = GaussScaling(y,a,b)
% This function scales the interval [-1,1] from GaussIntegration.m to [a,b]

tempx=(a*(1-y)+b*(1+y))/2;
Gauss_x=unique(tempx);

return
