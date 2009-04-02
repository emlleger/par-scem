function [Wb,Cb] = calccbwb(Beta);
% This function calculates the parameters for the exponential power density
% Equation [20] paper by Thiemann et al. WRR 2001, Vol 37, No 10, 2521-2535
A1 = gamma(3*(1+Beta)/2); 
A2 = gamma((1+Beta)/2); 
Cb = (A1/A2)^(1/(1+Beta));
Wb = sqrt(A1)/((1+Beta)*(A2^(1.5)));


