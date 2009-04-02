function [D,x] = reshuffle(C,SCEMPar);  
% Function performs reshuffling

counter = 1;
for qq=1:SCEMPar.q,
   for ii=1:SCEMPar.m
      x(counter,1:SCEMPar.n) = C(ii,1:SCEMPar.n,qq);
      D(counter,1:2) = [C(ii,SCEMPar.n+1,qq) counter];
      counter = counter + 1;
   end;
end;
D = -sortrows(-D,1);
