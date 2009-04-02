function C = partcomplexes(D,SCEMPar,x);
% Partition D into SCEMPar.q complexes

s = length(D);
for kk = 1:SCEMPar.q,
  idx = SCEMPar.q.*([1:SCEMPar.m]-1)+kk;
  C(1:SCEMPar.m,1:SCEMPar.n+1,kk) = [x(D(idx,2),1:SCEMPar.n) D(idx,1)];      		
end;																
