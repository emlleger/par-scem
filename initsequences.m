function Sequences = initsequences(D,SCEMPar,x);
% This function initializes the individual sequences

n = SCEMPar.n; q = SCEMPar.q;
% Sort the first row of D
D = sortrows(D,[1]);

for qq  = 1:SCEMPar.q,
   idx = SCEMPar.q.*([1:SCEMPar.m]-1)+qq;
   Sequences(1:SCEMPar.m,1:n+1,qq) = [x(D(idx,2),1:n) D(idx,1)];      		
end
