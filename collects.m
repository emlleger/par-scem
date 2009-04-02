function [Sequences,ParSet] = collects(SCEMPar,Sequences,newgen,ParSet);   
% Collect results in sequences

% Calculate size of Sequences
[NSeq,a,b] = size(Sequences);

for kk = 1:SCEMPar.q,
    Sequences(NSeq+1,1:SCEMPar.n+1,kk) = newgen(kk,1:SCEMPar.n+1);
end;

ParSet = [ParSet;newgen]; 
   