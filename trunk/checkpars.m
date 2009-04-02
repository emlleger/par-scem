function [T] = checkpars(NewPars,ParRange);
% Function checks whether new vector with parameter values is within boud or not
% T = 0 => parameter set is feasible; T = 1 => parameter set is infeasible

minn = ParRange.minn; maxn = ParRange.maxn;

T = 0; 
[m,n] = size(NewPars);
for ii = 1:n,
  if NewPars(1,ii) < minn(1,ii),
    T = -1;
  end;
  if NewPars(1,ii) > maxn(1,ii),
    T = -1;
  end;
end;
