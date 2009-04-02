function [C,Sequences,ParSet] = sem(C,Sequences,ParSet,ParRange,SCEMPar,Measurement,ModelName,Extra,option,Iter);
% The sequence evolution metropolis algorithm for Markov chain evolution
%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

for bb = 1:SCEMPar.L,
    % First generate new candidate point for each complex and sequence
    [newpar,Ratio] = offmetro(Sequences,C,SCEMPar,ParRange,ModelName,Measurement,Extra,option);
    
    % ------------------------------- PARALLEL ---------------------------------------------------
    CopyOrNot = 0; [pset,dummy] = DoParallel(newpar); pset(:,2) = [1:SCEMPar.q]';
    % --------------------------------------------------------------------------------------------

    % Now compute the likelihood of the new points
    %[pset] = computedensity(newpar,SCEMPar,Measurement,ModelName,Extra,option);
    % Now apply the acceptance/rejectance rule
    [C,newgen] = metropolis(pset,newpar,SCEMPar,C,Sequences,Ratio,Measurement,Extra,option);
    % Now collect results in Sequences and total file
    [Sequences,ParSet] = collects(SCEMPar,Sequences,newgen,ParSet);
    % reset negen
    newgen = [];
end;
