function [ParSet,GR,Sequences] = scem_ua(SCEMPar,ParRange,Measurement,ModelName,Extra,option);
% ------- Shuffled Complex Evolution Metropolis (SCEM-UA) global optimization algorithm --------%
%                                                                                               %
% The SCEM-UA algorithm is a multi-chain general purpose stochastic optimization algorithm that %
% uses adaptive MCMC sampling to provide an efficient search of the parameter space and         %
% estimate the most likely parameter combination and its underlying posterior distribution      %
% within a single optimization (sampling) run. The SCEM-UA method first starts with             %
% stategically placing s different points in the feasible parameter space using Latin Hypercube %
% sampling. The posterior density of each of the initial s points is calculated, and q Markov   %
% chains (or sequences) are launched starting from the q points of s that exhibit the highest   %
% density. These Markov chains are used to independently explore the search space. After        %
% initialization, new (or candidate) points in each chain (or sequence) are generated using a   %
% multivariate normal proposal distribution with mean centered at the current value of the      %
% chain, and covariance estimated from the set of points that define the current of k complexes %
% The Metropolis acceptance rule is then evaluated to check whether the chains should jump to   %
% their candidate points or whether they should stay at their current value. Similarly, the     %
% Metropolis rule is evaluated to check whether the newly generated points should replace       %
% randomly chosen members of D. The MCMC evolution of the various chains is repeated until      %
% the R-statistic of Gelman and Rubin (1992) indicates convergence to a stationary posterior    %
% distribution. Various examples have shown that SCEM-UA is generally efficient for posterior   %
% exploration. However, the method does not maintain detailed balance at every single step in   %
% the chain. We therefore recommend using the recently developed DiffeRential Evolution         %
% Adaptive Metropolis algorithm, which maintains balance and is fully ergodic. This method      %
% provides an EXACT estimate of the posterior probability density function and is especially    %
% designed to efficiently handle multimodality, nonlinearity and high-dimensionality.           %
%                                                                                               %
% SCEM-UA developed by Jasper A. Vrugt, University of Amsterdam                                 %
%                                                                                               %
% Version 0.5: June 2002                                                                        %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Initialize two important variables
Iter = SCEMPar.s; 								



% Sample s points in the parameter space using latin hypercube
x = lhsu(ParRange.minn,ParRange.maxn,SCEMPar.s);

% Calculate the parameters in the exponential power density function of Box and Tiao
[SCEMPar.Wb,SCEMPar.Cb] = calccbwb(SCEMPar.Gamma);

% Step 2: Parallel implementation: compute posterior density???
% ------------------------------- PARALLEL ---------------------------------------------------   
[p,dummy] = DoParallel(x); pset = p(:,1); pset(:,2) = [1:SCEMPar.s]';
% --------------------------------------------------------------------------------------------  

% Save the output in the array ParSet
ParSet = [x,pset(1:end,1)];

% Sort the points in order of decreasing posterior probability
D = -sortrows(-pset,1);

% Initialize starting points of sequences 
Sequences = initsequences(D,SCEMPar,x);

% Compute convergence
[Rstat] = gelman(Sequences,SCEMPar); GR = [Iter Rstat];

% Iterate until a predefined number of function evaluations has been readched ...
while (Iter < SCEMPar.ndraw),

    % Partition D into c complexes A
    C = partcomplexes(D,SCEMPar,x);

    % Evolve each sequence with the Sequence Evolution Metropolis (SEM) algorithm
    [C,Sequences,ParSet] = sem(C,Sequences,ParSet,ParRange,SCEMPar,Measurement,ModelName,Extra,option,Iter);
    
    % Reshuffle the points
    [D,x] = reshuffle(C,SCEMPar);
    
    % Compute convergence
    [Rstat] = gelman(Sequences,SCEMPar);
    
    % Update the Iteration
    Iter = Iter + SCEMPar.L*SCEMPar.q;
    
    % Now save the convergence statistic
    GR = [GR; Iter Rstat];

    % Save all the information so far
    eval(strcat('save -binary ','tempState',num2str(Iter)))

end;

% Write all the output
save -binary outSCEM


