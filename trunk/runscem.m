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
% SCEM-UA orignally developed by Jasper A. Vrugt                                                %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution           %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
% And used in:                                                                                  %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, S. Sorooshian, Effective and efficient    %
%       algorithm for multiobjective optimization of hydrologic models, Water Resour. Res., 39  %
%       (8), 1214, doi:10.1029/2002WR001746, 2003.                                              %
%                                                                                               % 
%   Vrugt J.A., S.C. Dekker, W. Bouten, Identification of rainfall interception model           %
%       parameters from measurements of throughfall and forest canopy storage, Water Resour.    %
%       Res., 39 (9), 1251, doi:10.1029/2003WR002013, 2003.                                     % 
%                                                                                               %
%   Vrugt J.A., G. Schoups, J.W. Hopmans, C. Young, W.W. Wallender, T. Harter, W. Bouten,       %
%       Inverse modeling of large-scale spatially distributed vadose zone properties using      %
%       global optimization, Water Resour. Res., 40, W06503, doi:10.1029/2003WR002706, 2004.    %
%                                                                                               % 
%   Vrugt J.A. W. Bouten, H.V. Gupta, and J.W. Hopmans, Toward improved identifiability of soil % 
%       hydraulic parameters: on the selection of a suitable parametric model, Vadose Zone      %
%       Journal, 2, 98 - 113, 2003.                                                             %
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by differential evolution with         %
%       self-adaptive randomized subspace sampling, International Journal of Nonlinear          %
%       Sciences and Numerical Simulation, In Press. 2008.                                      %
%                                                                                               %
%   ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%       input uncertainty in hydrologic modeling: Doing hydrology backwards using Markov        %
%       chain Monte Carlo, Water Resour. Res., In Press. 2008.                                  %
%                                                                                               %
% MATLAB code written by Jasper A. Vrugt, Center for NonLinear Studies (CNLS)                   %
%                                                                                               %
% Version 0.5: June 2002                                                                        %
%                                                                                               %
% ----------------------------------------------------------------------------------------------%

% How many nodes are being used
available_nodes = 6;

% Define global variables
global tag = 42; global nodes = -1; global rank = -1; global world = -1;



% Initialize global variables
global SCEMPar; global ModelName; global Extra; global Measurement; 

% Define problem specific variables
SCEMPar.n = 5;                          % Dimension of the problem
SCEMPar.q = 5;                          % Number of complexes
SCEMPar.s = 50;                         % Population size
SCEMPar.ndraw = 100; %5000;                   % Maximum number of function evaluations
SCEMPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
  
% Calculate the number of points in each complex, m and the number of offspring generated per iteration, L
SCEMPar.m = floor(SCEMPar.s./SCEMPar.q); 
SCEMPar.L = max([1,floor(SCEMPar.m./5)]);

% Define the parameter ranges
ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; ParRange.maxn = [500 2.00 0.99 0.10 0.99];

% Load the Leaf River data
load bound.txt;

% Then read the boundary conditions -- two years only
Extra.MaxT = 795; 

% Define the PET, Measured Streamflow and Precipitation.
Extra.PET = bound(1:Extra.MaxT,5); Extra.Precip = sum(bound(1:Extra.MaxT,6:9),2);

Extra.Jump = (2.4/sqrt(SCEMPar.n))^2; Extra.T = 1e7;

% Define the measured streamflow data
Measurement.MeasData = bound(65:Extra.MaxT,4); Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);

% Define the modelname 
ModelName = 'hymod';

option = 3; % SSE evaluation

% Define total number of nodes 
% nodes = available_nodes;

% Store the current directory
Extra.pwd = pwd;

% Initialize the MPI network
[rank,nodes,world] = initmpi(caste,available_nodes);

% Run the master and slaves
if (caste == master)
   [ParSet,GR,Sequences] = scem_ua(SCEMPar,ParRange,Measurement,ModelName,Extra,option);
else
   scem_slave      
end;
printf("I reached finalize(%d)", rank);
MPI_Finalize;
