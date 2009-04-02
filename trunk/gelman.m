function [R_stat] = gelman(Sequences,SCEMPar);
% Calculates the R-statistic convergence diagnostic
% ---------------------------------------------------------------------------------
% For more information please refer to: Gelman, A. and D.R. Rubin, 1992. 
% Inference from Iterative Simulation Using Multiple Sequences, 
% Statistical Science, Volume 7, Issue 4, 457-472.
%
% Written by Jasper A. Vrugt
% Los Alamos, August 2007
% ---------------------------------------------------------------------------------

% Compute the dimensions of Sequences
[nrX,nrY,m] = size(Sequences);

% Discard the first 50% of the sequence
n = max(floor(0.5 * nrX),1); Sequences = Sequences(n+1:end,1:SCEMPar.n,1:m);

if (n < 10),
    % Set the R-statistic to a large value
    R_stat = -2 * ones(1,SCEMPar.n);
else
    % Step 1: Determine the sequence means
    meanSeq = mean(Sequences); meanSeq = reshape(meanSeq(:),SCEMPar.n,m)';
    
    % Step 1: Determine the variance between the sequence means 
    B = n * var(meanSeq);
    
    % Step 2: Compute the variance of the various sequences
    for zz = 1:SCEMPar.q,
        varSeq(zz,:) = var(Sequences(1:end,1:end,zz));
    end;
    %varSeq = var(Sequences); varSeq = reshape(varSeq(:),SCEMPar.n,m)';
    
    % Step 2: Calculate the average of the within sequence variances
    W = mean(varSeq);
    
    % Step 3: Estimate the target mean
    mu = mean(meanSeq);
    
    % Step 4: Estimate the target variance
    sigma2 = ((n - 1)/n) * W + (1/n) * B;
    
    % Step 5: Compute the R-statistic
    R_stat = sqrt((m + 1)/m * sigma2 ./ W - (n-1)/m/n);
    
end;

