function [pset] = computedensity(x,SCEMPar,Measurement,ModelName,Extra,option);
% This function calculates the posterior density

MeasData = Measurement.MeasData;
Sigma = Measurement.Sigma;
Cb = SCEMPar.Cb;
Wb = SCEMPar.Wb;

[NrParSets,r] = size(x);

% For each sampled point compute function value and calculate posterior density ...
for ii=1:NrParSets,
    % Call model to generate simulated data
    evalstr = ['SimData = ',ModelName,'(x(ii,:),Extra);'];
    eval(evalstr);

    if option == 1,
        % Model directly computes posterior density
        p = SimData;
        pset(ii,1:2) = [p ii];
    end;

    if option == 2,
        N = length(MeasData(:));
        Err = (MeasData(:)-SimData(:));
        % Equation 5 from the manual
        % In situations with large values for N (many points), these values become very small...
        % taking the logarithm is a good way of handling these small values...
        lnp = N.*log(Wb./Sigma) - Cb.*(sum((abs(Err./Sigma)).^(2/(1+SCEMPar.Gamma))));
        pset(ii,1:2) = [lnp ii];
    end;

    if option == 3,
        % Determine the error residual
        Err = (MeasData(:)-SimData(:));
        % Part of equation 8 from the manual
        SSE = sum(abs(Err).^(2/(1+SCEMPar.Gamma))); 
        pset(ii,1:2) = [-SSE ii]; % Minus is added so minimum M corresponds to highest likelihood
    end;
    
    if option == 4, % Model directly computes log posterior density
        pset(ii,1:2) = [SimData ii]; 
    end;

    if option == 5, % Similar as 3, but now weights with the Measurement Sigma
        % Determine the error residual
        Err = (Measurement.MeasData(:)-SimData(:));
        % Derive the sum of squared error
        SSE = sum(abs(Err).^(2/(1+SCEMPar.Gamma)));
        % And retain in memory
        pset(ii,1:2) = [-SSE ii]; 
    end;

end;
