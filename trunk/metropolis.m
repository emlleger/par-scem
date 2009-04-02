function [C,newgen] = metropolis(pset,newpar,SCEMPar,C,Sequences,Ratio,Measurement,Extra,option)
% Metropolis rule for acceptance or rejection

Threshold = Extra.T;
% Calculate size of Sequences
[NSeq,a,b] = size(Sequences);
% Calculate number of observations
MeasData = Measurement.MeasData; N = length(MeasData(:));

for kk = 1:SCEMPar.q,

    % Sort the m points in Ckk in order of decreasing posterior density
    C(1:end,1:SCEMPar.n+1,kk)=-sortrows(-C(1:end,1:SCEMPar.n+1,kk),SCEMPar.n+1);
    % Now track last member of Markov Chain
    b = Sequences(NSeq,1:SCEMPar.n,kk); fq = Sequences(NSeq,SCEMPar.n+1,kk);
    OF = fq; uq = b;
    % Now set fm and the parameter set
    fm = pset(kk,1); new = newpar(kk,1:SCEMPar.n);
    % METROPOLIS HASTINGS selection step
    Z = rand;

    if option == 1, % Direct probability evaluation
        if (fm/fq) > Z,
            uq = new; OF = fm;
            %Randomly replace member of Ckk with new...
            C(1,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        elseif (Ratio(kk,1) > Threshold) & (OF>C(SCEMPar.m,SCEMPar.n+1,kk))
            C(SCEMPar.m,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        end;
    end;

    if option == 2 || option == 4, % Lnp probability evaluation
        qval = fm - fq;
        if qval >= log(Z),
            uq = new; OF = fm;
            C(1,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        elseif (Ratio(kk,1) > log(Threshold)) & (OF>C(SCEMPar.m,SCEMPar.n+1,kk))
            C(SCEMPar.m,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        end;
    end;

    if option == 3, % SSE probability evaluation
        if (fm/fq).^(-N.*(1+SCEMPar.Gamma)./2) > Z;
            uq = new; OF = fm;
            C(1,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        elseif (Ratio(kk,1) > Threshold) & (OF>C(SCEMPar.m,SCEMPar.n+1,kk))
            C(SCEMPar.m,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        end;
    end;
    
    if option == 5, % SSE probability evaluation weighted with sigma
        if  exp(-0.5*(-fm + fq)/Measurement.Sigma^2) > Z;            
            uq = new; OF = fm;
            C(1,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        elseif (Ratio(kk,1) > Threshold) & (OF>C(SCEMPar.m,SCEMPar.n+1,kk))
            C(SCEMPar.m,1:SCEMPar.n+1,kk) = [uq(1,1:SCEMPar.n) OF];
        end;
    end;    
    
    % Save the new location
    newgen(kk,1:SCEMPar.n+1) = [uq(1,1:SCEMPar.n) OF];
end;
