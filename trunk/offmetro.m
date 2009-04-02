function [new,Ratio] = offmetro(Sequences,C,SCEMPar,ParRange,ModelName,Measurement,Extra,option);
% Generates offspring using Markov Chain Monte Carlo sampling

MeasData = Measurement.MeasData; N = length(MeasData(:));
% Initialize JumpRate and T - value
JumpRate = Extra.Jump; Threshold = Extra.T;
% Calculate size of Sequences
[NSeq,a,b] = size(Sequences);

% Loop over the individual sequences
for kk = 1:SCEMPar.q,

    % Compute the covariance of complex k;
    COVCOMPLEX = cov(C(1:end,1:SCEMPar.n,kk));
    % Compute the mean likelihood in the complex
    PCom = mean(C(1:SCEMPar.m,SCEMPar.n+1,kk));
    % Compute the mean likelihood over the last m points of the sequence
    PSeq = mean(Sequences(max(NSeq-SCEMPar.m,1):NSeq,SCEMPar.n+1,kk));
    % Now track last member of Markov Chain
    b = Sequences(NSeq,1:SCEMPar.n,kk); 
    % Now start looping to generate new parameter combination that falls within desired hypercube
    accept = -1;

    while accept == -1,
        ru = randn(SCEMPar.n,1);
        
        if option == 1,
            if (PCom/PSeq) > Threshold,
                new(kk,1:SCEMPar.n) = real(mean(C(1:SCEMPar.m,1:SCEMPar.n))+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            else
                new(kk,1:SCEMPar.n) = real(b+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            end;
            Ratio(kk,1) = C(1,SCEMPar.n+1)/C(SCEMPar.m,SCEMPar.n+1);
        end;
        
        if option == 2 || option == 4,
            if (PCom-PSeq)> log(Threshold),
                new(kk,1:SCEMPar.n) = real(mean(C(1:SCEMPar.m,1:SCEMPar.n))+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            else
                new(kk,1:SCEMPar.n) = real(b+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            end;
            Ratio(kk,1) = C(1,SCEMPar.n+1)-C(SCEMPar.m,SCEMPar.n+1);
        end;
        
        if option == 3,
            if (PCom/PSeq).^(-N.*(1+SCEMPar.Gamma)./2) > Threshold
                new(kk,1:SCEMPar.n) = real(mean(C(1:SCEMPar.m,1:SCEMPar.n))+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            else
                new(kk,1:SCEMPar.n) = real(b+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            end;
            Ratio(kk,1) = (C(1,SCEMPar.n+1)/C(SCEMPar.m,SCEMPar.n+1)).^(-N.*(1+SCEMPar.Gamma)./2);
        end;
        
        if option == 5,
            if exp(-0.5*(-PCom + PSeq)/Measurement.Sigma^2) > Threshold            
                new(kk,1:SCEMPar.n) = real(mean(C(1:SCEMPar.m,1:SCEMPar.n))+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            else
                new(kk,1:SCEMPar.n) = real(b+(sqrtm(JumpRate*COVCOMPLEX)*ru)');
            end;
            Ratio(kk,1) = (C(1,SCEMPar.n+1)/C(SCEMPar.m,SCEMPar.n+1)).^(-N.*(1+SCEMPar.Gamma)./2);
        end;
        
        accept = checkpars(new(kk,1:SCEMPar.n),ParRange);        
        
    end;

end;