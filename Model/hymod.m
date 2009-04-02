function [OF] = hymod(Pars,Measurement,Extra); 
% Runs the HYMOD model

% Define the rainfall
PET = Extra.PET; Precip = Extra.Precip; MaxT = Extra.MaxT;
% Define the parameters
cmax = Pars(1); bexp = Pars(2); alpha = Pars(3); Rs = Pars(4); Rq = Pars(5);
% HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
x_loss = 0.0;
% Initialize slow tank state
x_slow = 2.3503/(Rs*22.5);
% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0; t=1; outflow = [];
% Time loop
while t < (MaxT + 1),

    Pval = Precip(t,1); PETval = PET(t,1);
    
    % Compute excess precipitation and evaporation
    %[UT1,UT2,x_loss] = excess(x_loss,cmax,bexp,Pval,PETval);
    
    % -------------------------------------------------
    ct_prev = cmax*(1 - ((1-((bexp+1)*(x_loss)/cmax)).^(1/(bexp+1))));
    UT1 = max((Pval-cmax + ct_prev),0.0);
    Pval = Pval - UT1;
    dummy = min(((ct_prev + Pval)/cmax),1);
    xn = (cmax/(bexp+1))*(1 - ((1-dummy).^(bexp+1)));
    UT2 = max(Pval-(xn-x_loss),0);
    evap = min(xn,PETval);
    xn = xn-evap; x_loss = xn;
    % -------------------------------------------------
    
    % Partition UT1 and UT2 into quick and slow flow component
    UQ = alpha*UT2 + UT1; US = (1-alpha)*UT2;
    
    % Route slow flow component with single linear reservoir
    inflow = US; %[x_slow,outflow] = linres(x_slow,inflow,outflow,Rs); QS = outflow;
    
    % Define x_slow and outflow
    x_slow = (1-Rs)*x_slow + (1-Rs)*inflow; outflow = (Rs/(1-Rs))*x_slow; QS = outflow;

    % Route quick flow component with linear reservoirs
    inflow = UQ;

    % Then do three times for quickflow reservoirs
    x_quick(1) = (1-Rq)*x_quick(1) + (1-Rq)*inflow; outflow = (Rq/(1-Rq))*x_quick(1); inflow = outflow;
    x_quick(2) = (1-Rq)*x_quick(2) + (1-Rq)*inflow; outflow = (Rq/(1-Rq))*x_quick(2); inflow = outflow;
    x_quick(3) = (1-Rq)*x_quick(3) + (1-Rq)*inflow; outflow = (Rq/(1-Rq))*x_quick(3); inflow = outflow;

    %while k < 4,
    %   [x_quick(k),outflow] = linres(x_quick(k),inflow,outflow,Rq); inflow = outflow;
    %   k = k+1;
    %end;
    % Compute total flow for timestep
    output(t,1) = QS + outflow;
    t = t + 1;
end;

SimRR = 22.5 * output(65:MaxT,1);

% Compute the objective function
OF = sum((SimRR - Measurement.MeasData).^2);

