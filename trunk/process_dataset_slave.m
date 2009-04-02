function process_dataset_slave()
  global world; global tag; global nodes; global rank; global SCEMPar; global Extra; global ModelName; global Measurement; global scratch_dir;

	 % And go to the appropriate directory
    evalstr = strcat('cd ',scratch_dir); eval(evalstr); 

    [minfo mstat] = MPI_Probe(MPI_ANY_SOURCE,tag,world);
    [minfo dbls] = MPI_Get_elements(mstat,[]);
    N = zeros(1,dbls);
    [minfo mstat] = MPI_Recv(N,MPI_ANY_SOURCE,tag,world);
    % Get the parameter sets
    x = reshape(N, dbls/SCEMPar.n, SCEMPar.n);
    % x is a vector of parameter sets
    [nrsets  nrpars] = size(x);   

    % Then evaluate each combination of x
    % out is a vector of model evaluations per each parset
    for i=1:nrsets
        evalstr = ['OF = ',ModelName,'(x(i,1:nrpars),Measurement,Extra);']; eval(evalstr);
        out(i,:) = OF;
    endfor

    % Return the result   
    [dummy nrelts] = size(OF);
    package = reshape(out, 1, nrsets*nrelts);
    MPI_Send(package,0,tag,world);
    clear out;

endfunction