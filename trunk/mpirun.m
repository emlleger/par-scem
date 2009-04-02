function [ModPred] = mpirun(Pop,Nelem,nodes);
% Function that takes existing population and distributes over pre-specified number of computers 
% 
% Input -- Pop:         Population with parameter values
% 
% Output - ModPred:     The model output corresponding to each row of Pop

global tag; global world; 

% First determine size of the population
[Npop,n] = size(Pop);

 printf("\nPop size %d,%d%d\n",Npop,n, nodes);
% First make individual packages of Pop and send to slaves.
for slave = 1:nodes-1
    [lhs,rhs] = divvy(Npop,nodes-1,slave);
    package = reshape(Pop(lhs:rhs,:),1,(rhs-lhs+1)*n);
    MPI_Send(package,slave,tag,world);
end;

% Receive from slaves.
for slave = 1:nodes-1
    [lhs rhs] = divvy(Npop,nodes-1,slave);
    [minfo mstat] = MPI_Probe(slave,tag,world);
    [minfo dbls] = MPI_Get_elements(mstat,[]);
    N = zeros(1,dbls);
    [minfo mstat] = MPI_Recv(N,slave,tag,world);
    ModPred(lhs:rhs,:) = reshape(N,rhs-lhs+1,Nelem);
end
% ###############################################################
