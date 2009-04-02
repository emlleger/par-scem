function [rank,nodes,world] = initmpi(caste,available_nodes);
% Initializes the MPI environment

global master; global slave; global Extra;

% Slave messages
slave_terminate = 42;
slave_computedensity = 43;
slave_SEM = 44;

MPI_Init;
[info flag] = MPI_Initialized;
[info rank] = MPI_Comm_rank(MPI_COMM_SELF);

slave_file = strcat(Extra.pwd , '/slave.m');


% Spawn and merge communicators.
switch(caste)
    case master
        [info children errs] = ...
            MPI_Comm_spawn('octave', ...
            {"-q", slave_file}, ...
            available_nodes-1, MPI_INFO_NULL, 0, MPI_COMM_SELF);
        [info world] = MPI_Intercomm_merge(children, 0);
    case slave
        [info parent] = MPI_Comm_get_parent;
        [info world] = MPI_Intercomm_merge(parent, 0);
    otherwise
        error("Caste unassigned %d\n", caste);
endswitch;

% Get nodes and rank.
[info rank] = MPI_Comm_rank(world); 
[info nodes] = MPI_Comm_size(world);
[info pname] = MPI_Get_processor_name;
printf("Caste: %d Rank: %d/%d(%s)\n", caste, rank, nodes, pname);
