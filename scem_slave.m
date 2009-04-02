% Wait for data, process it and return it.

% Termination?
function scem_slave()
	global world; global tag; global nodes; global rank; global SCEMPar; global Extra; global ModelName; global Measurement;
	global scratch_dir = strcat('/Users/admin/workspace/ParSCEM/scem_ua',num2str(rank));
	global model_dir = '$MODELDIR ';
	
	% Copy files from the model directory to the scratch space of the node
  	evalstr = strcat('mkdir ',scratch_dir); unix(evalstr);
  	
    evalstr = strcat('cp -pfr ',model_dir, scratch_dir); unix(evalstr);
    
	process_dataset_slave; 	
	evalstr = strcat('cd ',Extra.pwd); eval(evalstr); 
 
 Iter = SCEMPar.s;
 while (Iter < SCEMPar.ndraw)
 for bb = 1:SCEMPar.L
	process_dataset_slave; 	
	evalstr = strcat('cd ',Extra.pwd); eval(evalstr); 
 end;
 Iter = Iter + SCEMPar.L*SCEMPar.q;
 endwhile
endfunction

% eof.
