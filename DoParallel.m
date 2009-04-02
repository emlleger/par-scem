function [p,log_p] = DoParallel(x);
% Distributes the parameter combinations and runs the model on each individual node

global nodes; 

% Then evaluate each individual of x
 [OF] = mpirun(x,1,nodes);
% [l,c]= size(OF);
% printf("\nOF size %d,%d\n",l,c);

for ii = 1:size(x,1),
   % Retain in memory
   p(ii,1) = [-OF(ii,1)]; log_p(ii,1) = -0.5 * OF(ii,1);
end;
