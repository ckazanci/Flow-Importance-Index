% Computes "flow importance measure" given the Stoichiometric matrix

function flow_measure(S)

n = size(S,1); % number of compartments
m = size(S,2); % number of flows

feasible = combnk(1:m,n);
Nfeasible = size(feasible,1);
Npossible = 0;
NpossibleF = zeros(1,m);

for i=1:Nfeasible 
  if rank(S(:,feasible(i,:))) == n
    Npossible = Npossible + 1;
    conditions( Npossible ) = cond(S(:,feasible(i,:)));
    working_flow_set = setdiff(1:m,feasible(i,:));
    NpossibleF(working_flow_set) = NpossibleF(working_flow_set) + 1;
  end  
end
disp('======================================');
disp(['Number of compartments : ' num2str(n)]);
disp(['Number of flows        : ' num2str(m)]);
disp(['Number of flow groups of size ' num2str(m-n) '  : ' num2str(Nfeasible)]);
disp(['Number of successful flow groups : ' num2str(Npossible)])
disp('Flow measures:');
disp([ [1:m]'  NpossibleF'/Npossible] );

hist(conditions,40);
