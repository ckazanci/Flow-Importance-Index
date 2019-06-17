% Computes "flow importance measure" given the Stoichiometric matrix
% You can run as:
%   meas = importance_measure(SM);
% after running EcoNet_Results.m

function meas = importance_measure_RREF(S)

tic();         % Initialize the internal timer.
n = size(S,1); % number of compartments
m = size(S,2); % number of flows

Nfeasible = nchoosek(m,n);
feasible = zeros(Nfeasible,1);
Npossible = 0;
NpossibleF = zeros(1,m);
NpossibleFcond = zeros(1,m);
conditions = zeros(1,m);

[Npossible,NpossibleF,NpossibleFcond,conditions] = ...
    searchFlows(Npossible,NpossibleF,NpossibleFcond,conditions,...
                S,rref(S),[],[],...
                1,[1:n]*0,1:m,m);
%for i=1:Nfeasible 
%  if rank(S(:,feasible(i,:))) == n
%    Npossible = Npossible + 1;
%    conditions( Npossible ) = cond(S(:,feasible(i,:)));
%    working_flow_set = setdiff(1:m,feasible(i,:))
%    NpossibleF(working_flow_set) = NpossibleF(working_flow_set) + 1;
%    NpossibleFcond(working_flow_set) = NpossibleFcond(working_flow_set) + 1/conditions( Npossible );
%  end  
%end


totalTime = toc();

disp('======================================');
disp(['Time required for calculation : ' num2str(totalTime) ' seconds']);
disp(['Number of compartments : ' num2str(n)]);
disp(['Number of flows        : ' num2str(m)]);
disp(['Number of flow groups of size ' num2str(m-n) '  : ' num2str(Nfeasible)]);
disp(['Number of successful flow groups : ' num2str(Npossible)])
disp('Flow measures:');

meas = [ NpossibleF'/Nfeasible NpossibleFcond'/Nfeasible ]; 

figure(2)
disp([ [1:m]' meas ] );

hist(conditions,40);
