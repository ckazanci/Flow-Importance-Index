% Computes "flow importance measure" given the 
% Stoichiometric matrix S, flow rates that are already 
% determined (fk: f known), and flow rates that are very 
% difficult to determine (fu: f unknown).

% After running EcoNet_Results.m, type:
% meas = importance_measure(SM, [1 2],[3 4]);

% This code supersedes importance_measure.m 
% in that output of both are identical:
% meas = importance_measure(SM);
% meas = importance_measure(SM, [],[]);

function meas = importance_measure_partial(S,fk,fu)

n = size(S,1); % number of compartments
k = size(S,2); % number of flows
Nfk = size(fk,2);
Nfu = size(fu,2);

% flow indices that are not "known or unknown"
k_ind = setdiff(1:k,[fk fu]); % actual flow indices
Ksubsets = combnk(k_ind,n-Nfu);
nKsubsets = size(Ksubsets,1); % number of k subsets
% number of subsets of size k inlcuding flow i
nKsubsets_i = nKsubsets*(k-n-Nfk)/(k-Nfu-Nfk);

Nfeas = 0; % number of feasible k subsets 
% number of feasible k subsets inlcuding flow i
Nfeas_i = zeros(1,k); 
% Conditions sum of feasible k subsets including i
sum_cond = zeros(1,k); 
% Sum of inverse conditions of feasible k sub...
sum_invcond = zeros(1,k); 

for i = 1:nKsubsets 
  if rank(S(:,[Ksubsets(i,:) fk])) == n
    Nfeas = Nfeas + 1;
    conditions( Nfeas ) = cond(S(:,[Ksubsets(i,:) fk]));
    for j = setdiff(k_ind,Ksubsets(i,:))
      Nfeas_i(j) = Nfeas_i(j) + 1;
      sum_cond(j) = sum_cond(j) + conditions( Nfeas );
      sum_invcond(j) = sum_invcond(j)+1/conditions( Nfeas );
    end
  end  
end
format compact
disp('======================================');
disp(['Number of compartments : ' num2str(n)]);
disp(['Number of flows        : ' num2str(k)]);
disp(['Number of flow groups  : ' num2str(nKsubsets)]);
disp(['including a specific flow : ' num2str(nKsubsets_i)]);
disp(['Number of feasible groups : ' num2str(Nfeas)]);
disp('Flow measures (-1:known -2:unknown):');
meas1 = Nfeas_i'/nKsubsets_i;
meas2 = sum_invcond'.*sum_cond'./Nfeas_i'/nKsubsets_i;
meas = [ meas1 meas2 ]; 
meas(fk,:) = -1; % set known flows to -1
meas(fu,:) = -2; % set unknown flows to -2
disp([ k_ind' meas(k_ind,:) ] );
disp(['   max   : ' num2str(max(meas(k_ind,:)))]);
disp(['   min   : ' num2str(min(meas(k_ind,:)))]);
disp(['   mean  : ' num2str(mean(meas(k_ind,:)))]);
disp(['   std   : ' num2str(std(meas(k_ind,:)))]);
figure;
set(gcf,'Position',[500 500 1000 400]);
subplot(1,2,1) % histogram of condition numbers 
hist(conditions,25);
subplot(1,2,2) % importance measure, with or w/o cond num
minM = min(min(meas(k_ind,:)));
maxM = max(max(meas(k_ind,:)));
plot(meas1(k_ind,:),meas2(k_ind,:),'o'); hold on
plot([0 1],[0 1],':')
axis([minM*.98 maxM*1.02 minM*.98 maxM*1.02])
xlabel('without condition number')
ylabel('with condition number')
% close all
