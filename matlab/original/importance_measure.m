% Computes "flow importance measure" given the 
% Stoichiometric matrix S
% After running EcoNet_Results.m, type:
% meas = importance_measure(SM);

function meas = importance_measure(S)

n = size(S,1); % number of compartments
k = size(S,2); % number of flows

Ksubsets = combnk(1:k,n);
nKsubsets = size(Ksubsets,1); % number of k subsets
% number of subsets of size k inlcuding flow i
nKsubsets_i = nKsubsets*(k-n)/k;

Nfeas = 0; % number of feasible k subsets 
% number of feasible k subsets inlcuding flow i
Nfeas_i = zeros(1,k); 
% Conditions sum of feasible k subsets including i
sum_cond = zeros(1,k); 
% Sum of inverse conditions of feasible k sub...
sum_invcond = zeros(1,k); 

for i = 1:nKsubsets 
  if rank(S(:,Ksubsets(i,:))) == n
    Nfeas = Nfeas + 1;
    conditions( Nfeas ) = cond(S(:,Ksubsets(i,:)));
    for j = setdiff(1:k,Ksubsets(i,:))
      Nfeas_i(j) = Nfeas_i(j) + 1;
      sum_cond(j) = sum_cond(j) + conditions( Nfeas );
      sum_invcond(j) = sum_invcond(j) + 1/conditions( Nfeas );
    end
  end  
end
% Nfeas_i
% sum_cond
% sum_invcond
format compact
disp('======================================');
disp(['Number of compartments : ' num2str(n)]);
disp(['Number of flows        : ' num2str(k)]);
disp(['Number of flow groups  : ' num2str(nKsubsets)]);
disp(['including a specific flow : ' num2str(nKsubsets_i)]);
disp(['Number of feasible groups : ' num2str(Nfeas)]);
disp('Flow measures:');
meas1 = Nfeas_i'/nKsubsets_i;
meas2 = sum_invcond'.*sum_cond'./Nfeas_i'/nKsubsets_i;
meas = [ meas1 meas2 ]; 
disp([ [1:k]' meas ] );
disp(['   max   : ' num2str(max(meas))]);
disp(['   min   : ' num2str(min(meas))]);
disp(['   mean  : ' num2str(mean(meas))]);
disp(['   std   : ' num2str(std(meas))]);
figure;
set(gcf,'Position',[500 500 1000 400]);
subplot(1,2,1) % histogram of condition numbers 
hist(conditions,25);
subplot(1,2,2) % importance measure, with or w/o cond num
minM = min(min(meas));
maxM = max(max(meas));
plot(meas1,meas2,'o'); hold on
plot([0 1],[0 1],':')
axis([minM*.98 maxM*1.02 minM*.98 maxM*1.02])
xlabel('without condition number')
ylabel('with condition number')
% close all

