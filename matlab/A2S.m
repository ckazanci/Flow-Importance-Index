function S = A2S(A)
% S = A2S(A) computes the Stoichiometric matrix of an ecosystem 
% model using its adjacency matrix, where the last row and column
% represents the environment. For an n-compartment model, A is an  
% (n+1)x(n+1) square matrix where: 
% - A(i,j) represents the flow from compartment j to i, for 0<i,j<n+1 
% - A(i,n+1) represents the environmental input into compartment i. 
% - A(n+1,j) represents the environmental output from compartment j. 

if ( nargin ~= 1 )
  help A2S
  return
end

if ( size(A,1) ~= size(A,2) ) 
  fprintf('\n Adjacency matrix has to be a square matrix!\n\n');
  return
end

if ( length(find( A~=0 & A~=1 )) ~= 0 ) 
  fprintf('\n Adjacency matrix can only take values of 0 and 1!\n\n');
  return
end
  
if ( length(find( diag(A)==1 )) ~= 0  ) 
  fprintf('\n Self loops are not allowed!\n\n');
  return
end

n = size(A,1)-1; % number of compartments

k = 1;
for i=1:n
  for j=1:n
    if A(i,j) == 1
      S(j,k) = -1;
      S(i,k) = 1;
      k = k + 1;
    end
  end
  if A(i,n+1) == 1
    S(i,k) = 1;
    k = k + 1;
  end
  if A(n+1,i) == 1
    S(i,k) = -1;
    k = k + 1;
  end
end
