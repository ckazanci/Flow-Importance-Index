% model2dotv( values, NM, SM)
% This Matlab script is named after the bash script named model2dot
% and adds "v"alues to edges. The output is "model.eps".

% The other two variables needed are:
% NM: Compartment names
% SM: Stoichiometric matrix
% These two variable names can be kept as they are if one calls this function
% after running "EcoNet_Results.m".

% To enumerate flows, after running EcoNet results, run:
% model2dotv((1:size(SM,2))',NM,SM)
% To see the numerical values of importance measures, run:
% model2dotv(meas(:,1),NM,SM)


function model2dotv( values, NM, SM )

Ndigits = 2;
  
% EcoNet_Results;
n = size(SM,1);  % #compartments
m = size(SM,2);  % #flows

% Edge width adjustment
if size(values,1) ~= m
  disp(['Size of magnitudes is = ' num2str(size(values,1))])
  disp(['But number of flows is = ' num2str(m)])
end

% linear
v = round(10^Ndigits*values)/10^Ndigits;

pp{1} = 'digraph model {';
pp{2} = 'graph [bgcolor=transparent]';
pp{3} = 'node [shape=point color=red];';

fType = sum(SM);
n_in_out = size( find(fType ~= 0) ,2);
for i=1:n_in_out
  pp{3+i} = [ 'o' num2str(i) ];
end

pp{4+n_in_out} = 'node [shape=ellipse style=filled fillcolor=cadetblue1 fontname=Helvetica];'; % lightblue2
pp{5+n_in_out} = 'edge [fontname=Helvetica fontcolor=blue3];';

ctr_in_out = 1;
for i=1:m
  if fType(i) == 0
    L = NM{find(SM(:,i)==-1)};
    R = NM{find(SM(:,i)==1)};
    pp{5+n_in_out+i} = [ L ' -> ' R ' [label=<<table cellpadding="3" border="0" cellborder="0"><tr><td>' num2str(v(i)) '</td></tr></table>>];'];
  elseif fType(i) == 1
    L = [ 'o' num2str(ctr_in_out) ];
    ctr_in_out = ctr_in_out + 1;
    R = NM{find(SM(:,i)==1)};
    pp{5+n_in_out+i} = [ L ' -> ' R ' [label=<<table cellpadding="3" border="0" cellborder="0"><tr><td>' num2str(v(i)) '</td></tr></table>>];'];
  elseif fType(i) == -1
    L = NM{find(SM(:,i)==-1)};
    R = [ 'o' num2str(ctr_in_out) ];
    ctr_in_out = ctr_in_out + 1;    
    pp{5+n_in_out+i} = [ L ' -> ' R ' [label=<<table cellpadding="3" border="0" cellborder="0"><tr><td>' num2str(v(i)) '</td></tr></table>>];'];
  end 
end

pp{6+n_in_out+m} = '}';

% store to file

fid = fopen('model.dot', 'wt');
for i=1:length(pp)
  fprintf( fid, '%s\n', pp{i});
end
fclose(fid);

!dot -Tpng model.dot -o model.png
!dot -Tps model.dot -o model.ps

