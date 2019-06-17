% model2dots( magnitudes, levels, NM, SM, linlog )
% This Matlab script is named after the bash script named model2dot
% and is basically a "s"caled version of that. It creates a network diagram 
% with edges in varying thickness, based on a scale vector "magnitudes" and 
% number of thickness grades "levels", with using a linear scale ("linlog=1")
% or logarithmic scale ("linlog=2"). The output is "model.ps".

% The other two variables needed are:
% NM: Compartment names
% SM: Stoichiometric matrix
% These two variable names can be kept as they are if one calls this function
% after running "EcoNet_Results.m".

% model2dots(meas(:,1),5,NM,SM,1)

function model2dots( magnitudes, levels, NM, SM, linlog )

Ndigits = 2;
% EcoNet_Results;
n = size(SM,1);  % #compartments
m = size(SM,2);  % #flows

% Edge width adjustment
if size(magnitudes,1) ~= m
  disp(['Size of magnitudes is = ' num2str(size(magnitudes,1))])
  disp(['But number of flows is = ' num2str(m)])
end

% indentify regular flows from known and unknown ones 
% by checking negative measures
ind = find( magnitudes > 0 ); 

% linear
v = magnitudes;
vmax = max(v(ind));
vmin = min(v(ind));
vspace = (vmax-vmin)/levels;
% trick to prevent the max value be not 1 more than what it should be
vtemp = v;
vtemp(find(v==vmax))= vmax - vspace/10;
widths = 1 + floor((vtemp-vmin)/vspace);
 
% exponential
u = log(magnitudes);
umax = max(u(ind));
umin = min(u(ind));
uspace = (umax-umin)/levels;
% trick to prevent the max value be not 1 more than what it should be
utemp = u;
utemp(find(u==umax)) = umax - uspace/10;
lwidths = 1 + floor((utemp-umin)/uspace);

% Display limits:
limlin = vmin:vspace:vmax;
limlog = exp(umin:uspace:umax);
disp([limlin' limlog'])


pp{1} = 'digraph model {';
pp{2} = 'graph [bgcolor=transparent, pad=0.25]';
pp{3} = 'node [shape=point color=red];';

fType = sum(SM);
n_in_out = size( find(fType ~= 0) ,2);
for i=1:n_in_out
  pp{3+i} = [ 'o' num2str(i) ];
end

pp{4+n_in_out} = 'node [shape=ellipse style=filled fillcolor=lightblue2];';
pp{5+n_in_out} = 'node [fontname=Helvetica];';

w = ones(m,1);
if linlog == 1
  w = widths;
elseif linlog == 2
  w = lwidths;
end
w(find(magnitudes==-1))=-1;
w(find(magnitudes==-2))=-2;
disp([magnitudes w]);

ctr_in_out = 1;
for i=1:m
  if fType(i) == 0
    L = NM{find(SM(:,i)==-1)};
    R = NM{find(SM(:,i)==1)};
    if w(i) == -1
      pp{5+n_in_out+i} = [ L ' -> ' R ' [color="blue" style="dashed" penwidth=2];'];
    elseif w(i) == -2
      pp{5+n_in_out+i} = [ L ' -> ' R ' [color="red" style="dashed" penwidth=2];'];
    else
      pp{5+n_in_out+i} = [ L ' -> ' R ' [penwidth=' num2str(w(i)) '];'];
    end
  elseif fType(i) == 1
    L = [ 'o' num2str(ctr_in_out) ];
    ctr_in_out = ctr_in_out + 1;
    R = NM{find(SM(:,i)==1)};
    if w(i) == -1
      pp{5+n_in_out+i} = [ L ' -> ' R ' [color="blue" style="dashed" penwidth=2];'];
    elseif w(i) == -2
      pp{5+n_in_out+i} = [ L ' -> ' R ' [color="red" style="dashed" penwidth=2];'];
    else
      pp{5+n_in_out+i} = [ L ' -> ' R ' [penwidth=' num2str(w(i)) '];'];
    end
  elseif fType(i) == -1
    L = NM{find(SM(:,i)==-1)};
    R = [ 'o' num2str(ctr_in_out) ];
    ctr_in_out = ctr_in_out + 1;  
    if w(i) == -1
      pp{5+n_in_out+i} = [ L ' -> ' R ' [color="blue" style="dashed" penwidth=2];'];
    elseif w(i) == -2
      pp{5+n_in_out+i} = [ L ' -> ' R ' [color="red" style="dashed" penwidth=2];'];
    else  
      pp{5+n_in_out+i} = [ L ' -> ' R ' [penwidth=' num2str(w(i)) '];'];
    end
  end 
end

pp{6+n_in_out+m} = '}';

% store to file

fid = fopen('model.dot', 'wt');
for i=1:length(pp)
  fprintf( fid, '%s\n', pp{i});
end
fclose(fid);

% compose the legend:

lg{1} = 'digraph {';
lg{2} = 'rankdir=LR';
lg{3} = 'node [shape=plaintext]';
lg{4} = 'graph [splines=ortho]';
lg{5} = 'subgraph cluster_01 { ';
lg{6} = 'label = "Legend"; style=invis;';
lg{7} = 'key [label=<<table border="0" cellpadding="1" cellspacing="0" cellborder="0">';

lims = round(10^Ndigits*limlin)/10^Ndigits;

for i=1:levels
  lg{7+i} = ['<tr><td align="right" port="i' num2str(i) '">' num2str(lims(i)) ' - ' num2str(lims(i+1)) ' </td></tr>'];
end

lg{8+levels} = '</table>>]';
lg{9+levels} = 'key2 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">';
		
for i=1:levels
  lg{9+levels+i} = [' <tr><td port="i' num2str(i) '">&nbsp;</td></tr>'];
end

lg{10+2*levels} = '</table>>]';

for i=1:levels
  lg{10+2*levels+i} = ['key:i' num2str(i) ':e -> key2:i' num2str(i) ':w [penwidth=' num2str(levels+1-i) ']'];
end

lg{11+3*levels} = '} }';

fid2 = fopen('legend.dot', 'wt');
for i=1:length(lg)
  fprintf( fid2, '%s\n', lg{i});
end
fclose(fid2);

!dot -Tps model.dot -o model.ps
!ps2eps -l -B -s b0 -c -n -f model.ps 
!rm model.ps
!dot -Tpng model.dot -o model.png

!dot -Tps legend.dot -o legend.ps
!ps2eps -l -B -s b0 -c -n -f legend.ps 
!rm legend.ps
!dot -Tpng legend.dot -o legend.png

%subplot(2,1,1); hist(widths,levels);
%title('Linear distribution of edge widths')
%subplot(2,1,2); hist(lwidths,levels);
%title('Logarithmic distribution of edge widths')
