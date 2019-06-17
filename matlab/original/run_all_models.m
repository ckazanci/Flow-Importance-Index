models = {'_3comp','_oyster','_lake_findley','_silver_springs','_cone_spring','_marion_lake','_lake_wingra','_English_channel','_somme_estuary',''}
Nmodels = size(models,2)

for i=1:Nmodels
  disp(models{i});
  run([ 'EcoNet' models{i} ]);
  meas = importance_measure(SM);
  model2dots(meas(:,1),7,NM,SM,1)
  system(['mv model.ps Scaled' models{i} '_imp1.ps']);
  system(['mv model.png Scaled' models{i} '_imp1.png']);
  meas = importance_measure(SM);
  model2dots(meas(:,2),7,NM,SM,1);
  system(['mv model.ps Scaled' models{i} '_imp2.ps']);
  system(['mv model.png Scaled' models{i} '_imp2.png']);
end
