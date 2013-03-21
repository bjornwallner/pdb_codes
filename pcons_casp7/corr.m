[target,casp,new]=textread('tomat','%s %f %f');
targets=textread('target_list','%s');
C=[];
for i=1:length(targets)
  t=find(strcmp(target,targets{i}));
  
  tmp=corrcoef([casp(t) new(t)]);
  C=[C;tmp(2)];
  
  
end
