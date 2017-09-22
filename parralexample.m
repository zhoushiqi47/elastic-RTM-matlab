


%parfor并行方式计算


c2=1;

parfor ii = 1:5

  c2 = ii;
  a(ii)=2*c2

end
a


