function yhat = evalspline(TN,tfeaturez,n,m)

[Nt, dt]=size(tfeaturez);

In= m+n;

bs = bspline([0:n+1]);
M = flipud(bs.coefs)';

knotdist = 1/m;

indexes = floor(tfeaturez/knotdist)+1;
indexes(indexes>m)= m;

inputs = (tfeaturez/knotdist)-indexes+1;

for i=1:dt
bn = inputs(:,i).^[n:-1:0]*M;

ut{i} = zeros(Nt,In);
for ii=1:Nt
   ut{i}(ii,indexes(ii,i):indexes(ii,i)+n) = bn(ii,:);
end
end

yhat = ones(Nt,1);
for i=1:size(TN.core,2)
    yhat=dotkron(yhat,ut{i})*reshape(TN.core{i},[TN.sz(i,1)*TN.sz(i,2),TN.sz(i,3)]); 
end


end

