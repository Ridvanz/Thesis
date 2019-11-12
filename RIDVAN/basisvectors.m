function [un] = basisvectors(featurez,n,m)

[N, d]=size(featurez); 

In= m+n;

bs = bspline([0:n+1]);
M = flipud(bs.coefs)';

knotdist = 1/m;

indexes = floor(featurez/knotdist)+1;
indexes(indexes>m)= m;

inputs = (featurez/knotdist)-indexes+1;
for i=1:d
bn = inputs(:,i).^[n:-1:0]*M;

un{i} = zeros(N,In);
for ii=1:N
   un{i}(ii,indexes(ii,i):indexes(ii,i)+n) = bn(ii,:);
end
end

end

