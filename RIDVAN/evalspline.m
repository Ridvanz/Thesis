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

yhat=zeros(Nt,1);

G = TN.core;

for i = 1:length(G)
Gsize = [Nt size(G{i},1)  size(G{i},3)];
tempz = reshape(ut{i}*unfold(G{i},2), Gsize);
V{i} = permute(tempz, [2 3 1]);
end


for jj=1:Nt   
f = (V{1});
f=1;
for i = 1:length(G)
f = f*V{i}(:,:,jj);
end
yhat(jj)=f;
end

end

