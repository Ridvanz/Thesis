function P = reshapeP(P)

for j=1:size(P.n,1)
P.core{j} = reshape(P.core{j},P.n(j,:));
end

end

