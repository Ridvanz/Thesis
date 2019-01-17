function [G, A, R] = hosvd(X)

% High Order SVD, Tucker decomposition

l = length(size(X));
A = cell(l,1);
G = X;

for n = 1:l
    Xu = unfold(X,n);
    R(n)= rank(Xu);
    [U,~,~] = svd(Xu); 
    A{n} = U(:,1:R(n));
[~, G] = kproduct(G,A{n}',n);
end

end

