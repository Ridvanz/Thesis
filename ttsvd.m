function [G, rk] = ttsvd(A,e)

% Tensor Train decomposition with svd, given a tolerance e

Sizes = size(A);
d = ndims(A);
delta = (e/(sqrt(d-1)));

Atest = inproduct(A,A);
rk = zeros(d,1);
G  = cell(d,1);
rk(1)=1;

C=A;

for i = 1:d-1
    
C = reshape(C,[rk(i)* Sizes(i) numel(C)/(rk(i)* Sizes(i))]);
[U,S,V, rk(i+1)] = trunsvd(C,delta, Atest);
G{i}= reshape(U, [rk(i), Sizes(i), rk(i+1)]);
C = S*V';

end

G{d} = C;

end