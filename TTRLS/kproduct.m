function [Y, Ys] = kproduct(X,U,n)

% K mode product of tensor with matrix

N = length(size(X));
O = [n setdiff(1:N,n)];
B = permute(X,O);
D = size(B);
Xn = reshape(B, D(1), []);

Yn = U*Xn;

D(1)= size(U,1);
Yr = reshape(Yn, D);
Ys = ipermute(Yr,O);

Y = squeeze(Ys);

end

