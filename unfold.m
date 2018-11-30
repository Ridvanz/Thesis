function [U, O] = unfold(A,n)

N = length(size(A));
O = [n setdiff(1:N,n)];
B = permute(A,O);
d = size(B,1);
U = reshape(B, d, []);
end
