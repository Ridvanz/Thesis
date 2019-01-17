function [C] = contract(A,B, d1, d2)

% Contract tensors A and B along dimension d1 and d2


[U1, O1] = unfold(A,d1);

[U2, O2] = unfold(B,d2); 

Asize = size(A);
Bsize = size(B);

temp = reshape(U1'*U2,[Asize(O1(2:end)) Bsize(O2(2:end))]);

C = squeeze(temp);


end

