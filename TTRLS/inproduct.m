function S=inproduct(A,B)
    
%  Inproduct of two tensors with equal dimensions

% C = dot(A,B);
% S = sum(C,'all');

S = A(:)'*B(:);

end