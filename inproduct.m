function S=inproduct(A,B)
    
% C = dot(A,B);
% S = sum(C,'all');

S = A(:)'*B(:);

end