function S=inproduct(A,B)

if nargin==0
    error('Please input vectors');
else
    
% C = dot(A,B);
% S = sum(C,'all');

S = A(:)'*B(:);

end
end