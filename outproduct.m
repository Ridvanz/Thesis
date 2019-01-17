function K=outproduct(varargin)

%  Outer product of an arbitraty number of vectors

if nargin==0
    error('Please input vectors');
else
    
K= 1;

sizes = zeros(size(varargin));

for i = 1:length(varargin)
vector = varargin{i};

if ~iscolumn(vector)
vector = vector';
end
    
sizes(i)=length(vector);
K = kron(vector, K);
end

K = reshape(K,sizes);

end
end