function [un] = basisvectors(featurez,n,m)

%[un] = basisvectors(featurez,n,m)
% -------------
% Construct basis vectors from data samples
% 
% un        =   B-spline basis vectors
% 
% featurez  = Matrix of training data, rows are samples, columns are different features
% 
% n         = degree of B-spline
% 
% m         = number of knots intervals

[N, d]=size(featurez); 

In= m+n;

bs = bspline([0:n+1]);  %Generate b-spline model
M = flipud(bs.coefs)';  %Get the basis Matrix of size (n+1, n+1)

knotdist = 1/m;         %Get the distance between knots

indexes = floor(featurez/knotdist)+1;       %Calculate in which knot intervals the data samples fall.
indexes(indexes>m)= m;                      %Correction for when a data sample equals the upper limit (one).    

inputs = (featurez/knotdist)-indexes+1;     %Map the datapoints to the knot intervals.

% The complexity of this operation is O(N*d*(n+1)^2) = O(N*d*n^2)
for i=1:d
bn = inputs(:,i).^[n:-1:0]*M;   % Construct the nonzero elements of the b-spline basis vectors using the matrix form.

un{i} = zeros(N,In);
for ii=1:N
   un{i}(ii,indexes(ii,i):indexes(ii,i)+n) = bn(ii,:); %Store them in the correct location within the basis vector. 
end
end

end

