function [simoutput] = simnarx(TN,input,output,inlags,outlags,n,m)

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

siminput = input;
simoutput = output;
numz = numel(outlags);
beginz = max(outlags(end),inlags(end)) +1;
for j = beginz:length(input)
simfeaturez = [simoutput(j-outlags(end:-1:2))' siminput(j-inlags(1:end))'];
simfeaturez(simfeaturez<0)=0;
simfeaturez(simfeaturez>1)=1;
simoutput(j) = evalspline(TN,simfeaturez,n,m);
end

end

