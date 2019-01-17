function grad = gradidi(a,n)

batchsize = length(a);
grad = [zeros(batchsize,1) (a.^[0:n-1]).*[1:n]];

end