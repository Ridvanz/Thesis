function [G, A, error1, Xr] = hooi(X, A, R, e, maxiter)

l = length(size(X));
iter=0;
error1= 1;
error2= 0;

while(iter<maxiter && abs(error1-error2)>0.001)
iter=iter+1;
    
G=X;
    
for i = 1:l

o = [setdiff(1:l,i)];
    
Y = X;

    for n = o

    [~, Y] = kproduct(Y,A{n}',n);

    end

Yu = unfold(Y,i);
[U,~,~] = svd(Yu);
A{i} = U(:,1:R(i));

[~, G] = kproduct(G, A{i}',i);

end

Xr=G;

for j = 1:l

[~, Xr] = kproduct(Xr,A{j},j);

end

Xdiff = X-Xr;

error2=error1;
error1 = inproduct(Xdiff,Xdiff);
    
    if (error1<e)
        disp(iter)
        disp(error1)
        break
    end
    
if (mod(iter,100)==0)    
    disp(iter)
    disp(error1)
end

end

end

