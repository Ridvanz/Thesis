function [G, A, e] = hooi(X, A, R, e)

l = length(size(X));
iter=0;

while(true&&iter<10000)
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
error = inproduct(Xdiff,Xdiff);
    
    if (error<e)
        break
    end
    
if (mod(iter,100)==0)    
    disp(iter)
end

end

end

