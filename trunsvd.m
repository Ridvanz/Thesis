function [Utrun,Strun,Vtrun,r] = trunsvd(C,e)

[U,S,V] = svd(C);
sings = diag(S).^2;
cumuls = cumsum(sings)/sum(sings);
r = find(cumuls>1-e, 1);

Utrun = U(:,1:r);
Strun = S(1:r,1:r);
Vtrun = V(:,1:r);

end