function [optima,MORDOR] = findorder(Z)

sz= size(Z);
Bgy =ceil(max(Z,[],'all'))+1;
MORDOR=[];

BAD = Z+Bgy*eye(sz(1));
for r = 1:sz(1)
    
    DAD = BAD;
    row = r;
    ORDOR = r;
    
for j=1:sz(1)-1
[M,col] = min(DAD(row,:));

distas(r,j)=M;
ORDOR = [ORDOR  col];

DAD(row,:) = Bgy;
DAD(:,row) = Bgy;
row=col;
end

MORDOR = [MORDOR ; ORDOR];
end

[~,bread] = min(sum(distas,2));
optima = MORDOR(bread,:)';

end

