function f = contractX(G,C)

d= size(G.n,1);

% if G.n(1,3)>1
    Gsize = G.n(:,[1 end]);
% else
%     Gsize = G.n(:,[1 end]);
% end


% [G.n(i,1) G.n(i,end)]; 

for i = 1:d
Gt.core{i} = reshape(C(i,:)*unfold(G.core{i},2), Gsize(i,:));
end

f = Gt.core{1};
for j = 2:d
% f = reshape(f,1, []) *squeeze(G.core{j});
f = f*(Gt.core{j});

end

end

