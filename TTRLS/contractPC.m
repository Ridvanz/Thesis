function P = contractPC(P,C,mode)

% contractPC
%mode = 0 for right mode contraction, 1 for left

P.n(:,3)=[];

dd= size(P.n,1);

Psize=P.n(:,[1 2 end]);

for i = 1:dd
ind=3-(P.n(i,1)==1) - mode;
P.core{i} =reshape(C(i,:)*unfold(squeeze(P.core{i}),ind),Psize(i,:));
end


end

