function [X] = tcontract(G)

% Input is a cell array of tensor cores G
% does not work yet when there are TT ranks of 1

G1 = permute(G{1},[3 2 1]);
Gsizes= [1 size(G1)];    
Left = reshape(G1, [Gsizes(1)*Gsizes(2) Gsizes(3)])';

i_order = [Gsizes(3)];

for i = 2:length(G)

Gi = permute(G{i},[3 2 1]);
Gsizes= size(Gi);    
G_unfold = reshape(Gi, [Gsizes(1)*Gsizes(2) Gsizes(3)])';
    
contract = Left * G_unfold;

i_order = [i_order Gsizes(2)]; 
order = [i_order(1:i-1) Gsizes(1) i_order(end)];

i_permute = [1:i+1];
i_permute([end-1 end])= i_permute([end end-1]);

resh = reshape(contract, order);
Left = reshape(permute(resh, i_permute), [], Gsizes(1));

end

X = squeeze(resh);

end

