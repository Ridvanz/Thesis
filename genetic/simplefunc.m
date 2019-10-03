function fitness = simplefunc(chrom, objects, points, par, W)

warning('off','all')

L=chrom(1:3,:);

lum=zeros(1,size(points,2));
totalcost=0;

% ------------------------------------------------------

for l=1:size(L,2)
    
    lamp = L(1:2,l);
    power = L(3,l);
% ----------------------------------------
for k=1:size(points,2)
    
    for i=1:size(objects,2)
o1 = objects{i}.p1;
o2 = objects{i}.p2;
    
b = lamp-o1 ;
A = [lamp-points(:,k) o2-o1] ;

T(:,k,i) = A\b;

     end
end
% ----------------------------------------

isittrue = sum((0<=T) & (T<=1),1)==2;
M = max(isittrue,[],3);
indexes = find(~M);
lightpoints=  points(:,indexes);

dist = sum((lightpoints-lamp).^2);

illu= power * exp(-((dist)/(par.h)));
lum(indexes)= lum(indexes) + illu; 

cost = par.price + par.powercost * power; %cost per year
totalcost = totalcost + cost;

end

fitness =  (1-W)*sum((lum-par.targetlum).^2) +  W * totalcost;

end
