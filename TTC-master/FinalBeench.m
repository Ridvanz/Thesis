clear all
close all
clc

load('Benchmark_EEG_medium')
load('Benchmark_EEG_small')

PAR = size(EEGdata,1);          % #participants 
REA = size(EEGdata{2}.angle,1); % #realizations
PER = size(EEGdata{2}.angle,2); % #periods
SAM = size(EEGdata{2}.angle,3); % #samples

rng(1);
%% 
% Preprocessing Data

ds = 8

for i = 1:PAR
    for j = 1:REA
out = squeeze(EEGdata{i}.comp(j,:,:))';
in  = squeeze(EEGdata{i}.angle(j,:,:))';

outlierz = isoutlier(mean(out)) | isoutlier(var(out));

out=smoothdata(out(:,~outlierz),'gaussian',ds);
in=smoothdata(in(:,~outlierz),'gaussian',ds);

outbounds = (max(abs(out))<4);

EEGclean{i}.out{j}=out(1:ds:end,outbounds)/8+0.5;
EEGclean{i}.in{j}=in(1:ds:end,outbounds)/0.12+0.5;
    end
end
%%


% average over periods and downsample, see Section II-B of paper
for Sidx=1:PAR
    smalldata.u(Sidx,:,:)=squeeze(mean(EEGdata{Sidx}.angle(:,:,1:ds:end),2)/0.12+0.5);
    smalldata.y(Sidx,:,:)=squeeze(mean(EEGdata{Sidx}.comp(:,:,1:ds:end),2))/8+0.5;
end
%%

plot(squeeze(smalldata.u(1,1,:)))
plot(EEGclean{1}.in{1})

plot(squeeze(smalldata.y(1,1,:)))
plot(mean(EEGclean{1}.out{1}(:,:)'))

plot((EEGclean{1}.out{1}))
%%
partici = 1 

input = [];
output = [];
    for j = 1:REA
input = [EEGclean{partici}.in{j} input];
output = [EEGclean{partici}.out{j} output];
    end
    
tinput = squeeze(smalldata.u(partici,:,:))';
toutput = squeeze(smalldata.y(partici,:,:))';
%%
inlags=[8 7 6 5];
for l = 1:length(inlags)
u(:,l) = reshape(input(end-inlags(l)-200:end-inlags(l),:), [],1);
end

outlags=[0 1 2 3 4];
for l = 1:length(outlags)
y(:,l) = reshape(output(end-outlags(l)-200:end-outlags(l),:), [],1);
end
% -----------------------

for l = 1:length(inlags)
tu(:,l) = reshape(tinput(end-inlags(l)-200:end-inlags(l),:), [],1);
end

for l = 1:length(outlags)
ty(:,l) = reshape(toutput(end-outlags(l)-200:end-outlags(l),:), [],1);
end

% autocorr(u(:,l),80)
% autocorr(y(:,l),80)
% 
% parcorr(u(:,l),20)
% parcorr(y(:,l),20)

% heatmap(corr([u y]));
%%

trainind = randperm(size(u,1),250000);
testind = setdiff([1:size(u,1)],trainind);

testind = [20000:20200 60000:60200 10000:10200 130000:130200 170000:170200 200000:200200 240000:240200];
trainind = setdiff([1:size(u,1)],testind);

featurez = [u(trainind,:) y(trainind,2:end)];
zeta = y(trainind,1);

tfeaturez = [u(testind,:) y(testind,2:end)];
yt = y(testind,1);

% tfeaturez = [tu ty(:,2:end)];
% yt = ty(:,1);

%%
[N, d]=size(featurez); 

nn = 2;  %degree B-spline
knotintervals = 2;
nnn= knotintervals+nn;

bs = bspline([0:nn+1]);
M = flipud(bs.coefs)';

knotdist = 1/knotintervals;

indexes = floor(featurez/knotdist)+1;
indexes(indexes>knotintervals)= knotintervals;

inputs = (featurez/knotdist)-indexes+1;
for i=1:d
bn = inputs(:,i).^[nn:-1:0]*M;

un{i} = zeros(N,nnn);
for ii=1:N
   un{i}(ii,indexes(ii,i):indexes(ii,i)+nn) = bn(ii,:);
end
end

%%

R = nnn^2;
R=15;
r=R*ones(1,d+1);

% %  r(d-1)=nnn^2;
r(1) = 1;   r(d+1)=1;
r(2)=nnn;  r(d)=nnn;
r(3)=2*nnn;  r(d-1)=2*nnn;
r(4)=3*nnn;  r(d-2)=3*nnn;
% r(5)=3*nnn+1;  r(d-3)=3*nnn+1;
% r(6)=3*nnn+2;  r(d-4)=3*nnn+2;
% r(7)=3*nnn+3;  r(d-5)=3*nnn+3;
r(3)=nnn*3;  r(d-1)=nnn*3;
 r(d-2)=nnn*4;

 r = [1 4 8 12 16 12 8 4 1];
% r = [1 4 8 12 8 4 1];

disp(r)
n=nnn*ones(1,d);

init=cell(1,d);
    for i=1:d
        init{i}.core=randn(r(i),n(i),r(i+1));
        init{i}.core=(init{i}.core/(norm(init{i}.core(:))));
        init{i}.r(1)=r(i);
        init{i}.n = n(i);
        init{i}.r(2)=r(i+1);
    end

%%
                    

for i=1:d
    n(i)=size(un{i},2);
end
  

for i=1:d
    TN.core{i}=init{i}.core;
    TN.n(i,:)=[init{i}.r(1),init{i}.n,init{i}.r(2)];
end

dof = sum(prod(TN.n,2));
Pcount = (nnn)^d;
disp(dof/Pcount)

Vp=cell(1,d);
Vm=cell(1,d);

Vm{1}=ones(N,1);
Vp{d}=ones(N,1);

% initialize right-orthonormal cores with prescribed TN ranks
for i=d:-1:2
    Vp{i-1}=dotkron(Vp{i},un{i})*reshape(permute(TN.core{i},[3 2 1]),[r(i+1)*n(i),r(i)]); 
end

difforder=2;
P = diff(eye(nnn),difforder);
PP = P'*P;
%%
res1=[];
res2=[];
MAXITR = d %ceil(1.5*d);
bb=1;
itr=1;                          % counts number of iterations
ltr=1;                          % flag that checks whether we sweep left to right
sweepindex=1;                   % index that indicates which TT core will be updated


nselect = floor(logspace(log10(dof/N),0,MAXITR)*N);

plot(nselect);

%%
itr=1;  
MAXITR =  d;
nselect=N*ones(MAXITR,1);


%%
gamma=0.0001
lambda = 0.01;
% gamma=1*10^-2;
 tic

% while itr<2 || ((e(itr) < e(itr-1)) && (itr < MAXITR) && e(itr) > THRESHOLD)
while (itr < MAXITR )
%     ---------------------updateTT-----------------;
        
dataselect = randperm((N),nselect(itr));
% dataselect = 1:N;
    %         % first construct the linear subsystem matrix
        A=dotkron(Vm{sweepindex}(dataselect,:),un{sweepindex}(dataselect,:),Vp{sweepindex}(dataselect,:));
%         g=pinv(A'*A)*(A'*y(dataselect,:));
     
%         -------------------difference penalty 

        for j=1:d
        Csize = TN.n(j,:);
        
        Dm = reshape(permute(TN.core{j}, [2 1 3]), [Csize(2) Csize(1)*Csize(3)]);
        mDDm = reshape(Dm'*Dm, [Csize(1) Csize(3) Csize(1) Csize(3)]);
        DD{j} = reshape(permute(mDDm,[1 3 2 4]), [Csize(1)*Csize(1) Csize(3)*Csize(3)]);
        PD = P*Dm;
        DPPD = reshape(PD'*PD, [Csize(1) Csize(3) Csize(1) Csize(3)]);
        DWD{j} = reshape(permute(DPPD,[1 3 2 4]), [Csize(1)*Csize(1) Csize(3)*Csize(3)]);
        eyez{j}= reshape(eye(Csize(1)), [Csize(1)^2 1]);
        eyep{j}= reshape(eye(Csize(2)), [Csize(2)^2 1]);
        end
        eyez{d+1}=1;
        
        for p= 1:d
            Dsize = TN.n(sweepindex,:);   
                if sweepindex==p
                    D1 = eyez{sweepindex};
                    D2 = PP(:);
                    D3 = eyez{sweepindex+1};
                elseif sweepindex<p
                    D1 = eyez{sweepindex};  
                    D2 = eyep{sweepindex};   
                    D3= DWD{p}*eyez{p+1};           
                    for it=(p-1):-1:(sweepindex+1)               
                        D3 = DD{it}*D3;
                    end                               
                elseif sweepindex>p   
                    D1= eyez{p}'*DWD{p};
                    for it=(p+1):(sweepindex-1)               
                        D1 = D1*DD{it};
                    end
                    D1=D1';
                    D2 = eyep{sweepindex};   
                    D3= eyez{sweepindex+1}; 
                end
          
            WW = kron(D3, kron(D2, D1));
            Wtemp = permute(reshape(WW, [Dsize(1) Dsize(1) Dsize(2) Dsize(2) Dsize(3) Dsize(3)]), [1 3 5 2 4 6]);
            W{p} = reshape(Wtemp, [prod(Dsize) prod(Dsize)]); 
        end
        WWW = W{1};
        for s =2:d
            WWW = WWW + W{s};
        end
        
        Wopslag{itr} = W;
%         ----------------------- 
%         D = eye(r(sweepindex)*(n(sweepindex))*r(sweepindex+1));         
        AA = A'*A;
        bb = bb+1;
        
      lambda = sqrt(trace(AA)/trace(WWW));

       
        Tracez(bb) = trace(AA);
         Wcez(bb) = trace(WWW);
         
%          tic
        g=pinv(AA + lambda*WWW)*(A'*zeta(dataselect,:));
%         mister(bb) = toc;
%         tic
%         [g,~] = pcg((AA + lambda*WWW),(A'*zeta(dataselect,:)),10^-6,100,[],[],TN.core{sweepindex}(:));
        
%         sister(bb) = toc;
        
        
        if ltr
            % left-to-right sweep, generate left orthogonal cores and update vk1
%             TN.core{sweepindex}=reshape(g,[r(sweepindex),n(sweepindex),r(sweepindex+1)]);
            [Q,R]=qr(reshape(g,[r(sweepindex)*(n(sweepindex)),r(sweepindex+1)])); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex+1)),[r(sweepindex),n(sweepindex),r(sweepindex+1)]);
            TN.core{sweepindex+1}=reshape(R(1:r(sweepindex+1),:)*reshape(TN.core{sweepindex+1},[r(sweepindex+1),(n(sweepindex+1))*r(sweepindex+2)]),[r(sweepindex+1),n(sweepindex+1),r(sweepindex+2)]);
            Vm{sweepindex+1}=dotkron(Vm{sweepindex},un{sweepindex})*reshape(TN.core{sweepindex},[r(sweepindex)*n(sweepindex),r(sweepindex+1)]); % N x r_{i}
        else
            % right-to-left sweep, generate right orthogonal cores and update vk2
%             TN.core{sweepindex}=reshape(g,[r(sweepindex),n(sweepindex),r(sweepindex+1)]);
            [Q,R]=qr(reshape(g,[r(sweepindex),(n(sweepindex))*r(sweepindex+1)])'); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex))',[r(sweepindex),n(sweepindex),r(sweepindex+1)]);
            TN.core{sweepindex-1}=reshape(reshape(TN.core{sweepindex-1},[r(sweepindex-1)*(n(sweepindex-1)),r(sweepindex)])*R(1:r(sweepindex),:)',[r(sweepindex-1),n(sweepindex-1),r(sweepindex)]);
            Vp{sweepindex-1}=dotkron(Vp{sweepindex},un{sweepindex})*reshape(permute(TN.core{sweepindex},[3 2 1]),[r(sweepindex+1)*n(sweepindex),r(sweepindex)]); % N x r_{i-1}    

        end
    
   %     ---------------------updatesweep-----------------; 
     
                if ltr
                    sweepindex=sweepindex+1;
                    if sweepindex== d                
                        ltr=0;
                        timer=toc;
                    end
                else
                    sweepindex=sweepindex-1;
                    if sweepindex== 1                
                        ltr=1;
                        timer=toc;
                    end
                end
    
    
        % only check residual after 1 half sweep
        if (sweepindex==d) || (sweepindex==1) % half a sweep
            
            res1(itr)=norm(A*g-zeta(dataselect,:))^2/(nselect(itr)); % check residual
            res2(itr)=1*(g'*g);
            res3(itr)=1*(g'*WWW*g);
           
            if (itr>4)
%                 if (res1(itr)<1) || (abs(res1(itr)-res1(itr-2))<1)
%                  disp(res1(itr))
%                  break;
%                  end
            end
            
            itr=itr+1; %update iteration
            
            disp(["iteration:" itr timer])
        end   
end  
%%
plot(res1)
plot(res2)
plot(res3)

plot(Wcez)
plot(Tracez)
plot(Tracez./Wcez)

for k = 1:d
    for cud = 1:MAXITR-1
plottester(k,cud) = trace(Wopslag{cud}{k});
end
end
plot(plottester')

% axis([-inf inf 0 100])

figure;
hold on
tic
 plot(pinv(AA + lambda*WWW)*(A'*zeta(dataselect,:)));
   toc; tic 
 plot(pcg(   (AA + lambda*WWW),(A'*zeta(dataselect,:)),10^-6,1000));
toc; tic 
 plot(cgs(   (AA + lambda*WWW),(A'*zeta(dataselect,:)),10^-6,1000));
toc
hold off

figure;
hold on
% plot(mister)
% plot(sister)
hold off
%% 
% Evaluate Test data

[Nt, dt]=size(tfeaturez);

indexes = floor(tfeaturez/knotdist)+1;
indexes(indexes>knotintervals)= knotintervals;

inputs = (tfeaturez/knotdist)-indexes+1;

for i=1:dt
bn = inputs(:,i).^[nn:-1:0]*M;

ut{i} = zeros(Nt,nnn);
for ii=1:Nt
   ut{i}(ii,indexes(ii,i):indexes(ii,i)+nn) = bn(ii,:);
end
end
%%
% [N, d]=size(featurez);

yhat=zeros(Nt,1);

G = TN.core;

for i = 1:length(G)
Gsize = [Nt size(G{i},1)  size(G{i},3)];
tempz = reshape(ut{i}*unfold(G{i},2), Gsize);
V{i} = permute(tempz, [2 3 1]);
end


for jj=1:Nt   
f = (V{1});
f=1;
for i = 1:length(G)
f = f*V{i}(:,:,jj);
end
yhat(jj)=f;
end

%%


erboi = (yhat-yt);
histogram((erboi)')
VAF = 1-var(erboi)/var(yt)
MSE = immse(yhat,yt)

figure
hold on
plot(yt)
plot(yhat)
hold off
axis([0 200 -inf inf])

[sortout,Iss]=sort(yhat);
figure
hold on
plot(yt(Iss))
plot(sortout,'Linewidth',3)
hold off


%%
meanVAF = mean(VAF)



% plot(erboi)
pererror = sum(abs(erboi))/Nt
relerror = sum(abs(erboi)/sum(abs(yt)))

FY = mag2db(abs(fft(yhat)));
FZ = mag2db(abs(fft(yt)));

figure;
hold on
plot(FY(1:floor(end/2)))
plot(FZ(1:floor(end/2)))
hold off

% ploterrcorr(erboi')
% autocorr(erboi')


%%

% smoothed= smoothdata(output,'sgolay',20); 
% figure;
% hold on
% plot(output)
% plot(refer)
% % plot(smoothed)
% hold off

% Mdl = arima(3,4,3)
% EstMdl = estimate(Mdl,zeta)
%%
% for i = 1:90
% testers(i,:) = forecast(EstMdl,3,y1t(i+50:i+140));
% end
%%
% 
% figure
% hold on
% plot(y1t)
% plot(51:142,testers(:,3))
% hold off