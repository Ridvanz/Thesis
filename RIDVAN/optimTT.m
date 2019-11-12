function [TN,Vm,Vp,res1,res2] = optimTT(TN,Vm,Vp,un,zeta,MAXITR,nselect,lambda,difforder)

res1=[];
res2=[];
itr=1;                          % counts number of iterations
ltr=1;                          % flag that checks whether we sweep left to right
sweepindex=1;                   % index that indicates which TT core will be updated

[N, ~]=size(un{1}); 
d=size(TN.sz,(1)); 
I = TN.sz(:,2);
r = [TN.sz(:,1);1];
P = diff(eye(I(1)),difforder);
PP = P'*P;

tic
while (itr < MAXITR )
%     ---------------------updateTT-----------------;
        %Select random batch of data
        dataselect = randperm((N),nselect(itr));  
        
        %Construct regressor A
        A=dotkron(Vm{sweepindex}(dataselect,:),un{sweepindex}(dataselect,:),Vp{sweepindex}(dataselect,:));

        %Construct difference penalty matrices

        for j=1:d
        Csize = TN.sz(j,:);
        
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
            Dsize = TN.sz(sweepindex,:);   
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
        
        AA = A'*A;
        
        WWW = trace(AA)/trace(W{1})*W{1};
        for s =2:d
            WWW = WWW + trace(AA)/trace(W{s})*W{s};
        end
        
%       D = eye(r(sweepindex)*(I(sweepindex))*r(sweepindex+1));         
        

        %Solve linear subsystem 
%         if itr>MAXITR-1
             g=pinv(AA + lambda*WWW)*(A'*zeta(dataselect,:));
%         else
%             [g,~] = pcg((AA + lambda*WWW),(A'*zeta(dataselect,:)),1/nselect(itr),1000,[],[],TN.core{sweepindex}(:));
%         end
        
        if ltr
%             left-to-right sweep, generate left orthogonal cores and update Vm
            [Q,R]=qr(reshape(g,[r(sweepindex)*(I(sweepindex)),r(sweepindex+1)])); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex+1)),[r(sweepindex),I(sweepindex),r(sweepindex+1)]);
            TN.core{sweepindex+1}=reshape(R(1:r(sweepindex+1),:)*reshape(TN.core{sweepindex+1},[r(sweepindex+1),(I(sweepindex+1))*r(sweepindex+2)]),[r(sweepindex+1),I(sweepindex+1),r(sweepindex+2)]);
            Vm{sweepindex+1}=dotkron(Vm{sweepindex},un{sweepindex})*reshape(TN.core{sweepindex},[r(sweepindex)*I(sweepindex),r(sweepindex+1)]); 
        else
%             right-to-left sweep, generate right orthogonal cores and update Vp
            [Q,R]=qr(reshape(g,[r(sweepindex),(I(sweepindex))*r(sweepindex+1)])'); 
            TN.core{sweepindex}=reshape(Q(:,1:r(sweepindex))',[r(sweepindex),I(sweepindex),r(sweepindex+1)]);
            TN.core{sweepindex-1}=reshape(reshape(TN.core{sweepindex-1},[r(sweepindex-1)*(I(sweepindex-1)),r(sweepindex)])*R(1:r(sweepindex),:)',[r(sweepindex-1),I(sweepindex-1),r(sweepindex)]);
            Vp{sweepindex-1}=dotkron(Vp{sweepindex},un{sweepindex})*reshape(permute(TN.core{sweepindex},[3 2 1]),[r(sweepindex+1)*I(sweepindex),r(sweepindex)]);  
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
            res2(itr)=1*(g'*WWW*g);
           
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



end

