function [result,FV,YY,SV,Wv]=PFSC(X,s,alpha,beta)
% X is the data
% s is the true class label.


B=cell(size(X));
for ii=1:size(X,1)
    B{ii}= inv((X{ii}*X{ii}')+beta*eye(size(X{1},1)));
end

%initialize parameter

[v1,v2]=size(X);
FV=cell(v1,v2);
SV=cell(v1,v2);
Wv=1/v1*ones(v1,1);
P=ones(v1,v1);
Q=ones(v1,1);
G=ones(v1,1);
LV=cell(v1,1);
c=length(unique(s));
n=size(X{1},1);
Fv=randn(n,c);
Fv= orth(Fv);
 Sv=eye(n);


Aeq=ones(1,v1);
Beq=1;
lb=zeros(v1,1);
    
for num = 1:v1
  FV{num}=Fv;
  SV{num}=Sv;
end
% YY = rand(n,c);
% YY = orth(YY);
YY=FV{1};
% begin iteration
for i=1:200

    Yold=YY;
%update SV
    for num = 1:v1
        Sv=SV{num};
        xv=X{num};
        xv=xv';
        Fv=FV{num};
        %update Sv
        parfor ij=1:n
            d=distance(Fv,n,ij);
            xx=xv'*xv
            Sv(:,ij)=B{num}*(xx(ij,:)' - alpha/4*d'); 
        end

%        Sv(find(Sv<0))=0;
       Sv=Sv-diag(diag(Sv));
        Sv=(Sv+Sv')/2;
        D = diag(sum(Sv));
        Lv = D-Sv;
        SumFj=zeros(n,n);
        for nn=1:v1
            Fj=FV{nn};
            if(nn ~= num)
              SumFj=SumFj+(Fj*Fj');
            end
        end
        % calculate Fv 
        Mv=alpha*Lv-2*(YY*YY')+Wv(num)*eye(n)-2*SumFj;
        [Fv, temp, ev]=eig1(Mv, c, 0); 
        
        FV{num}=Fv;
        SV{num}=Sv;
        LV{num}=Lv;
    end
  %calculate P(i,j)
    for ii=1:v1
        for jj=1:v1
            Fii=FV{ii};
            Fjj=FV{jj};
            P(ii,jj) = 2*trace(Fii*Fii'*(Fjj*Fjj'));
        end
    end
    %calculate Qi
    for ii=1:v1 
        G(ii) =(norm(X{ii}'-X{ii}'*SV{ii},'fro'))^2+alpha * trace(FV{ii}'*LV{ii}*FV{ii}+beta*trace(SV{ii}'*SV{ii}));
        Q(ii)=-1*G(ii)+ 2*trace(YY*YY'*FV{ii}*FV{ii}');
        Q(ii)=-1*Q(ii);
    end
    %update Wv
    Wv= quadprog(P,Q,[],[],Aeq,Beq,lb);
    
    %update Y
    A=zeros(n);
    for num=1:v1
        Fv=FV{num};
        A=2*Wv(num)*(Fv*Fv')+A; 
    end
    AA = eye(n)-A;
    [YY, temp, ev]=eig1(AA, c, 0);

   if i>5 &&((norm(YY-Yold)/norm(Yold))<1e-5)
       break
   end 
end



for ij=1:10
actual_ids= kmeans(YY, c, 'emptyaction', 'singleton', 'replicates', 1, 'display', 'off');
[res(ij,:)] = Clustering8Measure(actual_ids,s);
end


result(1,1)=mean(res(:,1));result(1,2)=std(res(:,1));
result(2,1)=mean(res(:,2));result(2,2)=std(res(:,2));
result(3,1)=mean(res(:,3));result(3,2)=std(res(:,3));
result(4,1)=mean(res(:,4));result(4,2)=std(res(:,4));
result(5,1)=mean(res(:,5));result(5,2)=std(res(:,5));
result(6,1)=mean(res(:,6));result(6,2)=std(res(:,6));
end
function [all]=distance(F,n,ij)
  for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
  end
end   
       
    
    
        
        