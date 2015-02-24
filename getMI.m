function I = getMI(X,K,varargin)
    
d=size(X{1},1); %dimension
% K=20; %neighbours

% Volume
V=pi^(d/2)/gamma(d/2+1);

Nq=length(X);
% fX=cell(Nq,Nq);
hRS=nan(1,Nq); 
poolX=[];
for s1=1:Nq
%     for s2=1:Nq            
%         if s1==s2 %same well            
            [~,DistK]=annquery(X{s1},X{s1},K+1);
            Dk=DistK(K+1,:)';

            Nsamp=size(X{s1},2);
            fX=(K)./(Nsamp*V*Dk.^d);                 

%             hRS(s1)=-sum(log2(max(fX,eps)))/length(Dk);
            hRS(s1)=-sum(log2(fX(fX>eps)))/nnz((fX>eps));
            
            poolX=[poolX,X{s1}];

%         else   %different well            
%             [~,DistK]=annquery(X{s2},X{s1},K);
%             Dk=DistK(K,:)';
% 
%             Nsamp=size(X{s2},2);
%             fX{s1,s2}=(K)./(Nsamp*V*Dk.^d);   
%         end
%     end
end
       
[~,DistK]=annquery(poolX,poolX,K+1);
Dk=DistK(K+1,:)';

Ntotal=size(poolX,2);
fX=(K)./(Ntotal*V*Dk.^d);                 

% hR=-sum(log2(max(fX,eps)))/length(Dk);
hR=-sum(log2(fX(fX>eps)))/nnz((fX>eps));

if numel(varargin)>0
    Qstart=varargin{1};
else
    Qstart=ones(Nq,1);
    Qstart=Qstart/sum(Qstart);
    I=hR-hRS*Qstart;  
end
