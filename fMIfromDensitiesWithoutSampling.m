function [I,Hr_s,Hr] = fMIfromDensitiesWithoutSampling(Fcond,q)
% Fcond is be a N x N cell array of conditional probability densities
%       the cell i,j has the probability to get the observations seen in i given 
%       stimulus j. 
% q is the probability for inputs

Nq = numel(q); 
assert(size(Fcond,1)==size(Fcond,2),'F should be square')
assert(size(Fcond,1)==Nq,'Fcond should be size of Nq x Nq'); 

%% First calculate the non-conditional probability densitites
Q = repmat(q(:)',Nq,1); 
FtimesQ = cellfun(@(f,q) f(:)*q,Fcond,num2cell(Q),'uniformoutput',0); 

F = cell(Nq,1); 
for i=1:Nq
    F{i} = sum(cat(2,FtimesQ{i,:}),2); 
end
Fself = Fcond(eye(Nq)==1); 

%% calculate Entropies for conditional case (easy)
fH = @(F) -sum(log2(F(F>eps)))/nnz(F>eps); 
Hr_s = cellfun(fH,Fself);


Wq = cellfun(@(q,f) ones(size(f))*q,num2cell(q),F,'uniformoutput',0); 
Wq = cat(1,Wq{:}); 
Fall = cat(1,F{:});
Fall(Fall<=eps)=nan; 


Wq = Wq./nansum(Wq); 
Hr = -nansum(log2(Fall).*Wq); 

I = Hr-sum(q.*Hr_s);


%% calcualte entropies for complete. 
% We have the acurate densities, F but since we have a biased sample
% of them, we cann't use the actual sample we have. But we could bootstrap
% ourselves from this mess as follows: 


% 
% SZ  = 3000; 
% Iter = 1000; 
% % Hbt = nan(Iter,1); 
% cnt=0; 
% for iq = 1:numel(q)
%     cnt=cnt+ceil(SZ*q(iq)); 
% end
% Fall2 = nan(cnt,Iter);  
% for it = 1:Iter
%     Fall = cell(numel(q),1); 
%     for iq=1:numel(q)
%         ix = randi(numel(F{iq}),ceil(SZ*q(iq)),1); 
% %         ix = randsample(numel(F{iq}),ceil(SZ*q(iq)),true); 
%         Fall{iq}=F{iq}(ix); 
%     end
%     Fall2(:,it) = cat(1,Fall{:});
% %     Fall = cat(1,Fall{:});
% %     Hbt(it) = fH(Fall);
% end
% Fall2(Fall2<eps)=nan; 
% % Hr = mean(nanmean(-log2(Fall2))); 
% Hr = nanmean(-log2(Fall2(:))); 
% % Hr = mean(Hbt); 
% I = Hr-sum(q.*Hr_s);
% 
% % %% delete me later
% Fall = cat(1,F{:});
% Hr = fH(Fall); 
% I = Hr-sum(q.*Hr_s);
