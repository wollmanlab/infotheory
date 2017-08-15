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

%% create weights for a weighted average of the contribution of the non-conditional densities
%Wq = cellfun(@(q,f) ones(size(f))*q,num2cell(q),F,'uniformoutput',0); 
Wq = cellfun(@(q,f) ones(size(f))*q/length(ones(size(f))),num2cell(q),F,'uniformoutput',0); %LeeLab_SG_line modified to match supplement
Wq = cat(1,Wq{:}); 
Fall = cat(1,F{:});
Fall(Fall<=eps)=nan; 


%Wq = Wq./nansum(Wq); % LeeLab_SG_Commented out original line
Hr = -nansum(log2(Fall).*Wq); 

I = Hr-sum(q.*Hr_s);
