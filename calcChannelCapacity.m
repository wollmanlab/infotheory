function [Imax,qopt,Ijack] = calcChannelCapacity(Esb,varargin)

arg.jackknife = false;
arg.jackreplicas = 1; 
arg.festim = {}; 
arg.qopt = []; 
arg.k=10; 
arg = parseVarargin(varargin,arg);

%% jackknife if needed
if arg.jackknife && isempty(arg.qopt)
    jackProb = repmat(linspace(0.6,0.95,20),1,arg.jackreplicas);
    Ijack = nan(size(jackProb));
    parfor i=1:numel(jackProb)
        IX = cellfun(@(e) randsample(size(e,2),ceil(jackProb(i)*size(e,2))),Esb,'uniformoutput',0);
        Ejack = cellfun(@(e,ix) e(:,ix),Esb,IX,'uniformoutput',0);
        Ijack(i) = calcChannelCapacity(Ejack);
    end
    b = regress(Ijack(:),[ones(numel(jackProb),1) 1./jackProb(:)]);
    Imax = b(1); 
    qopt = []; 
    return
end

%% get the probability densities if they were not supplied
if isempty(arg.festim)
    Nq = numel(Esb);
    Festim = cell(Nq);
    for i=1:Nq
        Xi = Esb{i};
        parfor j=1:Nq
            if j~=i
                Festim{i,j} = kNNdensityEstimateUsingAnotherWell(Esb{j},Xi);
            else
                Festim{i,j} = kNNdensityEstimateUsingAnotherWell(Esb{j},Xi);
            end
        end
    end
else
    % arg.festim could be matrix Nq x Nq of estimation or a row vectpr of
    % prob distribution (or gmm) objects with the method pdf
    Nq=size(arg.festim,2); 
    if isa(arg.festim{1},'ProbDist') || isa(arg.festim{1},'gmdistribution')
        Festim = cell(Nq); 
        for i=1:Nq
            for j=1:Nq
                Festim{i,j} = arg.festim{j}.pdf(Esb{i}); 
            end
        end
        
    else % assume that arg.festim are actual densities
        Festim = arg.festim;
    end
end

%% Optimize Channel capacity

fun = @(q) -fMIfromDensitiesWithoutSampling(Festim,q);
opts = optimoptions(@fmincon,'Display','Iter');
q0 = ones(Nq,1)/Nq;

if ~isempty(arg.qopt)
    qopt = arg.qopt(:); 
else
    qopt= fmincon(fun,q0,[],[],ones(1,Nq),1,zeros(Nq,1),ones(Nq,1),[],opts);
end
Imax = -fun(qopt);

% try
%     opts = psoptimset('Display','Iter');
%     qopt= patternsearch(fun,q0,[],[],ones(1,Nq),1,zeros(Nq,1),ones(Nq,1),[],opts);
% catch e
%     warning('couldn''t run patternsearch (error: %s), using fminsearchcon instead',e.message)
%     opts2 = optimset('Display','Iter');
%     qopt = fminsearchcon(fun,q0,zeros(Nq,1),ones(Nq,1),ones(1,Nq),1,[],opts2);
% end
% 
% Imax = -fun(qopt);