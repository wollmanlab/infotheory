function f = kNNdensityEstimateUsingAnotherWell(Xw,Xq,varargin)

%% set input arguments; 
arg.k = 10; 
arg.kdtree = []; 
arg = parseVarargin(varargin,arg); 
k = arg.k; 

[r,T]=size(Xw);

%% deal with the case that Xq is empty
if ~exist('Xq','var') || isempty(Xq) || isequal(Xw,Xq)
    Xq=Xw; 
    k=k+1; % to remove the selfie case
end

if isempty(arg.kdtree)
    % to speed things a bit can use annquery => but mex file can cause Matlab
    % to crash
    [~,D]=annquery(Xw,Xq,k);
elseif isa(arg.kdtree,'KDTreeSearcher')
    [~,D] = knnsearch(arg.kdtree,Xq','K',k); 
else
    error('If you suppply a kdtree it must be of KDTreeSearcher class!'); 
end

D=D(k,:);

% [~,D] = knnsearch(Xw',Xq','K',k); 
% D=D(:,k); 

NDsphere = pi.^(r/2)./gamma(r/2+1); 
f = (k-1)./D.^r/T/NDsphere; 
f=f(:); 
