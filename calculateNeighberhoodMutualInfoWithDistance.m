function [MI,MIvec,dSpacing] = calculateNeighberhoodMutualInfoWithDistance(Mx,D,varargin)

%% parse varargin
arg.n = 2:25; 
arg.sz=1000;
arg.woundedge=225; 
arg.dmx = 800;
arg.cellradius = 15;
arg = parseVarargin(varargin,arg); 
%%
N=arg.n; % number of gridpoint to divide space into
dSpacing = nan(size(N)); 
SZ=arg.sz; 
WoundEdge=arg.woundedge; 
Dmx = arg.dmx; 
cellRadius = arg.cellradius;
% resthreshhold = 0.0175; 
% 
%% resample cells into neighberhoods 
Ncell = nan(N,1); 
Grd = cell(N,1); 
MIvec = nan(N,1); 
for iN = 1:numel(N)
    %%
    dGrid = linspace(WoundEdge,Dmx,N(iN)); 
    dSpacing(iN) = mean(diff(dGrid)); 
    Grd{iN}=dGrid(1:end-1)+diff(dGrid)/2;
    Ncell(iN) = ceil((dSpacing(iN)/cellRadius).^2); 
    M=cell(numel(dGrid)-1,1); 
    Qs = nan(size(M)); 
    for i=2:numel(dGrid)
        ixOrg = find(D>dGrid(i-1) & D<=dGrid(i));
        Qs(i-1)=numel(ixOrg); 
        ix = reshape(randsample(ixOrg,SZ*Ncell(iN),true),Ncell(iN),SZ);
        M{i-1}=Mx(ix); 
    end
    Qs=Qs./sum(Qs);
    M = cellfun(@mean,M,'uniformoutput',0);
    MIvec(iN) = calcChannelCapacity(M,'qopt',Qs(:)');
end
MI=max(MIvec); 

%%
