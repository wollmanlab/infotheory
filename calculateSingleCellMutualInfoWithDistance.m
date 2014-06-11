function [MI,MIvec,dGrid] = calculateSingleCellMutualInfoWithDistance(Mx,D,varargin)


%% set input parameters
arg.dgrid = linspace(20,150,24)';
arg.dmx = 800;
arg.woundedge = 225; 
arg.k = 8; 
arg = parseVarargin(varargin,arg); 

dGrid = arg.dgrid; 
Dmx = arg.dmx; 
WoundEdge=arg.woundedge; 
MIvec = nan(numel(dGrid),1); 

for iG=1:numel(dGrid)

    %% group into bins
    grp = nan(size(D));
    Dgrid = WoundEdge:dGrid(iG):Dmx;
    for i=2:numel(Dgrid)
        grp(Dgrid(i-1)<D & D<=Dgrid(i))=i-1;
    end
    
    %% create cell array for Ca by bin
    CaGrp = grp2cell(Mx(~isnan(grp))',grp(~isnan(grp)));
    CaGrp = cellfun(@(x) x(:,1:end)',CaGrp,'uniformoutput',0);
    
    %% distance distribution
    Qs = cellfun(@numel,CaGrp); 
    Qs = Qs./sum(Qs); 
    
    %% calcualte MI
    try
%         PD = cell(1,numel(CaGrp)); 
%         for igrp = 1:numel(CaGrp)
%            PD{igrp} = ProbDistUnivKernel(CaGrp{igrp}(:)); 
%         end
        MIvec(iG) = calcChannelCapacity(CaGrp,'qopt',Qs,'jackknife',false,'adaptivek',true); %,'festim',PD
    catch e
        warning(e.message)
        MIvec(iG)=nan; 
    end
    
end

%% regress to zero bin size
b = regress(MIvec(:),[ones(size(MIvec(:))) dGrid(:)]); 
MI=b(1); 
