function ErkSB = subSampleErkData(Erk,Terk,varargin)

if iscell(Erk)
    numImages = size(Erk{1},1);
    seconds = (numImages-1).*60;
    if numImages<60
        numImages = 60;
        seconds = 3540;
    end
else
    numImages = 60;
    seconds = 60*59;
end

arg.flt = sum(fspecial('gauss',7,2));
arg.tgrid = (0:numImages:seconds)';
arg.dim = 10; 
arg.tsample = []; % which minutes to sample from
arg = parseVarargin(varargin,arg); 

% Erk could be a matrix or a cell array
if iscell(Erk) && iscell(Terk)
    ErkSB = cellfun(@(e,t) subSampleErkData(e,t,arg),Erk,Terk,'uniformoutput',0);
    return
end

Tgrid = arg.tgrid; 
Dim = arg.dim; 
flt = arg.flt; 

Erk(:,max(isnan(Erk)))=[]; 
Erk = interp1(Terk(:),Erk,Tgrid); 
Erk = filtfilt(flt,1,Erk); 
if isempty(arg.tsample)
    ErkSB = interp1((0:59)',Erk,linspace(0,59,Dim+2));
    ErkSB = ErkSB(2:end-1,:);
else
    ErkSB = interp1((0:60)',Erk,arg.tsample);
end
    
    
