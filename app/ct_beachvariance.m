function [res,meta] = ct_beachvariance(obj,mobj,caserec)
%
%-------function help------------------------------------------------------
% NAME
%   ct_beachvariance.m
% PURPOSE
%   Compute the beach variability over a set of profiles
% USAGE
%   [res,meta] = ct_beachvariance(obj,mobj,caserec);
% INPUTS
%   obj - handle to CT_BeachAnalysis class object
%   mobj - handle to CoastalTools UI
%   caserec - ids pf selected profiles
% OUTPUT
%   res - tabular summary of results
%   meta - metadata for time and profiles used
% NOTES
%   uses functions in CT_BeachAnalysis.m and 
% SEE ALSO
%   flist = matlab.codetools.requiredFilesAndProducts('ct_beachvariance.m')'
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2026
%--------------------------------------------------------------------------
%
    if isempty(mobj.Inputs) || ~isfield(mobj.Inputs,'ctWaveParameters')
        warndlg('Site parameters not defined');
        return; 
    end
    site = mobj.Inputs.ctWaveParameters;

    muicat = mobj.Cases.Catalogue;
    proflist = muicat.CaseDescription(caserec);         %selected profles
    np = numel(proflist);

    inp = inputdlg({'Upper limit (mOD)','Lower limit(mOD)','Interval (m)',...
                     'Save results'},'Beach levels',1,{'6','0','0.5','0'});
                                    
    if isempty(inp), res = []; meta = []; return; end   %user cancelled

    zmax = str2double(inp{1});
    zmin = str2double(inp{2});
    zint = str2double(inp{3});
    zlevels = zmin:zint:zmax;
    nz = numel(zlevels);

    %get the composite time intervals for all profiles
    [ptime,~] = ProfileTimes(mobj,caserec);
    nt = numel(ptime);
    nv = 5;                            %number of variables
    res = cell(1,nv);                   %5 variables each cell [nt,np,nz]
    for iv = 1:nv,   res{iv} = zeros(nt,np,nz); end 

    hw = waitbar(0,'Processing');
    for k=1:nz
        %function returns pos={E,N,xD,BS,SLang} and struct of regression coefficients
        [temp,~] = getPositionAndRates(obj,mobj,caserec,ptime,...
                                                        zlevels(k),site);
        for iv = 1:nv
            res{iv}(:,:,k) = temp{iv};
        end
        waitbar(k/nz,hw)
    end

    meta.zdata = zlevels;
    meta.time = ptime;  meta.proflist = proflist;
    meta.issave = logical(str2double(inp{4}));    %save results flag
    close(hw)
end
%%
function [time,dsts] = ProfileTimes(mobj,caserec)
    %get the composite time intervals for all profiles and an array of the
    %datasets (dsts)
    npro = length(caserec);            
    time = [];           %time intervals in any of the surveys
    dsts(npro,1) = dstable;
    for j=1:npro
        dsts(j) = getDataset(mobj.Cases,caserec(j),1); %idset=1
        newtime = dsts(j).RowNames;
        time = unique(vertcat(time,newtime));
    end
end

