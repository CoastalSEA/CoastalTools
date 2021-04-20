classdef CT_BeachAnalysis < muiDataSet                  
%
%-------class help------------------------------------------------------
% NAME
%   CT_BeachAnalysis.m
% PURPOSE
%   Class for Model CT_BeachAnalysis to be run from CoastalTools
% SEE ALSO
%   muiDataSet
%   For full list of functions used:
%   flist = matlab.codetools.requiredFilesAndProducts('CT_BeachAnalysis.m')'
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties: 
        VolumeCalcs      %struct for results from volume calcs
        ShorelineCalcs   %struct for results from shoreline position calcs
        FitCoeffs   %used to store BMV model coefficients with volume data
                    %used to store regression results for shoreline position
        BVI         %profile indices from BVI analysis
    end
    
    properties (Hidden)
        leglimit = 12    %limit to number of rows in legend 
        ModelType        %model used for the particular instance
    end
    
    properties (Transient)
        ModelList = {'Volumes','Shoreline position','Shoreline',...
                     'BVI profile set','Shore profile','Dean profile'}          
    end
    
    methods (Access = private)
        function obj = CT_BeachAnalysis()                    
            %class constructor
        end
    end      
%%
    methods (Static)                 
        function runModel(mobj,src)
            %function to run a simple 2D diffusion model
            obj = CT_BeachAnalysis;    
            id_model = find(strcmp(obj.ModelList,src.Text));  
            if isempty(id_model), return; end
            dsp = modelDSproperties(obj,id_model);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Load beach profiles to use these tools');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); 
            %--------------------------------------------------------------
            % Model code
            %--------------------------------------------------------------
            % output - struct with summary of model results
            %   results - cell array of output to be saved
            %   proflist - profile id or x-y location data
            %   modeltime - datetime for model output
            %   metatxt - additional text to define conditions of run
            %   source - text to define input data sources
            switch id_model
                case 1   %Volumes for individual profile
                    output = volumesModel(obj,mobj);
                    mtype = 'Volumes_model';
                    if ~isempty(output)
                        obj.VolumeCalcs = output;
                        src = findobj(mobj.mUI.Tabs.Children,'Tag','Volumes'); 
                        headtxt = sprintf('Profile: %s\nFor: %s',...
                                      output.source,output.metatxt);
                        tablefigure(src,headtxt,output.table);
                    end
                case 2  %Shoreline position for individual profile
                    output = shorelinePositionModel(obj,mobj);
                    mtype = 'ShorePos_model';
                    if ~isempty(output)
                        obj.ShorelineCalcs = output.table;
                        src = findobj(mobj.mUI.Tabs.Children,'Tag','Shoreline'); 
                        headtxt = sprintf('Profile: %s\nFor: %s',...
                                      output.source,output.metatxt);
                        tablefigure(src,headtxt,output.table);
                    end
                case 3  %Shoreline (position extracted from set of profiles)
                    output = profiles2shoreline(obj,mobj);
                    mtype = 'SLine_model';
                    if ~isempty(output)
                        obj.FitCoeffs = output.table;                
                        src = findobj(mobj.mUI.Tabs.Children,'Tag','BVI');
                        headtxt = 'Shoreline fit coefficients';
                        tablefigure(src,headtxt,output.table);
                    end
                case 4 %beach vulnerability index
                    output = getBVindex(obj,mobj);
                    mtype = 'BVI_model';
                    if ~isempty(output)
                        obj.FitCoeffs = output.table;                
                        obj.BVI = output.bvitable;
                        src = findobj(mobj.mUI.Tabs.Children,'Tag','BVI'); 
                        headtxt = 'Beach vulnerability index table';
                        tablefigure(src,headtxt,output.bvitable);
                    end
                case 5  %plot an idealised profile and closure depths
                    src = findobj(mobj.mUI.Tabs.Children,'Tag','Profile'); 
                    bmvProfile(obj,mobj,src);
                    output = [];
                case 6  %plot the Dean profile
                    src = findobj(mobj.mUI.Tabs.Children,'Tag','Profile'); 
                    getDeanProfile(obj,mobj,src);
                    output = [];
            end
            if isempty(output) || isempty(output.results{1})
                return; %user cancelled or no results returned
            end       
            %--------------------------------------------------------------
            % Assign model output to a dstable using the defined 
            % dsproperties meta-data
            %--------------------------------------------------------------                  
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(output.results{:},'RowNames',output.modeltime,...
                                                    'DSproperties',dsp);
            %--------------------------------------------------------------
            % Save results
            %--------------------------------------------------------------                        
            %assign metadata about model
            obj.ModelType = obj.ModelList{id_model};
            dst.Source = sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                       obj.ModelType);
            dst.MetaData = output.metatxt;
            dst.UserData = output.proflist;
            %save results
            casename = {sprintf('%s-%s',obj.ModelType,output.source)};    
            setDataSetRecord(obj,muicat,dst,mtype,casename,true);
            getdialog('Run complete');
        end
    end
 %%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end   
%%
        function tabCalcs(obj,src)
            %generate table to display calcs output
            ht = src.Children;   %clear any existing tab content
            delete(ht);
            %prompt selection of compatible dataset for chosen tab
            idx = strcmp({obj(:).ModelType},src.Tag);            
            caselist = arrayfun(@(x) x.Data.Dataset.Description,obj,'UniformOutput',false);
            caselist(~idx) = [];   
            if isempty(caselist)
                warndlg(sprintf('Model ouput of type: %s not found',src.Tag));
                return;
            end

            if length(caselist)>1
                [classrec,ok] = listdlg('Name','CalcTab','PromptString','Select Case', ...                                 
                                     'ListSize',[200,80],'SelectionMode','single', ...                                 
                                     'ListString',caselist);
                if ok==0, return; end 
            else
                classrec = 1;
            end
            %display relevant table for selected tab
            switch src.Tag
                case  'Volumes'
                    output = obj(classrec).VolumeCalcs;
                    headtxt = sprintf('Profile: %s\nFor: %s',...
                                      output.source,output.metatxt);
                    tablefigure(src,headtxt,output.table);
                case 'Shoreline'
                    output = obj(classrec).ShorelineCalcs;
                    headtxt = sprintf('Profile: %s\nFor: %s',...
                                      output.source,output.metatxt);
                    tablefigure(src,headtxt,output.table);                    
                case 'BVI'
                    if strcmp(obj(classrec).ModelType,'Shoreline')
                        output = obj(classrec).FitCoeffs;
                        headtxt = 'Shoreline regression fit coefficients';
                        output_table = output.table;  
                    else
                        output = obj(classrec).BVI;
                        headtxt = 'Beach vulnerability index table';
                        output_table = output.bvitable;                    
                    end
                    tablefigure(src,headtxt,output_table);
            end
        end
%% ------------------------------------------------------------------------
% UI Functions called directly from main menu
%--------------------------------------------------------------------------
        function getShorelinePlot(obj,mobj)
            %plot the results generated by profiles2shoreline and stored 
            %as a set of shorelines over time
            [tsc_time,Xdata,Ydata,Loc,pid,vname,caserec] = selectShoreline(obj,mobj);
            if isempty(tsc_time), return; end   
            obj = getCase(mobj.Cases,caserec);
            time = cellstr(tsc_time);
            hf = figure('Name','Shorelines','Tag','PlotFig','Units','normalized');
            %move figure to top right
            hf.Position(1) = 1-hf.Position(3)-0.01;
            hf.Position(2) = 1-hf.Position(4)-0.12;
            nline = length(time);
            set(groot,'defaultAxesColorOrder',jet(nline))
            ax = axes;
            if isempty(Xdata)      %user requested Chainage plot
                lobj = mobj.Cases.DataSets.ctBeachProfileData;  %handle to beach profile data
                nint = subsampleprofiles(lobj,length(Loc));
                [Loc,Ydata] = checkDirection(obj,Loc,Ydata);
                plotShoreline_CH(obj,ax,time,Loc,Ydata,pid,vname,nint);
            else                   %user requested E&N plot
                plotProfileBaseline(obj,mobj,ax,Loc);
                plotShoreline_EN(obj,ax,time,Xdata,Ydata,pid);                
            end
        end
%%
        function getShoreTablePlot(~,mobj,src)
            %plot the indices for a BVI set analysis
            
            switch src
                case 'Rates plot'
                    type = {'SLine_model','BVI_model'};
                case 'BVI set plot'
                    type = {'BVI_model'};
            end
            muicat = mobj.Cases;
            promptxt = 'Select shoreline';

            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                        'CaseType',type,...       
                        'SelectionMode','single','ListSize',[250,200]);
            if ok<1, return; end
            obj = getCase(mobj.Cases,caserec);
            dst = obj.Data.Dataset;
            pid = dst.Description;
            Loc = dst.UserData;          
            
            hf = figure('Name','Shore change','Tag','PlotFig','Units','normalized');
            %move figure to top right
            hf.Position(1) = 1-hf.Position(3)-0.01;
            hf.Position(2) = 1-hf.Position(4)-0.12;
            ax = axes;
            lobj = mobj.Cases.DataSets.ctBeachProfileData;  %handle to beach profile data
            nint = subsampleprofiles(lobj,length(Loc));            
            switch src
                case 'Rates plot'
                    outtable = obj.FitCoeffs;
                    [Loc,outtable] = checkDirection(obj,Loc,outtable);
                    plotShoreRates_CH(obj,ax,Loc,outtable,pid,nint);
                case 'BVI set plot'
                    outtable = obj.BVI;
                    [Loc,outtable] = checkDirection(obj,Loc,outtable);
                    plotBVI_CH(obj,ax,Loc,outtable,pid,nint);
            end           
        end
%%
        function getProfileCentroidsPlot(~,mobj)
            %setup call to function to plot centroids
             muicat = mobj.Cases;
            promptxt = 'Select profile volumes data set';
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                        'CaseType',{'Volumes_model'},...                    
                        'SelectionMode','single','ListSize',[250,200]);
            if ok<1, return; end
            [dst,~] = getDataset(muicat,caserec,1);   %idset=1 
            %reformat time to short form for use as labels
            dst_time = cellstr(dst.RowNames,'dd-MM-yyyy');
            x = dst.x1;        %N_D x-centroid
            labels.xlabel = dst.VariableLabels{3};
            z = dst.z1;        %N-D z-centroid
            labels.ylabel = dst.VariableLabels{4};
            labels.title = sprintf('Profile No: %s ',dst.Description);
            labels.tformat = 'dd-mm-yyyy';
            labels.pname = 'Centroids';
            phaseplot(x,z,dst_time,labels);
            
            complex_vectors(x,z);
        end
%%
        function getProfileSpaceTimePlot(obj,mobj)
            %setup selection of profiles for plotting as a space-time plot
            muicat = mobj.Cases;
            if length(obj)>1
                type = {obj(:).ModelType};
                utype = unique(type);
                answer = questdlg('Which model','Space-time plot',utype{:},utype{1});
                idobj = strcmp(type,answer);
                obj = obj(idobj);                
            end
            caserec = getMinShorelineSet(obj,mobj,2);  %2=min number of profiles to select
            if isempty(caserec), return; end  %user cancelled
            profnames = muicat.Catalogue.CaseDescription(caserec);            
            j = 1;
            catdesc = muicat.Catalogue.CaseDescription;
            for i=1:length(obj)
                adst = obj(i).Data.Dataset;
                matchcell = @(x) contains(adst.Description,x);
                if any(cellfun(matchcell,profnames))
                    dsts(j) = adst;
                    selrecs(j) = find(strcmp(catdesc,adst.Description));                    
                    j = j+1;
                end
            end
            time = getProfileTimes(obj,mobj,caserec);
            
            %if there is more than one variable prompt user to select one
            varnames = dsts(1).VariableNames;            
            idvar = 1;
            if length(varnames)>1
                promptxt = 'Select the variable to plot:';
                [idvar,ok] = listdlg('Name','Space-time plot','SelectionMode','single',...
                                'PromptString',promptxt,'ListString',varnames);
                if ok<1, return; end
            end
            vardesc = adst.VariableDescriptions{idvar};
            %unpack data
            nrec = length(dsts);
            datavar = NaN(length(time),nrec);
            for i=1:nrec
                ptime = dsts(i).RowNames;  %time intervals for a profile
                idp = ismember(time,ptime);          %ids of ptime in time    
                datavar(idp,i) = dsts(i).(varnames{idvar});
            end
            
            plotSpaceTime(obj,time,profnames,datavar,vardesc)
        end   
    end
%%
%--------------------------------------------------------------------------
% Main functions for CT_BeachAnalysis
%--------------------------------------------------------------------------    
    methods (Access=private)
        function output = volumesModel(obj,mobj)
            %compute the volumes per metre shore (CSA) for a timeseries
            %of beach profiles.
            % INPUTS
            % x,z,t - beach profile data
            % xmin  - lower limit for volume calculation
            % zmin  - lower limit for volume calculation
            % FUNCTION CALLS (external) - section_centroid.m
            output = [];
            [dst_time,x,z,pid] = selectBeachProfile(obj,mobj);  
            if isempty(dst_time), return; end
            figtitle = sprintf('Volume Model for profile: %s',pid);
            promptxt = 'Accept Min/Max definition';
            tag = 'VolsFig';
            butnames = {'Yes','No'};
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames);
            h_ax = axes(h_plt);
            plotBeachProfiles(obj,h_ax,x,z);
            mxmn = getMaxMin(obj,x,z);
            %get input parameters from user
            mxmn = getSampleBox(obj,mxmn,h_but,h_ax);  
            %get profiles within the defined box
            [X,Z,xrange,zrange] = getSampledDataPoints(obj,x,z,mxmn);
            %get the dimensional volumes and centroids (isND=false)
            [V,X1,Z1] = getVolumeCentroids(obj,X,Z,xrange,zrange,false);
            %get the non-dimensional volumen and centroids (isND=true)
            [m0,x1,z1] = getVolumeCentroids(obj,X,Z,xrange,zrange,true);            
            
            %assign ouput
            slope = xrange./zrange;
            output.results = {V,m0,x1,z1,slope};
            output.modeltime = dst_time;
            output.metatxt = sprintf('x from %g to %g; z from %g to %g',...
                mxmn.xmin,mxmn.xmax,mxmn.zmin,mxmn.zmax);  
            output.source = pid;
            output.proflist = pid;
            
            %data to be written to Calcs Tab
            output.table = table(V,xrange,zrange,m0,X1,Z1,...
              'VariableNames',{'Volume','Xrange','Zrange','m0','X1','Z1'});      
            output.table.Properties.RowNames = cellstr(dst_time); 
            output.table.Properties.Description = ...
                sprintf('Volume and centroid data for %s',pid);
        end
%%
        function output = shorelinePositionModel(obj,mobj)
            %compute the shoreline position and slope for a timeseries
            %of beach profiles.
            % INPUTS
            % x,z,t - beach profile data
            % zlevel  - level used to calculate position and slope
            % delz - +/- z value to calculate slope over
            % minx - starting x to avoid rear slope low points
            % FUNCTION CALLS (external)
            % none
            output = [];
            [dst_time,x,z,pid] = selectBeachProfile(obj,mobj);
            if isempty(dst_time), return; end
            figtitle = sprintf('Shoreline position model for profile: %s',pid);
            promptxt = 'Accept Min/Max definition';
            tag = 'VolsFig';
            butnames = {'Yes','No'};
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames);
            h_ax = axes(h_plt);
            plotBeachProfiles(obj,h_ax,x,z);
            
            %get input parameters from user
            [xdist,slope,opt] = getShorelinePositions(obj,x,z,...
                                                    dst_time,h_but,h_ax);
            %assign ouput
            output.results = {xdist,abs(slope)};
            output.modeltime = dst_time;
            output.metatxt = sprintf('zlevel = %g; z +/- %g; xmin %g',...
                                           opt.zlevel,opt.delz,opt.xmin);             
            output.source = pid;
            output.proflist = pid;
            
            %write data to Calcs Tab
            output.table = table(xdist,abs(slope),...
                                    'VariableNames',{'Chainage','Slope'}); 
            output.table.Properties.RowNames = cellstr(dst_time); 
            output.table.Properties.Description = ...
                sprintf('Shoreline position and slope for %s at %g mOD',...
                                                          pid,opt.zlevel);
        end
%%
        function output = profiles2shoreline(obj,mobj)
            %sample the available profiles and extract the shoreline at some
            %elevation, return vectors of position and rates of change
            output = [];
            %prompt user to select profiles to be used (multiple profiles)
            caserec = getMinShorelineSet(obj,mobj,2);    %2=min number of profiles to select 
            if isempty(caserec), return; end
            
            prompt = {'Elevation that defines shoreline:'};
            zlevel = getShoreLevel(obj,prompt);
            if isempty(zlevel), return; end %user cancelled
            
            %get site data
            msgtxt = 'No site parameters defined';
            site =  getClassObj(mobj,'Inputs','ctWaveParameters',msgtxt);   
            if isempty(site), return; end
            
            time = getProfileTimes(obj,mobj,caserec);
            [pos,stats,meta] = getPositionAndRates(obj,mobj,caserec,time,zlevel,site);
            if isempty(pos), return; end
            %setup output for writing to stored dataset
            output.results = pos;
            output.modeltime = time;
            output.metatxt = sprintf('zlevel = %g; z +/- %g; xmin %g',...
                                                    zlevel,0.5,meta.xmin);             
            output.source = sprintf('using %g profiles at %.2f mOD',meta.npro,zlevel);
            %list of profiles used to construct shoreline
            output.proflist = meta.proflist;  

            %write data to Calcs Tab
            output.table = stats; 
            output.table.Properties.RowNames = meta.proflist; 
            output.table.Properties.Description = ...
                sprintf('Change in shoreline position, slope and angle at %g mOD',...
                                                          zlevel);   
            output.table.Properties.VariableUnits = {'m','m/yr','-','m',...
                                                'm','m/yr','-','m',...
                                                'deg','deg/yr','-','deg'};
        end
%%
        function output = getBVindex(obj,mobj)
            %set up call to beachvulnerability and handle output
            %setup output for writing to stored dataset
            output = [];
            %prompt user to select profiles to be used (multiple profiles)
            caserec = getMinShorelineSet(obj,mobj,2);  %2=min number of profiles to select 
            if isempty(caserec), return; end
            npro = length(caserec);
            [mnpos,mnstats,itable,meta] = ct_beachvulnerability(obj,mobj,caserec);
            if isempty(itable), return; end
            
            %setup output for writing to stored dataset
            output.results = mnpos;
            output.modeltime = meta.time;
            output.metatxt = {sprintf('zlevel = %g; z +/- %g',0,0.5)};                                                                 
            output.source = sprintf('using %d profiles at %.1f mOD',npro,0);
            %list of profiles used to construct shoreline
            output.proflist = meta.proflist;
            
            %write data FitCoeffs property (not shown on Calcs tab)
            output.table = mnstats; 
            output.table.Properties.RowNames = meta.proflist; 
            output.table.Properties.Description = ...
                sprintf('Change in shoreline position, slope and angle at %d mOD',...
                                                          0);   
            output.table.Properties.VariableUnits = {'m','m/yr','-','m',...
                                                'm','m/yr','-','m',...
                                                'deg','deg/yr','-','deg'};
            %write data to Calcs Tab
            output.bvitable = itable;
            output.bvitable.Properties.Description = ...
                sprintf('Profile indices for %d profiles with "exposure” period of N = %g years',npro,meta.Nyear);                       
        end
%%
        function bmvProfile(~,mobj,src)
            %plot the BMV profile with surf and shoal depths
            %get handle to model input parameters
            h_inp = getClassHandle(mobj,'Sim_BMVinput');
            bmvinp = mobj.(h_inp);
            if isempty(bmvinp)
                warndlg('BMV input parameters have not been defined');
                return;
            end
            %class instance for inshore wave data
            inwaveobj = mobj.(getClassHandle(mobj,'InWaveModel'));
            %retrieve Hsi from inshore wave data set
            [hsts,recnum] = getInWaveTS(inwaveobj,mobj,'Hsi');
            tpts = getWavePeriod(inwaveobj,mobj,recnum);
            
            idx = [15,16];     %subset of variables to be checked
            bmv = setBMVmodelParams(bmvinp,mobj,idx);
            %mean(Hs),std(Hs),mean(Tp)used in option 3 of ClosureOption
            %for Hallermeier zones
            bmv.meanWaves.mH = mean(hsts,'MissingData','remove');  
            bmv.meanWaves.sH = std(hsts,'MissingData','remove');  
            bmv.meanWaves.mT = mean(tpts,'MissingData','remove');
            
            prompt = 'Specify Hs & Tp or a percentile (e.g. 95)';
            answer = inputdlg(prompt,'Wave definition',1);
            if isempty(answer), return; end  %user aborted
            
            if strcmp(answer{1},'')
                %no value defined so use mean values
                Hs = bmv.meanWaves.mH; Tp = bmv.meanWaves.mT;
            else   
                inputval = str2num(answer{1}); %#ok<ST2NM>  %input can be vector
                if length(inputval)==2
                    %use user defined Hs and Tp
                    Hs = inputval(1); Tp = inputval(2);
                else
                    %use percentile wave ehight and indicative  wave period
                    Hs = prctile(hsts.Data,inputval);
                    %limiting steepness based of 1/18 gives coeff=3.4
                    %for spanish coast BMV, 2003, suggest coeff~=6
                    Tp = 6*sqrt(Hs);
                end
            end
                
            [Xr,Zr] = referenceProfile(bmv,Hs,Tp,[]);
            
            %get limit points on profile
            Kr = bmvFittingCoeffs(bmv,Hs,Tp);
            [ha_r,hr_r] = getProfileDepths(bmv,Hs,Tp,Kr);
            XZa(2) = bmv.Zlw-ha_r;
            XZa(1) = interp1(Zr,Xr,XZa(2)); 
            XZr(2) = bmv.Zlw-hr_r;
            XZr(1) = interp1(Zr,Xr,XZr(2)); 
            
            %display on plot
            mnwv = [bmv.meanWaves.mH,bmv.meanWaves.sH,bmv.meanWaves.mT];
            titletxt = sprintf('Shore profile for Hs=%.1f, Tp=%.1f',Hs,Tp);
            tabShoreProfile(bmv,src,Xr,Zr,titletxt,XZa,XZr,mnwv)           
        end
%%
        function getDeanProfile(~,mobj,src)
            %plot the Dean's profile using values specified for wave model
            %get handle to model input parameters
            inp = mobj.Inputs.ctWaveParameters;
            if isempty(inp)
                warndlg('Input parameters have not been defined');
                return;
            end
            zBC = inp.BeachCrestLevel;
            z1km = inp.BedLevelat1km;
            ubs = inp.UpperBeachSlope;
            
            if length(z1km)>1
                y1km = z1km(1);       %y and z have been specified
                z1km = z1km(2);
            else
                y1km = 1000;          %default distance
            end
            
            yint = 1;
            y = 1:yint:y1km;
            [yp,zp,~] = deanbeachprofile(y,zBC,[y1km,z1km],ubs,false); 
        
            ha = findobj(src.Children,'-depth',1);
            delete(ha)
            ax = axes(src,'Tag','ShoreProfile');
            
            plot(ax,yp,zp,'-k');
            hold on
            y1km = y1km+zBC*ubs;
            plot(ax,y1km,z1km,'og');
            txtstr = sprintf('z1km=%.1f (mOD)   ',z1km);
            text(ax,y1km,z1km,txtstr,'FontSize',6,'HorizontalAlignment','right')
            plot(ax,0,zBC,'og');
            txtstr = sprintf('   zBC=%.1f (mOD)',zBC);
            text(ax,0,zBC,txtstr,'FontSize',6)
            plot(ax,xlim, 0*[1 1],'--b');
            plot(ax,zBC*ubs,0,'or');
            xlabel('Distance (m)')
            ylabel('Elevation (mOD)')
            title(ax,sprintf('Dean profile from 0 to %.1f (mOD)',z1km));            
            hold off
        end 
            
%% ------------------------------------------------------------------------
% UI Functions called by beach models
%--------------------------------------------------------------------------
        function [dst_time,x,z,pid] = selectBeachProfile(~,mobj)
            %user selects a single beach profile - returns chainage and
            %elevation time series for profile
            muicat = mobj.Cases;
            promptxt = 'Select beach profile';
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                                'CaseClass',{'ctBeachProfileData'},...
                                'SelectionMode','single','ListSize',[250,200]);
            if ok<1
                warndlg('No selection made');
                dst_time = []; x = []; z = []; pid = [];
                return;
            end
            [dst,~] = getDataset(muicat,caserec,1);   %idset=1 
            dst_time = dst.RowNames;
            x = dst.Chainage;
            z = dst.Elevation;
            pid = dst.Description;
        end
%%
        function [dst_time,E,N,Loc,pid,vname,caserec] = selectShoreline(obj,mobj)
            %user selects a shoreline position data set - returns Eastings and
            %Northings time series
            vname = [];
            muicat = mobj.Cases;
            promptxt = 'Select shoreline';

            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                                'CaseType',{'SLine_model'},...    
                                'SelectionMode','single','ListSize',[250,200]);
            if ok<1
                warndlg('No selection made');
                dst_time = []; E = []; N = []; Loc = []; pid = [];
                return;
            end
            [dst,~] = getDataset(muicat,caserec,1);   %idset=1 
            [dst_time,E,N,Loc,pid,vname] = selectShorelineType(obj,dst);
        end
%%
        function [dst_time,E,N,Loc,pid,vname] = selectShorelineType(~,dst)
            %select whether to plot E&N or chainage and return data
            vname = [];
            dst_time = dst.RowNames;
            pid = dst.Description;
            Loc = dst.UserData;  %list of profile ids;
            type = questdlg('Plot as E&N or Profile location?','Shoreline',...
                              'E&N','Profile location','E&N');
            if strcmp(type,'E&N')            
                E = dst.Eastings;
                N = dst.Northings;            
            else
                E = [];
                vname = questdlg('Select variable to plot?','Shoreline',...
                               'Chainage','Slope','Orientation','Chainage');
                switch vname
                    case 'Chainage'
                        N = dst.Chainage; 
                    case 'Slope'
                        N = dst.Slope; 
                    case 'Orientation'
                        N = dst.SLang; 
                end
            end
        end
%%
        function caserec = getShorelineSet(obj,mobj)
            %prompt user to select a set of profiles 
%             h_prof = getClassHandle(mobj,'BeachProfileData');
            muicat = mobj.Cases.Catalogue;
            idxlist = strcmp(muicat.CaseClass,'ctBeachProfileData');
            promptxt = 'Select profiles to be used:';
            caselist = muicat.CaseDescription(idxlist)';
            classid = muicat.CaseID(idxlist);
            %convert class selection to caserec values
            nselect = length(classid);
            caserec = zeros(1,nselect);
            for i=1:nselect
                caserec(i) = find(muicat.CaseID==classid(i));
            end
            %check that profiles are in order along coast
            [~,~,idd] = shorelineProfileOrder(obj,mobj,caserec);
            caserec = caserec(idd);
            caselist = caselist(idd);
            [idp,ok] = listdlg('Name','Shoreline extraction', ...
                            'PromptString',promptxt,'ListString',caselist);                        
            if ok<1, caserec = []; return; end
            caserec = caserec(idp);
        end
%%
        function caserec = getMinShorelineSet(obj,mobj,nprof)
            %get user to select a minimum number of profiles defined by nprof
            caserec = [];
            ok = 0;
            while ok<1
                caserec = getShorelineSet(obj,mobj);   
                if isempty(caserec)
                    return;
                elseif length(caserec)>=nprof
                    ok = 1;
                else
                    hw = warndlg(sprintf('Select a minimum of %d profiles',nprof));
                    uiwait(hw);
                end
            end
        end       
%%
        function zlevel = getShoreLevel(~,prompt)
            %prompt user to define the elevation of the shoreline
            title = 'Profile shoreline';
            numlines = 1;
            default = {num2str(0)};    
            answer = inputdlg(prompt,title,numlines,default);
            if isempty(answer), zlevel=[]; return; end %user cancelled
            zlevel = str2double(answer);
        end
%%
        function [location,ydata] = checkDirection(~,location,ydata)
            %prompt user to use default direction of profiles or reverse            
            promptxt = sprintf('Profiles are from %s to %s, reverse?',...
                                                location{1},location{end});
            answer = questdlg(promptxt,'Shore direction','Accept','Reverse','Accept');            
            if strcmp(answer,'Reverse')
                location = location(end:-1:1);
                if istable(ydata)
                    ydata{:,:} = flipud(ydata{:,:});
                else                    
                    ydata = fliplr(ydata);  
                end
            end
        end
%%
        function time = getProfileTimes(~,mobj,caserec)
            %get the composite time intervals for all profiles
            npro = length(caserec);            
            time = [];           %time intervals in any of the surveys
            for j=1:npro
               [dst,~] = getDataset(mobj.Cases,caserec(j),1);  %idset=1 
               newtime = dst.RowNames;
               time = unique(vertcat(time,newtime));
            end
        end
%%
        function mxmn = getSampleBox(~,mxmn,h_but,h_ax)
            %get min and max dimensions of sample box to be used
            %mxmn  - initial estimates of min/max values
            %h_but - handle for button in acceptfigure
            %h_ax  - handle to plot axes used for profile plots
            prompt = {'Minimum X distance (m):','Minimum Z elevation (mOD):'};
            title = 'Define input parameters';
            numlines = 1;
            defaultvalues{1} = num2str(mxmn.xmin);
            defaultvalues{2} = num2str(mxmn.zmin);
            h1 = [];
            ok = 0;
            while ok<1
                useInp=inputdlg(prompt,title,numlines,defaultvalues);
                if isempty(useInp), return; end %user cancelled
                xmin = str2double(useInp{1}); %useInp.NumVals(1);
                xmax = mxmn.xmax;
                zmin = str2double(useInp{2}); %useInp.NumVals(2);
                zmax = mxmn.zmax;
                delete(h1)
                hold on
                h1 = plot(h_ax,[xmin,xmin,xmax,xmax,xmin],...
                               [zmin,zmax,zmax,zmin,zmin],'--r');
                hold off
                waitfor(h_but,'Tag');
                if ~ishandle(h_but) %this handles the user deleting figure window    
                    ok = 0;
                elseif strcmp(h_but.Tag,'Yes')
                    ok = 1; %box definition accepted                    
                else
                    %Prompt user to change definition
                    defaultvalues{1} = num2str(xmin);
                    defaultvalues{2} = num2str(zmin);
                    h_but.Tag = '';                     
                end 
            end
            mxmn.xmin = xmin; mxmn.xmax = xmax;
            mxmn.zmin = zmin; mxmn.zmax = zmax;
        end 
%%
        function mxmn = getMaxMin(~,x,z)
            %get the maximum and minimum values for the x and z data sets
            %mxmn - struct containing the 4 values
            mxmn.xmin = max(min(x,[],2));  %largest value of min of any profile
            mxmn.xmax = max(x(:));         %largest value across all profiles
            mxmn.zmin = max(min(z,[],2));
            mxmn.zmax = max(z(:));
        end        
%%
%--------------------------------------------------------------------------
% Geometry Functions called by beach models
%--------------------------------------------------------------------------
        function [pos,stats,meta] = getPositionAndRates(obj,mobj,...
                                                    caserec,time,zlevel,site)
            %for each time interval of selected profiles compute position, 
            %slope, orientation and associated rates of change                                            
            npro = length(caserec);   %number of profiles
            %check that profiles are in order along coast
            [~,~,idd] = shorelineProfileOrder(obj,mobj,caserec);
            caserec = caserec(idd);
            %now analyse each time interval
            if isscalar(zlevel), zlevel = repmat(zlevel,1,npro); end
            ntime = length(time);     %composite time for all profiles
            xD = NaN(ntime,npro); E = NaN(ntime,npro); N = E; BS = xD; %theta = xD;
            bE = E; bN = N;
            aX = zeros(npro,1); bX = aX; RsqX = aX; seX = aX;
            aM = aX; bM = aX; RsqM = aX; seM = aX;
%             aO = aX; bO = aX; RsqO = aX; seO = aX;
            proflist{npro} = [];
            for k=1:npro
                [dst,~] = getDataset(mobj.Cases,caserec(k),1);  %idset=1
%                 dst = getCaseDataSet(mobj.Cases,mobj,caserec(k));
                x = dst.Chainage;
                z = dst.Elevation;
                xmin = max(min(x,[],2));
                %find the distance to and slope at zlevel
                [xdist,slope] = slope_points(x,z,zlevel(k),0.5,xmin);    
                %ptime is in profile and time is unique set for all profiles   
                ptime = dst.RowNames;  %time intervals for a profile
                idp = find(ismember(time,ptime));   %ids of ptime in time           
                xD(idp,k) = xdist;
                BS(idp,k) = slope;               
                proflist{k} = mobj.Cases.Catalogue.CaseDescription{caserec(k)};
                %if the dataset includes E,N get the E,N of zlevel
                if any(strcmp('Eastings',dst.VariableNames))
                    [Es,Ns,Eb,Nb] = ENofShoreline(obj,dst,xdist);   
                    E(idp,k) = Es;
                    N(idp,k) = Ns; 
                    bE(idp,k) = Eb;
                    bN(idp,k) = Nb;
                end
%                 shoreang = shore_orientation(Es,Ns,Eb,Nb);
%                 theta(idp,k) = shoreang;
                %now get rate of change in distance and slope
                %eps(0) to avoid divide by zero in linear regression
                mtime = years(ptime-ptime(1)+eps(0)); 
                [aX(k),bX(k),RsqX(k)] = regression_model(mtime,xdist,'Linear');
                seX(k) = stderror(mtime,xdist,aX(k),bX(k),'Linear');
                [aM(k),bM(k),RsqM(k)] = regression_model(mtime,-slope,'Linear');
                seM(k) = stderror(mtime,-slope,aM(k),bM(k),'Linear');
%                 [aO(k),bO(k),RsqO(k)] = regression_model(mtime,shoreang,'Linear');
%                 seO(k) = stderror(mtime,theta(:,k),aO(k),bO(k),'Linear');
            end
            
            [theta,ostats] = getShoreOrientation(obj,E,N,bE,bN,time,site);
            
            %package output
            pos = {E,N,xD,BS,theta};
            stats = table(aX,bX,RsqX,seX,aM,bM,RsqM,seM,...
              'VariableNames',{'distIntercept','distSlope','distR2',...
              'distStdErr','slopeIntercept','slopeSlope','slopeR2',...
              'slopeStdErr'});
            stats = horzcat(stats,ostats);
            meta.xmin = xmin; meta.proflist = proflist; meta.npro = npro;
        end
%%
        function [theta,stats] = getShoreOrientation(~,Es,Ns,Eb,Nb,time,site)
            %use the Es,Ns coordinates of shoreline points for a set of profiles
            %to assign an orientation to each profile. the baseline
            %coordinates, Eb,Nb allow the correct direction to be identified
            %returns ts of angles for each profile and rates of change
            theta = NaN(size(Es));  %bE = NaN(1,size(Es,2));  bN = bE;
            shoreang = site.ShorelineAngle;
            for i=1:length(time)
%                 [sE,sN,idist] = sortENdata2line(Es(i,:),Ns(i,:));
%                 bE(idist) = Eb(i,:);
%                 bN(idist) = Nb(i,:);
                temp = shore_orientation(Es(i,:),Ns(i,:),Eb(i,:),Nb(i,:));
                invangle = mod(mean(temp,'omitnan')+180, 360);
                if invangle>shoreang-89 && invangle<shoreang+89
                    temp = mod(temp+180, 360);
                end
                theta(i,:) = temp;
%                 theta(i,:) = shore_orientation(Es(i,:),Ns(i,:),Eb(i,:),Nb(i,:));
            end
            %
            npro = size(theta,2);
            aO = zeros(npro,1); bO = aO; RsqO = aO; seO = aO;
            for k=1:npro
                %now get rate of change in orientation
                %eps(0) to avoid divide by zero in linear regression
                mtime = years(time-time(1)+eps(0)); 
                [aO(k),bO(k),RsqO(k)] = regression_model(mtime,theta(:,k),'Linear');
                seO(k) = stderror(mtime,theta(:,k),aO(k),bO(k),'Linear');
            end
             stats = table(aO,bO,RsqO,seO,...
                         'VariableNames',{'angleIntercept','angleSlope',...
                                          'angleR2','angleStdErr'});
        end
%%
        function [X,Z,xrange,zrange] = getSampledDataPoints(~,x,z,mxmn)
            %subsample the data set within the user defined window (box)
            % x, z define the profiles for each survey (rows are time)
            % mxmn defines the sample window of min and max values
            % X, Z are the profiles within the sample window
            % xrange, zrange the extent of the data for each profile
            nstep = size(x,1);
            xmin = mxmn.xmin;
            zmin = mxmn.zmin;
            X{nstep,1} = []; Z = X; 
            xrange = NaN(nstep,1); zrange = xrange;
            for i=1:nstep
                vx = x(i,:);
                vz = z(i,:);
                %remove any duplicates in x (for interpolation)
                [vx,ix] = uniquetol(vx,1,'DataScale',0.01); 
                vz = vz(ix);
                nrec = length(vx);

                %find points at each end of profile nearest to xmin and zmin
                [~,Ix] = min(abs(vx-xmin));
                [~,Iz] = min(abs(vz(Ix:end)-zmin));
                Iz = Iz+(Ix-1);   %needed to deal with troughs behind crest

                %interpolate to nearest point to origin if there is data
                if vx(Ix)<xmin
                    vz(Ix) = interp1(vx(Ix:Ix+1),vz(Ix:Ix+1),xmin);
                    vx(Ix) = xmin;
                elseif (Ix>1 &&  vx(Ix-1)<xmin) 
                    vz(Ix) = interp1(vx(Ix-1:Ix),vz(Ix-1:Ix),xmin);
                    vx(Ix) = xmin;
                end
                %
                if vz(Iz)<zmin
                    vx(Iz) = interp1(vz(Iz-1:Iz),vx(Iz-1:Iz),zmin);
                    vz(Iz) = zmin;
                elseif (Iz<nrec && ~isnan(vz(Iz+1)) && vz(Iz+1)<zmin)
                    vx(Iz) = interp1(vz(Iz:Iz+1),vx(Iz:Iz+1),zmin);
                    vz(Iz) = zmin;
                end

                %select data points above minimum values
                idv = find(vx>=vx(Ix) & vz>=vz(Iz));

                %distances from the respective offset origins
                xi = vx(idv)-xmin;
                zi = vz(idv)-zmin;                

                %handle nans in data and any negative values
                idx = isnan(xi) | xi<0 | zi<0;
                xi(idx) = 0;
                zi(idx) = 0;
                
                %data range: z value on x-axis and x value on z-axis
                %pre-assigned as NaNs. Nan values indicate that the x value
                %at minimum z is -ve, or the z value at minimum x is -ve
                if (vx(Iz)-xmin)>=0
                    xrange(i) = vx(Iz)-xmin;
                end
                %
                if (vz(Ix)-zmin)>=0
                    zrange(i) = vz(Ix)-zmin;
                end
                
                %store xi and zi for use in normalisation loop
                X{i} = xi;
                Z{i} = zi;                
            end
        end
 %%
         function [m0,x1,z1] = getVolumeCentroids(~,X,Z,xrange,zrange,isND)
            % get the non-dimensional volumes and centroid values for the profile
            % X, Z are the profiles within the sample window
            % xrange, zrange define extent of profile to be sampled
            %   - maximum values for any profile used to normalise data 
            % isND - flag, if true returns non-dimensional values
            % m0 is the volume and x1, z1 are the centroid moments            
            nstep = size(X,1);
            if isND
                maXrange = max(xrange);  %maximum range for any profile
                maZrange = max(zrange);
            else
                maXrange = 1; maZrange = 1;
            end
            m0 = zeros(nstep,1); x1 = m0; z1 = m0;
            for i=1:nstep %normalise x znd z based on range of all profiles
                xn = X{i}/maXrange;
                zn = Z{i}/maZrange;
                if length(xn)<2 || isnan(xrange(i)) || isnan(zrange(i))
                    m0(i) = NaN; x1(i) = NaN; z1(i) = NaN;
                else
                    %if non-dimensional values (multiply by xrange to scale)
                    [x1(i),z1(i),m0(i)] = section_centroid(xn,zn);
                end                
            end
         end
%%
        function [xdist,slope,options] = getShorelinePositions(obj,x,z,t,h_but,h_ax)
            %get the points at a user selected elevations to define the
            %shoreline position
            mxmn = getMaxMin(obj,x,z);
            xmin = mxmn.xmin; xmax = mxmn.xmax;
            
            prompt = {'Shoreline elevation (mOD):','Slope range:', ...
                'Start distance (m):'};   
            title = 'Define input parameters';
            numlines = 1;
            defaultvalues{1} = num2str(0);
            defaultvalues{2} = num2str(0.5);
            defaultvalues{3} = num2str(xmin);
            h1 = [];
            ok=0;
            hold on
            h2 = plot(h_ax,[xmin xmax],[0 0],'--r','DisplayName','z0');
            h3 = plot(h_ax,[xmin xmax],[0.5 0.5],':r','DisplayName','z0+offset');
            h4 = plot(h_ax,[xmin xmax],[-0.5 -0.5],':r','DisplayName','z0-offset');
            hold off
            datelist = datestr(t,'dd-mmm-yy');
            if length(datelist)<20
                legend(datelist);
            end
            nstep = size(x,1);
            %get input from user
            while ok<1
                useInp=inputdlg(prompt,title,numlines,defaultvalues);
                if isempty(useInp), return; end %user cancelled
                zlevel = str2double(useInp{1});
                delz = str2double(useInp{2});
                xmin = str2double(useInp{3});                               
                [xdist,slope,xs,zs,] = slope_points(x,z,zlevel,delz,xmin);
                delete(h1)
                delete(h2)
                delete(h3)
                delete(h4)
                
                hold on
                h2 = plot(h_ax,[xmin,xmax],[zlevel,zlevel],'--r','DisplayName','z0');
                h3 = plot(h_ax,[xmin xmax],[zlevel+delz zlevel+delz],':r','DisplayName','z0+offset');
                h4 = plot(h_ax,[xmin xmax],[zlevel-delz zlevel-delz],':r','DisplayName','z0-offset');
                for i=1:nstep
                    h1(i) = plot(h_ax,xs(:,i),zs(:,i),'*'); %#ok<AGROW>
                    % Exclude line from legend
                    set(get(get(h1(i),'Annotation'),'LegendInformation'),...
                                                'IconDisplayStyle','off'); 
                end
                hold off  
                
                waitfor(h_but,'Tag');
                if ~ishandle(h_but) %this handles the user deleting figure window    
                    ok = 0;
                elseif strcmp(h_but.Tag,'Yes')
                    ok = 1; %box definition accepted                    
                else
                    %Prompt user to change definition
                    defaultvalues{1} = num2str(zlevel);
                    defaultvalues{2} = num2str(delz);
                    defaultvalues{3} = num2str(xmin);
                    h_but.Tag = '';                     
                end 
                options.zlevel = zlevel; 
                options.delz = delz;
                options.xmin = xmin;
            end
        end
%%
        function [Es,Ns,Eb,Nb] = ENofShoreline(~,dst,xch)
            %find the eastings and northings from a timeseries of shoreline chainage
            %positions for a single beach profile. returns the shoreline
            %coordinates and the baseline coordinates.
            %- used in getPositionAndRates
            x = dst.Chainage;
            E = dst.Eastings;
            N = dst.Northings;  
            nrec = size(x,1);
            Es = zeros(nrec,1); Ns = Es; Eb = Es; Nb = Es;
            for i=1:nrec 
                xi = x(i,:);
                Ei = E(i,:);
                Ni = N(i,:);
                nanflags = isnan(Ei) | isnan(Ni);
                Ei(nanflags) = [];
                Ni(nanflags) = [];
                Etot = Ei(end)-Ei(1);
                Ntot = Ni(end)-Ni(1);
                offset = xi(1);  %the first point of the profile may not be at the baseline 
                if Etot<0, sgn = -1; else, sgn=1; end
                alp = atan(Ntot/Etot);
                Es(i,1) = Ei(1)+sgn*(xch(i)-offset)*cos(alp); %shoreline
                Ns(i,1) = Ni(1)+sgn*(xch(i)-offset)*sin(alp);
                [~,idminCh] = min(xi(:,1),[],'omitnan'); %minimum in first column
                Eb(i,1) = Ei(idminCh,1);
                Nb(i,1) = Ni(idminCh,1);
                
                
%                 Eb(i,1) = Ei(1);                              %baseline
%                 Nb(i,1) = Ni(1);
            end
        end
%%
        function [sortedE,sortedN,idd] = shorelineProfileOrder(~,mobj,caserec)
            %sort the profile order based on the E,N of the base point - min(Chainage)
            npro = length(caserec);   %number of profiles
            Es = NaN(1,npro); Ns = Es;
            for k=1:npro
                [dst,~] = getDataset(mobj.Cases,caserec(k),1);  %idset=1
                E = dst.Eastings;
                N = dst.Northings; 
                Ch = dst.Chainage;
                [~,idminCh] = min(Ch(:,1),[],'omitnan'); %minimum in first column
                Es(1,k) = E(idminCh,1);
                Ns(1,k) = N(idminCh,1); 
            end
            %sort profiles into alongshore order
            [sortedE,sortedN,idd] = sortENdata2line(Es,Ns);
        end

%--------------------------------------------------------------------------
% Plot Functions called by beach models
%--------------------------------------------------------------------------
        function plotBeachProfiles(~,h_ax,x,z)
            %plot the profiles
            xlabel('Chainage (m)');
            ylabel('Elevation (mOD');
            nstep = size(x,1);
            plot(h_ax,x(1,:),z(1,:));
            hold on
            for i=2:nstep
                plot(h_ax,x(i,:),z(i,:));
            end
            hold off      
        end
%%
        function plotShoreline_EN(obj,h_ax,time,Eastings,Northings,pid)   
            %plot shoreline positions as Eastings and Northings over time
            %create in axes handle, h_ax, a tab or standalone figure
            %Loc - list of profiles used                     
            nline = length(time);
            hold on
            for i=1:nline
                Ei = Eastings(i,:);
                Ni = Northings(i,:);
                [sortedE,sortedN] = sortENdata2line(Ei,Ni);
                %now update plot with a new line
                plot(h_ax,sortedE,sortedN,'DisplayName',time{i},...
                                       'ButtonDownFcn',@nameDisplay);                
            end
            if nline<obj.leglimit
                legend;
            end
            title(pid);
            h_ax.XLabel.String = 'Eastings (m)';
            h_ax.YLabel.String = 'Northings (m)';
            hold off
            set(groot,'defaultAxesColorOrder','remove') 
            
        end    
%%
        function plotShoreline_CH(obj,h_ax,time,Location,Variable,pid...
                                                            ,varname,nint)
            %plot variable as a function of profile location over time
            %create in axes handle, h_ax, a tab or standalone figure
            %pid is the selected shoreline dataset description
            %nint plots the shorelines at profile points at nint intervals
            nline = size(Variable,1);
            xints = 1:nint:length(Location);
            nint = length(xints);
            Xi = 1:nint;
            hold on
            for i=1:nline
                y = Variable(i,:);
                Yi = y(xints);
                %now update plot with a new line
                plot(h_ax,Xi,Yi,'.-','DisplayName',time{i},...
                                       'ButtonDownFcn',@nameDisplay);                
            end
            if nline<obj.leglimit
                legend;
            end
            title(pid);
            h_ax.XLabel.String = 'Location';
            h_ax.YLabel.String = varname;
            h_ax.XTick = 1:nint;
            h_ax.XTickLabel = Location(xints);
            h_ax.XTickLabelRotation = 45;
            hold off
            set(groot,'defaultAxesColorOrder','remove') 
        end
%%
        function plotShoreRates_CH(~,h_ax,Location,outtable,pid,nint)
            %plot rates of change as a function of profile location
            xints = 1:nint:length(Location);
            nint = length(xints);
            Xi = 1:nint;
            Yi = outtable.slopeSlope(Xi);
            Ei = outtable.slopeStdErr(Xi);
            
            s1 = subplot(2,1,1,h_ax);
            bar(Xi,Yi)
            hold on 
            er = errorbar(Xi,Yi,Ei,Ei);    
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            hold off
            s1.XTickLabel = {};
            distunits = outtable.Properties.VariableUnits{6};
            s1.YLabel.String = sprintf('Slope rate of change (%s)',distunits);
            s1.Position(2) = 0.56; s1.Position(4) = 0.36;
            legend('Rate of change','Std. Error')
            
            Yi = outtable.distSlope(Xi);
            Ei = outtable.distStdErr(Xi);
            s2 = subplot(2,1,2);
            bar(Xi,Yi)
            hold on 
            er = errorbar(Xi,Yi,Ei,Ei);    
            er.Color = [0 0 0];                            
            er.LineStyle = 'none'; 
            hold off            
            s2.XLabel.String = 'Location';            
            s2.XTick = 1:nint;
            s2.XTickLabel = Location(xints);
            s2.XTickLabelRotation = 45;   
            distunits = outtable.Properties.VariableUnits{2};
            s2.YLabel.String = sprintf('Rate of change (%s)',distunits);
            s2.Position(2) = 0.18; s2.Position(4) = 0.36;
            sgtitle(pid);
        end
%%
        function plotBVI_CH(~,h_ax,Location,outtable,pid,nint)
            %function to plot results of BVI analysis 
            xints = 1:nint:length(Location);
            nint = length(xints);
            Xi = 1:nint;
            fnames = outtable.Properties.VariableNames;
            marker = {'x-','d-','*-','s-','+-'};
            hold on            
            for i=1:width(outtable)-1
                Yi = outtable{:,i};
                plot(h_ax,Xi,Yi,marker{i},'MarkerFaceColor','auto',...
                                'DisplayName',fnames{i});
            end
            Yi = outtable{:,end};
            stem(h_ax,Xi,Yi,'MarkerFaceColor',[0.4660 0.6740 0.1880],...
                            'MarkerEdgeColor',[0.9290 0.6940 0.1250],...
                            'DisplayName',fnames{end});
            hold off
            h_ax.XLabel.String = 'Location';            
            h_ax.XTick = 1:nint;
            h_ax.XTickLabel = Location(xints);
            h_ax.XTickLabelRotation = 45;  
            h_ax.YLabel.String = 'Index';
            legend('Location','best')
            title(pid);
        end
%%
        function plotProfileBaseline(~,mobj,ax,Loc)
            %use the selevcted profiles to generate a base plot
            muicat = mobj.Cases;
            lobj = muicat.DataSets.ctBeachProfileData;  %handle to beach profile data
            %convert class selection to caserec values
            nselect = length(Loc);
            caserec = zeros(1,nselect);
            for i=1:nselect
                caserec(i) = find(strcmp(muicat.Catalogue.CaseDescription,Loc{i}));
            end
            plotProfileLocations(lobj,muicat,ax,caserec,true);
        end
%%
        function plotSpaceTime(~,time,pnames,data,varname)
            %plot the variation of a variable in time for multiple profiles
            %NB uses colormap generated using cbrewer.m and calls
            %inpaint_name.m
            xint = length(time);
            yint = length(pnames);
            xtext = 'Time';
            ytext = 'Profile number';
            legtxt = varname;
            titletxt = sprintf('Surface plot of %s',varname);
            ptype = 'contourf'; 
            %apply function from Forum to replace NaNs with interpolated values
            method = 3;  %results can be sensitive to fitting method used
            zt = inpaint_nans(data,method)'; 
            yt = 1:yint;
            
            %to normalise data use
            isnorm = questdlg('Normalise data?','Space-time plot',...
                                'No','Space','Time','No');
            if ~strcmp(isnorm,'No') 
                dim = 1; %use 1 to normalise over space and 2 over time
                if strcmp(isnorm,'Time')
                    dim = 2;
                end
                zt = ScaleVariable(zt,dim);
            end
            
            %reverse order of profiles if required
            qstr = sprintf('Profiles are from %s to %s',pnames{1},pnames{end});                
            revbut = questdlg(qstr,'Reverese order','Use','Reverse','Use');
            if strcmp(revbut,'Reverse')
                pnames = flipud(pnames);
                zt = flipud(zt);
            end        

            %convert time for surface plotting
            time = datetime(time);
            startyear = year(time(1)); 
            xt = startyear+years(time-datetime(startyear,1,1));   

            hfig = figure('Name','Results Plot','Units','normalized',...                
                          'Resize','on','HandleVisibility','on','Tag','PlotFig');
                      
            plotdir = questdlg('Time axis?','Space time plot','X','Y','X');
            if strcmp(plotdir,'X')
                muiPlots.xySurface(xt,yt,zt,xint,yint,xtext,ytext,...
                    legtxt,titletxt,ptype);
                figax = hfig.CurrentAxes; 
                hold(figax,'on')   
                figax.YTick = yt;
                figax.YTickLabel = pnames;
            else
                %to plot profiles on the x-axis use
                muiPlots.xySurface(yt,xt,zt',yint,xint,ytext,xtext,...
                    legtxt,titletxt,ptype);
                figax = hfig.CurrentAxes; 
                hold(figax,'on')   
                figax.XTick = yt;
                figax.XTickLabel = pnames;           
                figax.XTickLabelRotation = 45;
            end
            hold(figax,'off')
            %to fix the range of colorbar use
%             caxis([200,500]);  
            
            function var = ScaleVariable(inputvar,dim)
                vardims = fliplr(size(inputvar));
                %invert matrix if dim=2 to obtain values for each row 
                if dim==2
                    inputvar = inputvar';
                end
                for i = 1:vardims(dim)
                    subvar = inputvar(:,i);              
                    mvar = mean(subvar,'omitnan');
                    svar = std(subvar,'omitnan');
                    var(:,i) = (subvar-mvar)/svar;
                end
                %restore matrix if dim=2
                if dim==2
                    var = var';
                end
            end
        end    
    
%%
%--------------------------------------------------------------------------
% Define DSproperties for the various models
%--------------------------------------------------------------------------
        function dsp = modelDSproperties(~,id_model) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            switch id_model
                case 1              %Volumes
                    dsp.Variables = struct(...                       
                        'Name',{'Vol','m0','x1','z1','BS'},...
                        'Description',{'Beach volume','N_D zero moment',...
                                    'N_D x-centroid','N-D z-centorid',...
                                   'Beach slope (1:bs)'},...
                        'Unit',{'m^3','-','-','-','-'},...
                        'Label',{'Volume (m^3/m-run)','N-D volume',...
                                   'N-D x-centroid','N-D z-centroid',...
                                   'Beach slope (1:bs)'},...
                        'QCflag',repmat({'model'},1,5));
                case 2              %Shoreline position at a profile
                    dsp.Variables = struct(...                       
                        'Name',{'Chainage','Slope'},...
                        'Description',{'Shoreline position','Beach slope (1:bs)'},...
                        'Unit',{'m','-'},...
                        'Label',{'Shoreline position (m)','Beach slope (1:bs)'},...
                        'QCflag',{'model','model'});
                case 3              %Shoreline as a vector
                    dsp.Variables = struct(...                       
                        'Name',{'Eastings','Northings','Chainage',...
                                   'Slope','SLang'},...
                        'Description',{'Eastings','Northings','Shoreline position',...
                                   'Shoreline slope','Shoreline angle'},...
                        'Unit',{'m','m','m','-','degTN'},...
                        'Label',{'Easting (m)','Northing (m)',...
                                   'Shoreline position (m)',...
                                   'Shoreline slope (1:bs)',...
                                   'Shoreline angle (degTN)'},...
                        'QCflag',repmat({'model'},1,5));
                case 4              %Beach Vulnerability Index
                    dsp.Variables = struct(...                       
                        'Name',{'Eastings','Northings','Chainage',...
                                   'Slope','SLang'},...
                        'Description',{'Eastings','Northings','Shoreline position',...
                                   'Shoreline slope','Shoreline angle'},...
                        'Unit',{'m','m','m','-','degTN'},...
                        'Label',{'Easting (m)','Northing (m)',...
                                   'Shoreline position (m)',...
                                   'Shoreline slope (1:bs)',...
                                   'Shoreline angle (degTN)'},...
                        'QCflag',repmat({'model'},1,5));
            end
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});    
        end
    end  
end