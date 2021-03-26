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
    end
    
    properties (Transient)
        ModelList = {'Volumes','Shoreline position','Shoreline',...
                     'BVI profile set','Shore profile','Dean profile'}
        VolumeCalcs      %struct for results from volume calcs
        ShorelineCalcs   %struct for results from shoreline position calcs
        leglimit = 12    %limit to number of rows in legend         
    end
    
    methods (Access = private)
        function obj = CT_BeachAnalysis()                    
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj,src)
            %function to run a simple 2D diffusion model
            obj = CT_BeachAnalysis;    
            id_model = find(strcmp(obj.ModelList,src.Text));            
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
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
            switch id_model
                case 1   %Volumes for individual profile
                    output = volumesModel(obj,mobj);
                    if ~isempty(output)
                        obj(id_class).VolumeCalcs{irec} = output.table;
                        %SummaryTable call: mobj, table to display, table Tag 
                        %name, handle to tab or figure to use for table
                        src = findobj(mobj.GuiTabs.Children,'Tag','Volumes');
                        SummaryTable(mobj,output.table,'Volumes',src);
                    end
                case 2  %Shoreline position for individual profile
                    output = shorelinePositionModel(obj,mobj);
                    if ~isempty(output)
                        obj(id_class).ShorelineCalcs{irec} = output.table;
                        src = findobj(mobj.GuiTabs.Children,'Tag','Shoreline');
                        SummaryTable(mobj,output.table,'Shoreline',src);
                    end
                case 3  %Shoreline (position extracted from set of profiles)
                    output = profiles2shoreline(obj,mobj);
                    if ~isempty(output)
                        obj(id_class).ShorelineCalcs{irec} = output.table;
                        src = findobj(mobj.GuiTabs.Children,'Tag','Shoreline');
                        SummaryTable(mobj,output.table,'Shoreline',src);
                    end
                case 4 %beach vulnerability index
                    output = getBVindex(obj,mobj);
                    if ~isempty(output)
                        obj(id_class).ShorelineCalcs{irec} = output.bvitable;
                        src = findobj(mobj.GuiTabs.Children,'Tag','Shoreline');
                        SummaryTable(mobj,output.bvitable,'Shoreline',src);
                    end
                case 5  %plot an idealised profile and closure depths
                    src = findobj(mobj.GuiTabs.Children,'Tag','Profile');
                    bmvProfile(obj,mobj,src);
                    output = [];
                case 6  %plot the Dean profile
                    src = findobj(mobj.GuiTabs.Children,'Tag','Profile');
                    getDeanProfile(obj,mobj,src);
                    output = [];
            end
            if isempty(output) || isempty(output.results{1})
                return; %user cancelled or no results returned
            end       
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(output.results{:},'RowNames',output.modeltime,...
                                                    'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = sprintf('Class %s, using %s',metaclass(obj).Name,...
                                                obj.ModelList{id_model});
            dst.MetaData = output.metatxt;
            %save results
            setDataSetRecord(obj,muicat,dst,'model');
            getdialog('Run complete');
        end
 %%
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end       
    end
%%
%--------------------------------------------------------------------------
% Main functions for CT_BeachAnalysis
%--------------------------------------------------------------------------    
    methods
        function output = volumesModel(obj,mobj)
            %compute the volumes per metre shore (CSA) for a timeseries
            %of beach profiles.
            % INPUTS
            % x,z,t - beach profile data
            % xmin  - lower limit for volume calculation
            % zmin  - lower limit for volume calculation
            % FUNCTION CALLS (external) - section_centroid.m
            output = [];
            [tsc_time,x,z,pid] = selectBeachProfile(obj,mobj);  
            if isempty(tsc_time), return; end
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
            output.modeltime = tsc_time;
            output.metatxt = {sprintf('x from %g to %g; z from %g to %g',...
                mxmn.xmin,mxmn.xmax,mxmn.zmin,mxmn.zmax)};  
            output.source = pid;
            output.proflist = pid;
            
            %data to be written to Calcs Tab
            output.table = table(V,xrange,zrange,m0,X1,Z1,...
              'VariableNames',{'Volume','Xrange','Zrange','m0','X1','Z1'});      
            output.table.Properties.RowNames = cellstr(tsc_time); 
            output.table.Properties.Description = ...
                sprintf('Volume and centroid data for %s',pid);
        end
    end
%%    
%% --------------------------------------------------------------------------
% UI Functions called by beach models
%--------------------------------------------------------------------------
    methods (Access = private)
        function [tsc_time,x,z,pid] = selectBeachProfile(~,mobj)
            %user selects a single beach profile - returns chainage and
            %elevation time series for profile
            promptxt = 'Select beach profile';
            h_prof = getClassHandle(mobj,'BeachProfileData');
            [caserec,~,~,ok] = ScenarioList(mobj.Cases,h_prof,...
            'PromptText',promptxt,'ListSize',[300,100]);
            if ok<1
                warndlg('No selection made');
                tsc_time = []; x = []; z = []; pid = [];
                return;
            end
            [tsc,~] = getCaseDataSet(mobj.Cases,mobj,caserec); 
            tsc_time = datetime(getabstime(tsc));
            x = tsc.Chainage.Data;
            z = tsc.Elevation.Data;
            pid = tsc.Name;
        end
        
    end
    methods (Access = private) 
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
                        'Name',{'Vol',',m0','x1','z1','BS'},...
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