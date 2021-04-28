classdef Sim_YGOR < muiDataSet
%
%-------class help---------------------------------------------------------
% NAME
%   Sim_YGOR.m
% PURPOSE
%   Class description - Class for Sim_YGOR model to run as a CoastalTools app
% USAGE
%   obj = Sim_YGOR.runSimulation(mobj); %mobj is a handle to Main UI
%   obj = Sim_YGOR.runForecast(mobj);
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:  
        TabOutput    %output struct with results from model. Struct fields:
                     % table - values saved for display on tabs
                     % headtxt - source and metatext to be used as table header
    end
    
    properties (Transient)
        UIsel           %structure for the variable selection made in the UI
        UIset           %structure for the plot settings made in the UI
    end
    
    methods (Access = private)
        function obj = Sim_YGOR()                    
            %class constructor
        end
    end      
%%
    methods (Static)                
        function runDataFitting(gobj,mobj)
            %function to run a simple 2D diffusion model
            obj = Sim_YGOR;                            
            %check that the input data has been entered
            msgtxt = ('Shore position data required to run model');
            lobj = getClassObj(mobj,'Inputs','CT_BeachAnalysis',msgtxt);
            if isempty(lobj), return; end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); 
            %--------------------------------------------------------------------------
            % Model code 
            %--------------------------------------------------------------------------
            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;
            
            hw = waitbar(0, 'Loading data. Please wait');
            [wvdst,shdst] = getSimData(obj,muicat);
            if isempty(wvdst)
               warndlg('Failed to find overlapping data sets');
            end
            waitbar(0.4)
            %for some datasets the values are too small to work
            %well in the optimization function. Allow user to scale  
            shdvar = shdst.(shdst.VariableNames{1});
            if max(shdvar)<1
                [sgn,scale] = getYGORscaleInputs(obj);
                if isempty(sgn), return; end
                shd.v = sgn*shdvar*scale;
                if scale>1
                    label = shdst.VariableDescriptions{1};
                    label = sprintf('%s x%d',label,sgn*scale);
                    shdst.VariableDescriptions{1} = label;
                end
            end
            waitbar(0.6)
            %YGOR uses the mean values between profiles to estimate the fit
            %parameters, Function also returns numts which holds the number
            %of records found in each interval. If means have already been 
            %added to profile check they are from selected wave data
            %selected case object and data table
            [shobj,~,~] = getCase(muicat,obj.UIsel(2).caserec);   
            mnwv = checkIntervalData(obj,shobj,shdst,wvdst);
            waitbar(0.8)
            %interpolate wave data (which is usually 1 or 3 hourly) to the
            %time assigned to the beach profile (midnight or noon) and ADD
            %these points to the wave timeseries (YGORmodel uses the points
            %that match in both wave and profile timeseries)
            %OR take the average of a user defined preceding interval
            wvdst = getWavesAtProfile(obj,wvdst,shdst);
            waitbar(1)
            %run selected model 
            isfcast = false;
            coeffs = {};
            close(hw)
            %save('YGORmodeldata.mat','wvi','sh','mnwv','labels','coeffs','isfcast');
            
            [C,rmse,~] = simYGORmodel(wvdst,shdst,mnwv,coeffs,isfcast); 
            if isempty(C), return; end    %user aborted
            %--------------------------------------------------------------------------
            % Save estimated coefficients to Sim_YGORinput
            %--------------------------------------------------------------------------
            answer = questdlg('Update YGOR parameters?','YGOR model','Save','Quit','Quit');
            if strcmp(answer,'Quit'), return; end
            %SAVE to Sim_YGORinput
            lobj = getClassObj(mobj,'Inputs','Sim_YGORinput');
            vals = getProperties(lobj);            
            vals(3:7,1) = num2cell(C');
            vals{8} = rmse;
            mobj.Inputs.Sim_YGORinput = updateProperties(lobj,vals);
            getdialog('Run complete');
        end
%%
        function runSimulation(gobj,mobj)
            %function to run a simple 2D diffusion model
            obj = Sim_YGOR;                      
            dsp = modelDSproperties(obj);
            
            %check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('YGOR simulation parameters need to be specified');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); 
            coeffs = getSimParams(getClassObj(mobj,'Inputs','Sim_YGORinput'));
            %--------------------------------------------------------------------------
            % Model code 
            %--------------------------------------------------------------------------
            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;
            
            hw = waitbar(0, 'Loading data. Please wait');
            [wvdst,shdst] = getSimData(obj,muicat);
            if isempty(wvdst)
               warndlg('Failed to find overlapping data sets');
            end
            waitbar(0.6)

            if ~isempty(shdst)  
                %added to allow Sim_YGORmodel to be called in forecast mode
                %without needing to pass shoreline data 
                nrec = height(shdst.DataTable);
                wvb = zeros(nrec,1);
            else
                shdst = []; wvb = 0;
            end            
            
            waitbar(1)
            isfcast = true;
            close(hw);           
            [tt,ss,mtxt] = simYGORmodel(wvdst,shdst,wvb,coeffs,isfcast);
            answer = questdlg('Save the results?','YGOR model','Save','Quit','Quit');
            if strcmp(answer,'Quit'), return; end
            %--------------------------------------------------------------------------
            % Assign model output to a dstable using the defined dsproperties meta-data
            %--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(ss,'RowNames',tt,'DSproperties',dsp);
            %--------------------------------------------------------------------------
            % Save results
            %--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = mtxt;
            %save results
            setDataSetRecord(obj,muicat,dst,'model');
            getdialog('Run complete');            
        end
%%
        function preLoad(gobj,mobj)
            %extract data and add mean wave data to shoreline position data            
            %get the options selected by the user and id for results
            obj = Sim_YGOR; 
            %check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;
            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;
            
            hw = waitbar(0, 'Loading data. Long timeseries can take a while');
            [wvdst,shdst] = getSimData(obj,muicat);
            if isempty(wvdst)
               warndlg('Failed to find overlapping data sets');
            end
            waitbar(0.4)
            %YGOR uses the mean values between profiles to estimate the fit
            %parameters, Function also returns numts which holds the number
            %of records found in each interval.
            [mnwv,~] = getintervaldata(shdst,wvdst); %takes a long time
            waitbar(0.8)
            %selected case object and data table
            [cobj,~,~] = getCase(muicat,obj.UIsel(2).caserec);   
            dst = cobj.Data.Dataset;  %selected dataset
            waitbar(1)
            nvar = width(dst.DataTable)+1; 
            %use nvar to add count in case more than one dataset loaded
            newname = sprintf('Mean%s%d',wvdst.VariableNames{1},nvar);
            dst = addvars(dst,mnwv,'NewVariableNames',newname);
            dst.VariableDescriptions{nvar} = sprintf('Mean %s',wvdst.VariableDescriptions{1});
            dst.VariableUnits{nvar} = 'm';
            dst.VariableLabels{nvar} = wvdst.VariableLabels{1};
            dst.VariableQCflags{nvar} = 'model';
            cobj.Data.Dataset = dst;
            close(hw)
            delete(obj)
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private)
        function [wvdst,shdst] = getSimData(obj,muicat)
            %extract the data selected by the user
            %returns dstables of the selected data
            ok = 1;
            nvar = length(obj.UIsel);
            %initialise struct used in muiCatalogue.getProperty
            props(nvar) = setPropsStruct(muicat);
            for i=1:nvar
                %get the data for each variable
                if obj.UIsel(i).caserec>0
                    props(i) = getProperty(muicat,obj.UIsel(i),'dstable');
                    if isempty(props(i).data), ok = ok-1; end
                end
            end  

            wvdst = props(1).data;
            shdst = props(2).data;
            %get the overlapping portion for the two data sets
            [wvdst,shdst] = getoverlappingtimes(wvdst,shdst,true);
        end        
%%
        function  [sgn,scale] = getYGORscaleInputs(~)
            %get input from user needed for YGORmodel
            sgn = []; scale = [];
            prompt = {'Sign of variable (use -1 to reverse):',...
                      'Factor to scale variable:'};   
            title = 'Define input parameters';
            numlines = 1;
            defaultvalues{1} = num2str(1);
            defaultvalues{2} = num2str(1);
            useInp = inputdlg(prompt,title,numlines,defaultvalues);
            if isempty(useInp), return; end %user cancelled
            sgn = str2double(useInp{1});
            scale = str2double(useInp{s});
        end       
%%
        function mnwv = checkIntervalData(~,shobj,shdst,wvdst)
            %check if means have already been added to profile 
            %then check they are from selected wave data and get values
            shts_vars = shobj.Data.Dataset.VariableNames;          
            %use the shoreline dataset to see if mean wave values available
            idx = find(contains(shts_vars,'Mean'));
            if ~isempty(idx)                     %means available
                if length(idx)>1
                    %need to select which mean ts to use
                    listxt = shts_vars(idx);
                    answer = listdlg('PromptString','Select mean to use',...
                        'SelectionMode','single','ListString',listxt);
                    idx = idx(answer);
                end
                mnwv = shobj.Data.Dataset.(shts_vars{idx});
                
                %check if based on selected wave data
                % if strcmp(wvmetatxt,shtsc_mean.UserData.MetaData)
                %     mnwv = shtsc_mean.Data; %use existing mean values
                % else                        %derive new mean values
                %     [mnwv,~] = getintervaldata(shdst,wvdst);
                % end
            else                            %derive new mean values
                [mnwv,~] = getintervaldata(shdst,wvdst);
            end
        end
%%
        function newvdst = getWavesAtProfile(~,wvdst,shdst)
            %interpolate waves to time of profile, or get the mean values of
            %the waves for an interval proceeding the position survey dates
            prompt = 'Averaging interval in days (0 uses interpolation)';
            prompt = sprintf('Get values at times of survey dates\n%s',prompt);
            C = inputdlg(prompt,'Averaging interval',1,{'0'});
            if isempty(C) || strcmp(C,'0')
                %interpolate wave data to add values at times of profile
                %extract the wave data values
                wvtime = wvdst.RowNames;
                wvtable = wvdst.DataTable;
                wvvals = wvtable{:,1}; 
                shtime = shdst.RowNames;
                idx = ismember(shtime,wvtime);
                newtime = shtime(~idx);
                if ~isempty(newtime)
                    newvals = interp1(wvtime,wvvals,newtime,'linear','extrap');
                    newdst= dstable(newvals,'RowNames',cellstr(newtime),...
                                        'VariableNames',wvdst.VariableNames);         
                    newvdst = mergerows(wvdst,newdst);                    
                end
            else
                %find mean wave for waveint days before times of profile
                waveint = str2double(C); 
                newvdst = getpreceedingdata(wvdst,shdst,waveint);
            end 
        end  
%%        
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                     
                'Name',{'Chainage'},...
                'Description',{'Shoreline position'},...
                'Unit',{'m'},...
                'Label',{'Shoreline position (m)'},...
                'QCflag',{'model'});
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