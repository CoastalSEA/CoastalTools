classdef CoastalTools < muiModelUI                        
%
%-------class help---------------------------------------------------------
% NAME
%   CoastalTools.m
% PURPOSE
%   Main GUI for CoastalTools interface, which implements the 
%   muiModelUI abstract class to define main menus.
% SEE ALSO
%   Abstract class muiModelUI.m and tools provided in muitoolbox
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
% 
    properties  (Access = protected)
        %implement properties defined as Abstract in muiModelUI
        vNumber = '1.0'
        vDate   = 'Jan 2021'
        modelName = 'CoastalTools'                        
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = CoastalTools                        
            %constructor function initialises GUI
            obj = setMUI(obj);             
        end
    end
%% ------------------------------------------------------------------------
% Definition of GUI Settings
%--------------------------------------------------------------------------  
    methods (Access = protected)
        function obj = setMUI(obj)
            %initialise standard figure and menus    
            modelLogo = 'CoastalTools_logo.jpg';  %default splash figure
            %classes required to run model, format:           
            %obj.ModelInputs.<model classname> = {'Param_class1',Param_class2',etc}
            obj.ModelInputs.InWaveModel = {'ctWaveParameters','ctWaveData'};  
            obj.ModelInputs.OffWaveModel = {'ctWaveParameters','ctWaveData'};
            obj.ModelInputs.WindWaveModel = {'ctHindcastParameters','ctWindData'};  
            obj.ModelInputs.TidalAnalysis = {'ctWaterLevelData'};
            obj.ModelInputs.CT_WaveModels = {'ctWaveParameters','ctWaveModel'};  
            obj.ModelInputs.CT_BeachAnalysis = {'ctBeachProfileData'}; 
            obj.ModelInputs.CT_UserModel = {'ctWaveParameters'}; 
            
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'2D','3D','4D','2DT','3DT','4DT'}; %,'Profiles'};  
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Timeseries','Taylor','Intervals'};              
            
            initialiseUI(obj,modelLogo); %initialise menus and tabs                  
        end    
        
%% ------------------------------------------------------------------------
% Definition of Menu Settings
%--------------------------------------------------------------------------
        function menu = setMenus(obj)
            %define top level menu items and any submenus
            %MenuLabels can any text but should avoid these case-sensitive 
            %reserved words: "default", "remove", and "factory". If label 
            %is not a valid Matlab field name this the struct entry
            %is modified to a valid name (eg removes space if two words).
            %The 'gcbo:' Callback text triggers an additional level in the 
            %menu. Main menu labels are defined in sequential order and 
            %submenus in order following each brach to the lowest level 
            %before defining the next branch.         
                                                              % << Edit menu to suit model 
            MenuLabels = {'File','Tools','Project','Setup','Run',...
                                                        'Analysis','Help'};
            menu = menuStruct(obj,MenuLabels);  %create empty menu struct
            %
            %% File menu --------------------------------------------------
             %list as per muiModelUI.fileMenuOptions
            menu.File.List = {'New','Open','Save','Save as','Exit'};
            menu.File.Callback = repmat({@obj.fileMenuOptions},[1,5]);
            
            %% Tools menu -------------------------------------------------
            %list as per muiModelUI.toolsMenuOptions
            menu.Tools(1).List = {'Refresh','Clear all'};
            menu.Tools(1).Callback = {@obj.refresh, 'gcbo;'};  
            
            % submenu for 'Clear all'
            menu.Tools(2).List = {'Model','Figures','Cases'};
            menu.Tools(2).Callback = repmat({@obj.toolsMenuOptions},[1,3]);

            %% Project menu -----------------------------------------------
            menu.Project(1).List = {'Project Info','Cases','Export/Import'};
            menu.Project(1).Callback = {@obj.editProjectInfo,'gcbo;','gcbo;'};
            
            %list as per muiModelUI.projectMenuOptions
            % submenu for Scenarios
            menu.Project(2).List = {'Edit Description','Edit Data Set',...
                                    'Save Data Set','Delete Case','Reload Case',...
                                    'View Case Settings'};                                               
            menu.Project(2).Callback = repmat({@obj.projectMenuOptions},[1,6]);
            
            % submenu for 'Export/Import'                                          
            menu.Project(3).List = {'Export Case','Import Case'};
            menu.Project(3).Callback = repmat({@obj.projectMenuOptions},[1,2]);
            
            %% Setup menu -------------------------------------------------
            menu.Setup(1).List = {'Import Data','Site parameters',...
                                      'Model parameters','Data clean-up'};                                    
            menu.Setup(1).Callback = repmat({'gcbo;'},[1,4]);
            menu.Setup(1).Separator = {'off','off','off','on'}; %separator preceeds item
            
            % submenu for Import Data
            menu.Setup(2).List = {'Waves','Water levels','Winds',...
                                  'Beach profiles','Shorelines', ...
                                  'BlueKenue data','User dataset'}; 
            nitems = length(menu.Setup(2).List);
            menu.Setup(2).Callback = repmat({'gcbo;'},[1,nitems]);
            menu.Setup(2).Separator = {'off','off','off','off','off','on','off'};
            
            for j=1:nitems  %add standard submenu to all import menu items
            menu.Setup(j+2).List = {'Load','Add','Delete','Quality Control'};
            menu.Setup(j+2).Callback = repmat({@obj.loadMenuOptions},[1,4]);
            end
            % submenu for Site parameters
            offset = nitems+3;
            menu.Setup(offset).List = {'Wave propagation',...
                            'Wind-wave hindcast','Structure parameters'};
            menu.Setup(offset).Callback = repmat({@obj.runProps},[1,3]);
            % submenu for Model parameters
            menu.Setup(offset+1).List = {'YGOR simulation parameters',...
                            'BMV simulation parameters','Model constants'};                                   
            menu.Setup(offset+1).Callback = repmat({@obj.runProps},[1,3]);    
            % submenu for Data clean-up
            menu.Setup(offset+2).List = {'Concatenate two timeseries',...
                            'Resample timeseries','Patch timeseries',...                
                            'Trim timeseries','Delete multiple profiles',...                
                            'Edit or Delete profile in timeseries'};
            menu.Setup(offset+2).Callback = repmat({@obj.datacleanup},[1,6]);
            menu.Setup(offset+2).Separator = {'off','off','off','off','on','off'};
            
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Wave properties','Beach properties',...
                                'Tidal analysis','Derive Output',... 
                                'Simulation','Vulnerability','User Model'};
            menu.Run(1).Callback = {'gcbo;','gcbo;','gcbo;',@obj.runModel,...
                                    @obj.runModel,'gcbo;',@obj.runModel};
            menu.Run(1).Separator = {'off','off','off','on','on','off','on'};
            
            % submenu for Wave properties
            menu.Run(2).List = {'Deepwater waves','Nearshore Waves',...
                                'Wind-Waves','Wave Energy','Runup',...
                                'Littoral Drift','X-shore Transport',...
                                'Overtopping','Iribarren Number'}; 
            nitems = length(menu.Run(2).List);
            menu.Run(2).Callback = repmat({@obj.runWaveModel},[1,nitems]);
            
            % submenu for Beach properties
            menu.Run(3).List = {'Profiles','Shore change','Beach type',...
                                           'Shore profile','Dean profile'}; 
            menu.Run(3).Callback = {'gcbo;','gcbo;',@obj.runBeachAnalysis...
                            @obj.runBeachAnalysis,@obj.runBeachAnalysis};
            menu.Run(3).Separator = {'off','off','off','on','off'};
            
            menu.Run(4).List = {'Volumes','Shoreline position',...
                        'Location plot','Centroid plot','Space-time plot'}; 
            menu.Run(4).Callback = repmat({@obj.runBeachAnalysis},[1,5]);
            menu.Run(4).Separator = {'off','off','on','off','off'};
            
            menu.Run(5).List = {'Shoreline','Change plot','Rates plot'};
            menu.Run(5).Callback = repmat({@obj.runBeachAnalysis},[1,3]);
            menu.Run(5).Separator = {'off','on','off'};
            
            % submenu for Tidal analysis
            menu.Run(6).List = {'Analysis','Reconstruction'}; 
            menu.Run(6).Callback = repmat({@obj.runTides},[1,2]);
            
            % submenu for Beach vulnerability
            menu.Run(7).List = {'BVI site','BVI profile set','BVI set plot'};
            menu.Run(7).Callback = repmat({@obj.runModel},[1,3]);            
            
            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics'};
            menu.Analysis(1).Callback = repmat({@obj.analysisMenuOptions},[1,2]);
            
            %% Help menu --------------------------------------------------
            menu.Help(1).Callback = {@obj.Help}; %make model specific?
            
        end
        
%% ------------------------------------------------------------------------
% Definition of Tab Settings
%--------------------------------------------------------------------------
        function [tabs,subtabs] = setTabs(obj)
            %define main tabs and any subtabs required. struct field is 
            %used to set the uitab Tag (prefixed with sub for subtabs). 
            %Order of assignment to struct determines order of tabs in figure.
            %format for tabs: 
            %    tabs.<tagname> = {<tab label>,<callback function>};
            %format for subtabs: 
            %    subtabs.<tagname>(i,:) = {<subtab label>,<callback function>};
            %where <tagname> is the struct fieldname for the top level tab.
            tabs.Cases  = {'   Data  ',@obj.refresh};        
            tabs.Models = {'  Models  ',@obj.refresh};
            tabs.Site = {'  Site  ','gcbo;'};
            subtabs.Site(1,:) = {'  Waves  ',@obj.InputTabSummary};
            subtabs.Site(2,:) = {' Simulation ',@obj.InputTabSummary};
            tabs.Plot = {'  Q-Plot  ',@obj.getTabData};
            tabs.Calcs = {'  Calcs  ','gcbo;'};
            subtabs.Calcs(1,:) = {' Volumes ',@obj.InputTabSummary};
            subtabs.Calcs(2,:) = {' Shoreline ',@obj.InputTabSummary};
            subtabs.Calcs(3,:) = {' Profile ',@obj.InputTabSummary};
            tabs.Stats = {'   Stats   ','gcbo;'};
            subtabs.Stats(1,:) = {' Descriptive ',@obj.InputTabSummary};
            subtabs.Stats(2,:) = {' Extremes ',@obj.InputTabSummary};
        end
       
%%
        function props = setTabProperties(~)
            %define the tab and position to display class data tables
            %props format: {class name, tab tag name, position, ...
            %               column width, table title}
            % position and column widths vary with number of parameters
            % (rows) and width of input text and values. Inidcative
            % positions:  top left [0.95,0.48];    top right [0.95,0.97]
            %         bottom left [0.45, 0.48]; bottom rigth [0.45,0.97]
                                                             % << Edit input properties classnames 
            props = {...                                     % << Add additional inputs and adjust layout
                'ctWaveParameters','Waves',[0.95,0.48],{180,60},'Site data:';...
                'ctHindcastParameters','Waves',[0.95,0.97],{180,60},'Wind-wave parameters:';...
                'ctStructureInput','Waves',[0.56,0.97],{180,60},'Structure definition:';...
                'Sim_YGORinput','Simulation',[0.95,0.48],{180,60},'YGOR Forecast parameters:';...
                'Sim_BMVinput','Simulation',[0.95,0.98],{180,60},'BMV Model parameters:'};
        end    
 %%
        function setTabAction(~,src,cobj)
            %function required by GUIinterface and sets action for selected
            %tab (src)
            switch src.Tag                                    % << Edit match tab requirements
                case 'Plot' 
                     tabPlot(cobj,src);
                case 'Stats'
                     tabStats(cobj,src);    
            end
        end      
%% ------------------------------------------------------------------------
% Callback functions used by menus and tabs
%-------------------------------------------------------------------------- 
        %% File menu ------------------------------------------------------
        %use default menu functions defined in muiModelUI
            
        %% Tools menu -----------------------------------------------------
        %use default menu functions defined in muiModelUI
                
        %% Project menu ---------------------------------------------------
        %use default menu functions defined in muiModelUI           

        %% Setup menu -----------------------------------------------------
        function setupMenuOptions(obj,src,~)
            %callback functions for data input
            switch src.Text
                case 'Input parameters'                       % << Edit to call Parameter Input class
                    ParamInput_template.setParamInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Inputs');
                    InputTabSummary(obj,tabsrc);
                case 'Run parameters'                         % << Edit to call Data Import class
                    ParamInput_template.setParamInput(obj);  
                    %update tab display with input data
                    tabsrc = findobj(obj.mUI.Tabs,'Tag','Inputs');
                    InputTabSummary(obj,tabsrc);
                case 'Model Constants'
                    obj.Constants = editProperties(obj.Constants);
            end
        end  
%%
        function loadMenuOptions(obj,src,~)
            %callback functions to import data   
            mode = 'single';
            switch src.Parent.Text
                case 'Waves'
                    classname = 'ctWaveData';
                case 'Water levels'
                    classname = 'ctWaterLevelData';
                case 'Winds'
                    classname = 'ctWindData';
                case 'Beach profiles'
                    classname = 'ctBeachProfileData';  
                    mode = 'none';
                case 'Shorelines'
                    classname = 'ctShorelineData';  
                case 'BlueKenue data'
                    classname = 'ctBlueKenueData';  
                case 'User dataset'
                    classname = 'muiUserData';
            end
            %
            switch src.Text
                case 'Load'
                    fname = sprintf('%s.loadData',classname);
                    callStaticFunction(obj,classname,fname); 
                case 'Add'
                    useCase(obj.Cases,mode,classname,'addData');
                case 'Delete'
                    useCase(obj.Cases,mode,classname,'deleteData');
                case 'Quality Control'
                    useCase(obj.Cases,mode,classname,'qcData');
            end
            DrawMap(obj);
        end        
%%
        function runProps(obj,src,~)
            %callback functions for run and settings parameter input
            tabname = [];
            switch src.Text
                case 'Wave propagation'
                    ctWaveParameters.setParamInput(obj);
                    tabname = 'Waves';
                case 'Wind-wave hindcast'
                    ctHindcastParameters.setParamInput(obj);
                    tabname = 'Waves';
                case 'Structure parameters'
                    ctStructureInput.setParamInput(obj);
                    tabname = 'Waves';
                case 'YGOR simulation parameters'
                    Sim_YGORinput.setParamInput(obj);
                    tabname = 'Simulation';
                case 'BMV simulation parameters'
                    Sim_BMVinput.setParamInput(obj);
                    tabname = 'Simulation';
                case 'Model constants'
                    obj.Constants = editProperties(obj.Constants);
            end
            %
            if ~isempty(tabname)
                %update tab used for properties
                tabsrc = findobj(obj.mUI.Tabs,'Tag',tabname);
                InputTabSummary(obj,tabsrc);
            end                  
        end
%%
        function datacleanup(obj,src,~)
            %all cleanup options call the same function
            ct_data_cleanup(obj.Cases,src);
            DrawMap(obj);
        end
        %% Run menu -------------------------------------------------------
        function runModel(obj,src,~)
            %callback to run generic model functions
            switch src.Text
                case 'Derive Output'
                    obj.GuiManip = DataManip.getDataManipGui(obj);
                case 'Simulation'
                    obj.GuiSimulation = CT_Simulation.getCTSimGui(obj);
                case 'BVI site'  %vulnerability for single location
                    ct_sitevulnerability(obj);
                case 'BVI profile set'     %vulnerability for set of profiles
                    obj.h_CT_BeachAnalysis = CT_BeachAnalysis.runModel(obj,src);
                case 'BVI set plot'     %plots for vulnerability index
                    [isvalid,id_class,lobj] = check4instance(obj,5);
                    if ~isvalid, return; end  
                    getShoreTablePlot(lobj(id_class),obj,src.Text);
                case 'User Model'
                    %run wave model to create inshore wave data
                    callClassFunction(obj,'CT_UserModel','runUserModel');
            end
        end
%%
        function runWaveModel(obj,src,~)
            %callback functions to run wave models
            switch src.Text
                case 'Nearshore Waves'
                    %run wave model to create inshore wave data
                    callClassFunction(obj,'InWaveModel','runInWaveModel');
                case 'Deepwater waves'
                    %run wave model to create inshore wave data
                    callClassFunction(obj,'OffWaveModel','runOffWaveModel');
                case 'Wind-Waves'
                    %run wave model to create inshore wave data
                    callClassFunction(obj,'WindWaveModel','runWindWaveModel');    
                otherwise
                    obj.h_CT_WaveModels = CT_WaveModels.runModel(obj,src);   
            end
        end
%%
        function runTides(obj,src,~)
            %callback functions to run tidal analysis and reconstruction
            switch src.Text
                case 'Analysis'
                    %run tidal analysis to generate set of constituents
                    callClassFunction(obj,'TidalAnalysis','runTidalAnalysis');
                case 'Reconstruction'
                    %run utide to get reconstructed tidal signal
                    callClassFunction(obj,'TidalAnalysis','getTidalData');
            end
        end
%%
        function runBeachAnalysis(obj,src,~)
            %callback to run beach analysis functions
            %indices refer to model type in CT_BeachAnalysis
            %'Beach type','Volumes','Shoreline position','Shoreline',
            %'BVI profile set','Shore profile','Dean profile'
            switch src.Text    
                case 'Change plot'
                    [isvalid,id_class,lobj] = check4instance(obj,[4,5]);
                    if ~isvalid, return; end 
                    getShorelinePlot(lobj(id_class(1)),obj);
                 case 'Rates plot'
                    [isvalid,id_class,lobj] = check4instance(obj,[4,5]);
                    if ~isvalid, return; end 
                    getShoreTablePlot(lobj(id_class),obj,src.Text); 
                case 'Location plot'                    
                    lobj = obj.Cases.DataSets.ctBeachProfileData;
                    if isempty(lobj) %no data
                        warndlg('No data avaialable to plot')
                        return; 
                    end                     
                    getProfileLocationsPlot(lobj,obj.Cases)
                case 'Centroid plot'
                    [isvalid,id_class,lobj] = check4instance(obj,2);
                    if ~isvalid, return; end 
                    getProfileCentroidsPlot(lobj(id_class),obj);
                case 'Space-time plot'
                    [isvalid,id_class,lobj] = check4instance(obj,[1,2,3]);
                    if ~isvalid, return; end 
                    getProfileSpaceTimePlot(lobj(id_class),obj);
                otherwise
                    obj.h_CT_BeachAnalysis = CT_BeachAnalysis.runModel(obj,src);                     
            end
        end
%%
        function [isvalid,id_class,lobj] = check4instance(obj,ids)
            %check whether the instances in CT_BeachAnalysis or CT_WaveModels
            %defined by the ids vector exist or not. isvalid=true if exists.
            isvalid = false; id_class = [];
            count = 1;
            lobj = obj.h_CT_BeachAnalysis;
            if ~isempty(lobj)
                for i=1:length(ids)
                    [id_c,isnew] = getClassInstance(lobj,'ModelType',ids(i));
                    if ~isnew && ~isempty(lobj(id_c).mtsc)            
                        isvalid = true;
                        id_class(count) = id_c; %#ok<AGROW>
                        count = count+1;
                    end
                end 
            end
            %
            if ~isvalid
                warndlg('No data avaialable to plot')
            end
        end          
        
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,~,~)
            doc CoastalTools
        end 
%% ------------------------------------------------------------------------
% Overload muiModelUI.MapTable to customise Data and Model Tabs
%--------------------------------------------------------------------------  
        function MapTable(obj,ht)
            %create tables for Data and Model tabs - called by DrawMap
            % load case descriptions
            muicat = obj.Cases;
            idx = tabSubset(obj,ht.Tag);
            caseid = muicat.Catalogue.CaseID(idx);
            casedesc = muicat.Catalogue.CaseDescription(idx);
            caseclass = muicat.Catalogue.CaseClass(idx);

            cdata = {'0','Type','Description of individual cases','','','',''};
            irec = 1;
            for i=1:length(caseid)
                case_id = num2str(caseid(i));
                if ~isfield(muicat.DataSets,caseclass{i}) || ...
                                  isempty(muicat.DataSets.(caseclass{i}))
                    type = 'New';
                else
                    type = caseclass{i};
                    if contains(type,'CT_')
                        type = split(type,'_');
                        type = type{2};
                    elseif strcmpi(type(1:2),'ct')
                        type = type(3:end);
                    end
                end
                %
                dst = getDataset(muicat,i,1);
                range = getVarAttRange(dst,1,'Time');
                if ~isempty(dst.RowRange)
                    reclen = num2str(height(dst.DataTable));
                    stdate = datestr(range{1},'dd-mmm-yyyy'); %use to datestr to control ouput format
                    endate = datestr(dst.RowNames(end),'dd-mmm-yyyy');
                    qualcl = '';
                    if ~isempty(dst.VariableQCflags)
                        qualcl = dst.VariableQCflags{1};
                    end
                    cdata(irec,:) = {case_id,type,char(casedesc{i}),reclen,...
                                                     stdate,endate,qualcl};
                    irec = irec+1;
                end                
            end
            % draw table of case descriptions
            tc=uitable('Parent',ht,'Units','normalized',...
                'CellSelectionCallback',@obj.caseCallback,...
                'Tag','cstab');
            tc.ColumnName = {'ID','Data Type','Data Description','Tstep',...
                'Start','End','QC'};
            tc.RowName = {};
            tc.Data = cdata;
            tc.ColumnWidth = {20 80 180 55 80 80 20};
            tc.RowStriping = 'on';
            tc.Position(3:4)=[0.935 0.8];    %may need to use tc.Extent?
            tc.Position(2)=0.9-tc.Position(4);
        end   
    end
 
end    
    
    
    
    
    
    
    
    
    
    