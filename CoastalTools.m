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
        vNumber = '3.4'
        vDate   = 'Jan 2024'
        modelName = 'CoastalTools'                        
        %Properties defined in muiModelUI that need to be assigned in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = CoastalTools                        
            %constructor function initialises GUI
            isok = check4muitoolbox(obj);
            if ~isok, return; end
            %
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
            obj.ModelInputs.ctWaveModel = {'ctWaveParameters','ctWaveData'};  
            obj.ModelInputs.ctWindWaveModel = {'ctHindcastParameters','ctWindData'};  
            obj.ModelInputs.ctTidalAnalysis = {'ctWaterLevelData'};
            obj.ModelInputs.CT_WaveModels = {'ctWaveParameters','ctWaveModel'};  
            obj.ModelInputs.CT_BeachAnalysis = {'ctBeachProfileData'}; 
            obj.ModelInputs.CT_UserModel = {'ctWaveParameters'}; 
            obj.ModelInputs.Sim_YGOR = {'Sim_YGORinput'}; 
            obj.ModelInputs.Sim_BMV = {'Sim_BMVinput','CT_BeachAnalysis'};
            obj.ModelInputs.WRM_WaveModel = {''};
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'Timeseries','Profiles','Rose','2D','3D','2DT','3DT'}; %,'Profiles'};  
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
             %list as per callback function muiModelUI.fileMenuOptions
            menu.File.List = {'New','Open','Save','Save as','Exit'};
            menu.File.Callback = repmat({@obj.fileMenuOptions},[1,5]);
            
            %% Tools menu -------------------------------------------------
            %overlaod callback function toolsMenuOptions
            menu.Tools(1).List = {'Refresh','Clear all'};
            menu.Tools(1).Callback = {@obj.refresh, 'gcbo;'};  
            
            % submenu for 'Clear all'
            menu.Tools(2).List = {'Project','Figures','Data','Models'};
            menu.Tools(2).Callback = repmat({@obj.toolsMenuOptions},[1,4]);

            %% Project menu -----------------------------------------------
            menu.Project(1).List = {'Project Info','Cases','Export/Import'};
            menu.Project(1).Callback = {@obj.editProjectInfo,'gcbo;','gcbo;'};
            
            %list as per muiModelUI.projectMenuOptions
            menu.Project(2).List = {'Edit Description','Edit DS properties','Edit Data Set',...
                                    'Save Data Set','Delete Case','Reload Case',...
                                    'View Case Settings'};                                               
            menu.Project(2).Callback = repmat({@obj.projectMenuOptions},[1,7]);
            
            % submenu for 'Export/Import'                                          
            menu.Project(3).List = {'Export Case','Import Case'};
            menu.Project(3).Callback = repmat({@obj.projectMenuOptions},[1,2]);
            
            %% Setup menu -------------------------------------------------
            menu.Setup(1).List = {'Import data','Site parameters',...
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
            menu.Run(2).List = {'Deepwater Waves','Nearshore Waves',...
                                'Wind-Waves','Spectral Transfer',...
                                'Wave Power','Runup',...
                                'Littoral Drift','X-shore Transport',...
                                'Overtopping','Iribarren Number'}; 
            nitems = length(menu.Run(2).List);
            menu.Run(2).Callback = repmat({@obj.runWaveModel},[1,nitems]);
            sep = repmat({'off'},[1,10]);   sep{5} = 'on';
            menu.Run(2).Separator = sep;
            
            % submenu for Beach properties
            menu.Run(3).List = {'Profiles','Shore change','Beach type',...
                                'Crenulate Bay','Shore profile','Dean profile'}; 
            runbeach = repmat({@obj.runBeachAnalysis},[1,3]);               
            menu.Run(3).Callback = [{'gcbo;'},{'gcbo;'},{@obj.runWaveModel},...
                                                            runbeach(:)'];
            menu.Run(3).Separator = {'off','off','off','off','on','off'};
            
            menu.Run(4).List = {'Volumes','Shore position',...
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
            menu.Help.List = {'Documentation','Manual'};
            menu.Help.Callback = repmat({@obj.Help},[1,2]);
            
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
            tabs.Data  = {'   Data  ',@obj.refresh};        
            tabs.Models = {'  Models  ',@obj.refresh};
            tabs.Site = {'  Site  ',''};
            subtabs.Site(1,:) = {'  Waves  ',@obj.InputTabSummary};
            subtabs.Site(2,:) = {' Simulation ',@obj.InputTabSummary};
            tabs.Plot = {'  Q-Plot  ',@obj.getTabData};
            tabs.Calcs = {'  Calcs  ',''};
            subtabs.Calcs(1,:) = {' Volumes ',@obj.setTabAction};
            subtabs.Calcs(2,:) = {' Shoreline ',@obj.setTabAction};
            subtabs.Calcs(3,:) = {'  BVI  ',@obj.setTabAction};
            subtabs.Calcs(4,:) = {' Profile ',@obj.setTabAction};
            tabs.Stats = {'   Stats   ',''};
            subtabs.Stats(1,:) = {' Descriptive ',@obj.setTabAction};
            subtabs.Stats(2,:) = {' Extremes ',@obj.setTabAction};
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
        function setTabAction(obj,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src). Called by muiModelUI.
            msg = 'No results to display';
            switch src.Tag                             
                case 'Plot' 
                    if isa(cobj,'CT_BeachAnalysis')
                        tabPlot(cobj,src,obj);
                    else
                        tabPlot(cobj,src);
                    end
                case {'Descriptive','Extremes'}
                    cobj = getClassObj(obj,'mUI','Stats',msg);
                    if isempty(cobj), return; end
                    tabStats(cobj,src);    
                case {'Volumes','Shoreline','BVI'}
                    cobj = getClassObj(obj,'Cases','CT_BeachAnalysis',msg);
                    if isempty(cobj), return; end
                    tabCalcs(cobj,src);
                case 'Profile'    
                    %warns user if no data, otherwise displays whatever is
                    %already plotted using the profile model menu options
                    getClassObj(obj,'Cases','CT_BeachAnalysis',msg);
            end
        end      
%% ------------------------------------------------------------------------
% Callback functions used by menus and tabs
%-------------------------------------------------------------------------- 
        %% File menu ------------------------------------------------------
        %use default menu functions defined in muiModelUI
            
        %% Tools menu ----------------------------------------------------- 
        function toolsMenuOptions(obj,src,~)
            %callback function for Tools menu options
            switch src.Text
                case 'Project'
                    obj.clearModel;
                case 'Figures'
                    obj.clearFigures;
                case 'Data'
                    clearCases(obj,'Data');
                case 'Models'
                    clearCases(obj,'Models');
            end
        end         
        
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
                    useCase(obj.Cases,mode,{classname},'addData');
                case 'Delete'
                    if strcmp(mode,'none'), mode = 'single'; end
                    useCase(obj.Cases,mode,{classname},'deleteData');
                case 'Quality Control'
                    if strcmp(mode,'none'), mode = 'single'; end
                    useCase(obj.Cases,mode,{classname},'qcData');
            end
            ht = findobj(obj.mUI.Tabs.Children,'Tag','Data');
            DrawMap(obj,ht);
        end        
%%
        function runProps(obj,src,~)
            %callback functions for run and settings parameter input
            tabname = [];
            switch src.Text
                case 'Wave propagation'
                    ctWaveParameters.setInput(obj);
                    tabname = 'Waves';
                case 'Wind-wave hindcast'
                    ctHindcastParameters.setInput(obj);
                    tabname = 'Waves';
                case 'Structure parameters'
                    ctStructureInput.setInput(obj);
                    tabname = 'Waves';
                case 'YGOR simulation parameters'
                    Sim_YGORinput.setInput(obj);
                    tabname = 'Simulation';
                case 'BMV simulation parameters'
                    Sim_BMVinput.setInput(obj);
                    tabname = 'Simulation';
                case 'Model constants'
                    obj.Constants = setInput(obj.Constants);
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
                    obj.mUI.ManipUI = muiManipUI.getManipUI(obj);
                case 'Simulation'
                    if ~isfield(obj.mUI,'SimUI')
                        obj.mUI.SimUI = [];
                    end
                    obj.mUI.SimUI = CT_SimUI.getCTSimUI(obj);
                case 'BVI site'  %vulnerability for single location
                    ct_sitevulnerability(obj);
                case 'BVI profile set'     %vulnerability for set of profiles
                    CT_BeachAnalysis.runModel(obj,src);
                    ht = findobj(obj.mUI.Tabs.Children,'Tag','Models');
                    DrawMap(obj,ht);
                case 'BVI set plot'     %plots for vulnerability index
                    msgtxt = 'No data available to plot';                    
                    lobj = getClassObj(obj,'Cases','CT_BeachAnalysis',msgtxt);
                    if isempty(lobj), return; end  
                    getShoreTablePlot(lobj,src.Text);
                case 'User Model'
                    %run a model in a user defined model class
                    user_model(obj)
            end
        end
%%
        function runWaveModel(obj,src,~)
            %callback functions to run wave models
            switch src.Text
                case 'Nearshore Waves'
                    %run wave model to create inshore wave data
                    ctWaveModel.runModel(obj,true); %flag is for inshore case
                case 'Deepwater Waves'
                    %run wave model to create inshore wave data
                    ctWaveModel.runModel(obj,false);%flag is for offshore case
                case 'Wind-Waves'
                    %run wave model to create inshore wave data
                    ctWindWaveModel.runModel(obj);    
                case 'Spectral Transfer'
                    msgtxt =  'WaveRayModel not found. Required for this option.';
                    isok = initialise_mui_app('WaveRayModel',msgtxt,'wrm_functions');
                    if ~isok, return; end
                    WRM_WaveModel.runModel(obj);  
                otherwise
                    CT_WaveModels.runModel(obj,src);   
            end
            ht = findobj(obj.mUI.Tabs.Children,'Tag','Models');
            DrawMap(obj,ht);
        end
%%
        function runTides(obj,src,~)
            %callback functions to run tidal analysis and reconstruction
            switch src.Text
                case 'Analysis'
                    %run tidal analysis to generate set of constituents
                    ctTidalAnalysis.runTidalAnalysis(obj);
                case 'Reconstruction'
                    %run utide to get reconstructed tidal signal
                    ctTidalAnalysis.getTidalData(obj);
            end
        end
%%
        function runBeachAnalysis(obj,src,~)
            %callback to run beach analysis functions
            %indices refer to model type in CT_BeachAnalysis
            %'Beach type','Volumes','Shoreline position','Shoreline',
            %'BVI profile set','Shore profile','Dean profile'
            msgtxt = 'No data available to plot';
            switch src.Text
                case 'Location plot'
                    lobj = getClassObj(obj,'Cases','ctBeachProfileData',msgtxt);
                    if isempty(lobj), return; end
                    getProfileLocationsPlot(lobj,obj.Cases)
                case 'Change plot'
                    lobj = getClassObj(obj,'Cases','CT_BeachAnalysis',msgtxt);
                    if isempty(lobj), return; end
                    getShorelinePlot(lobj,obj);
                case 'Rates plot'
                    lobj = getClassObj(obj,'Cases','CT_BeachAnalysis',msgtxt);
                    if isempty(lobj), return; end
                    getShoreTablePlot(lobj,src.Text);
                case 'Centroid plot'
                    lobj = getClassObj(obj,'Cases','CT_BeachAnalysis',msgtxt);
                    if isempty(lobj), return; end
                    getProfileCentroidsPlot(lobj);
                case 'Space-time plot'
                    lobj = getClassObj(obj,'Cases','CT_BeachAnalysis',msgtxt);
                    if isempty(lobj), return; end
                    getProfileSpaceTimePlot(lobj,obj);
                otherwise
                    CT_BeachAnalysis.runModel(obj,src);
                    ht = findobj(obj.mUI.Tabs.Children,'Tag','Models');
                    DrawMap(obj,ht);
            end           
        end         

        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = CT_PlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------
        function Help(~,src,~)
            %menu to access online documentation and manual pdf file
            switch src.Text
                case 'Documentation'
                    if ~isdeployed
                        doc coastaltools   %must be name of html help file
                    end
                case 'Manual'
                    if isdeployed
                        %executable using the Matlab Runtime
                        [~, result] = system('path');
                        currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
                        winopen([currentDir,filesep,'CoastalTools manual.pdf'])
                    else
                        ct_open_manual;
                    end
            end
        end

        %% Check that toolboxes are installed------------------------------
        function isok = check4muitoolbox(~)
            %check that dstoolbox and muitoolbox have been installed
            fname = 'dstable.m';
            dstbx = which(fname);
        
            fname = 'muiModelUI.m';
            muitbx = which(fname);
        
            if isempty(dstbx) && ~isempty(muitbx)
                warndlg('dstoolbox has not been installed')
                isok = false;
            elseif ~isempty(dstbx) && isempty(muitbx)
                warndlg('muitoolbox has not been installed')
                isok = false;
            elseif isempty(dstbx) && isempty(muitbx)
                warndlg('dstoolbox and muitoolbox have not been installed')
                isok = false;
            else
                isok = true;
            end
        end        
%% ------------------------------------------------------------------------
% Overload muiModelUI.MapTable to customise Data and Model Tabs
%--------------------------------------------------------------------------  
        function MapTable(obj,ht)
            %create tables for Data and Model tabs - called by DrawMap
            % load case descriptions
            muicat = obj.Cases;
            caserec = find(tabSubset(obj,ht.Tag));
            caseid = muicat.Catalogue.CaseID(caserec);
            casedesc = muicat.Catalogue.CaseDescription(caserec);
            caseclass = muicat.Catalogue.CaseClass(caserec);

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
                dst = getDataset(muicat,caserec(i),1);
                %range = getVarAttRange(dst,1,'Time');
                range = dst.RowRange;
                reclen = num2str(height(dst.DataTable));
                qualcl = '';
                if ~isempty(dst.RowRange)                    
                    stdate = char(range{1},'dd-MMM-yyyy'); %control ouput format
                    endate = char(range{2},'dd-MMM-yyyy');                    
                    if ~isempty(dst.VariableQCflags)
                        qualcl = dst.VariableQCflags{1};
                    end
                    cdata(irec,:) = {case_id,type,char(casedesc{i}),reclen,...
                                                     stdate,endate,qualcl};
                    irec = irec+1;
                else
                    cdata(irec,:) = {case_id,type,char(casedesc{i}),reclen,...
                                                             '','',qualcl};
                    irec = irec+1;
                end                
            end
            
            if strcmp(ht.Tag,'Models')
                headers = {'ID','Model Type','Model Description','Tstep',...
                'Start','End','QC'};
            else
                headers = {'ID','Data Type','Data Description','Tstep',...
                'Start','End','QC'};
            end
            
            % draw table of case descriptions
            tc=uitable('Parent',ht,'Units','normalized',...
                'CellSelectionCallback',@obj.caseCallback,...
                'Tag','cstab');
            tc.ColumnName = headers;
            tc.RowName = {};
            tc.Data = cdata;
            tc.ColumnWidth = {25 100 200 50 75 75 25};
            tc.RowStriping = 'on';
            tc.Position(3:4)=[0.935 0.8];    %may need to use tc.Extent?
            tc.Position(2)=0.9-tc.Position(4);
        end   
    end
end    
    
    
    
    
    
    
    
    
    
    