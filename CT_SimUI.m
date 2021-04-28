classdef CT_SimUI < muiDataUI
%
%-------class help---------------------------------------------------------
% NAME
%   CT_SimUI.m
% PURPOSE
%   class to implement the muiDataUI for the shoreline simulation
%   models Sim_BMVmodel and Sim_YGORmodel
% USAGE
%   obj = CT_SimUI.getCTSimUI(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   muiDataUI.m
%   Sim_BMVfitting.m, Sim_BMVmodel.m, Sim_BMVinput.m
%   Sim_YGORmodel.m, Sim_YGORinput.m
% NOTES
%   naming convention for UIs is that the first 3 letters of the claas name
%   are not included in the handle name used in proporty mobj.mUI. 
%   eg: class muiPlotsUI + mobj.mUI.Plot, and CT_SimUI + mobj.mUI.SimUI
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------    
%
    properties (Transient)
        %Abstract variables for DataGUIinterface---------------------------
        
        %names of tabs providing different data accces options
        TabOptions = {'YGOR model','BMV model'};              
        %Additional variables for application------------------------------
        Tabs2Use         %number of tabs to include  (set in getPlotGui)  
    end
%%    
    methods (Access=protected)
        function obj = CT_SimUI(mobj)
            %initialise standard figure and menus
            guititle = 'Select Data for Shoreline Simulation';
            setDataUIfigure(obj,mobj,guititle);  %initialise figure 
        end
    end    
%%
    methods (Static)
        function obj = getCTSimUI(mobj)
            %this is the function call to initialise the UI and assigning
            %to a handle of the main model UI (mobj.mUI.****) 
            obj = [];
            if isempty(mobj.Cases.Catalogue.CaseID)
                warndlg('No data available to use');                
                return;
            elseif isa(mobj.mUI.SimUI,'CT_SimUI')
                obj = mobj.mUI.SimUI;
                if isempty(obj.dataUI.Figure)
                    %initialise figure 
                    guititle = 'Select Data for Simulation';
                    setDataUIfigure(obj,mobj,guititle);   
                    setDataUItabs(obj,mobj); %add tabs 
                else
                    getdialog('Simulation UI is open');
                end
            else
                obj = CT_SimUI(mobj);
                obj.Tabs2Use = {'YGOR model','BMV model'};
                setDataUItabs(obj,mobj); %add tabs                  
            end
        end
    end  
%%  
    methods (Access=protected)
        function setTabContent(obj,src)
            %setup default layout options for individual tabs
            %Abstract function required by muiDataUI
            itab = strcmp(obj.TabOptions,src.Tag);
            obj.TabContent(itab) = muiDataUI.defaultTabContent;          
            %customise the layout of each tab. Overload the default
            %template with a function for the tab specific definition
            switch src.Tag
                case 'YGOR model'
                    if ~license('test','GADS_Toolbox')
                        warndlg('Global Optimisation Toolbox needed for YGOR model');                      
                    end
                    setYGORtab(obj,src); 
                case 'BMV model'
                    if ~license('test','GADS_Toolbox')
                        warndlg('Optimisation Toolbox needed for BMV model');
                    end
                    setBMVtab(obj,src);
            end        
        end        
      
%%
        function setVariableLists(obj,src,mobj,caserec)
            %initialise the variable lists or values
            %Abstract function required by muiDataUI
%             muicat = mobj.Cases.Catalogue;
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            sel_uic = S.Selections;
            [caserec,idx] = subSelectCases(obj,src,mobj,caserec);            
            cobj = getCase(mobj.Cases,caserec);
            for i=1:length(sel_uic)                
                switch sel_uic{i}.Tag
                    case 'Case'
                        muicat = mobj.Cases.Catalogue;
                        sel_uic{i}.String = muicat.CaseDescription(idx);
                        sel_uic{i}.UserData = sel_uic{i}.Value; %used to track changes
                    case 'Dataset'
                        sel_uic{i}.String = fieldnames(cobj.Data);
                        sel_uic{i}.Value = 1; 
                    case 'Variable'     
                        ds = fieldnames(cobj.Data);
                        sel_uic{i}.String = cobj.Data.(ds{1}).VariableDescriptions;
                        sel_uic{i}.Value = 1;
                    case 'Type'
                        sel_uic{i}.String = S.Type;   
                end
            end        
            obj.TabContent(itab).Selections = sel_uic;
        end
%%
        function useSelection(obj,src,mobj)
            %make use of the selection made to create a plot of selected type
            %Abstract function required by muiDataUI
            simoption = obj.UIsettings.Type.String;
            switch src.Parent.Tag
                case 'YGOR model'
                    switch simoption
                        case 'Fit model'
                            Sim_YGOR.runDataFitting(obj,mobj)
                        case 'Simulation'
                            Sim_YGOR.runSimulation(obj,mobj)
                        case 'Preload mean forcing'
                            Sim_YGOR.preLoad(obj,mobj)
                    end
                case 'BMV model'
%                     warndlg('Not yet available')
%                     return;
                    switch simoption
                        case 'Fit model'
                            Sim_BMV.runDataFitting(obj,mobj)
                        case 'Simulation'
                            Sim_BMV.runSimulation(obj,mobj)  
                    end                   
            end
        end 
    end
%%
%--------------------------------------------------------------------------
% Additional methods used to control selection and functionality of UI
%--------------------------------------------------------------------------            
    methods (Access=private)     
        function [caserec,idx] = subSelectCases(~,src,mobj,caserec)
            %limit the variables displayed in the UI based on tab used
            muicat = mobj.Cases;
            idwv = find(strcmp(muicat.Catalogue.CaseClass,'ctWaveModel'));
            lobj = getClassObj(mobj,'Cases','CT_BeachAnalysis');
            if strcmp(src.Tag,'YGOR model')                
                idum = find(strcmp(muicat.Catalogue.CaseClass,'muiUserModel'));                
                caseidx = getClassInstances(lobj,'ModelType','ShorePosition');
                idps = caseRec(muicat,caseidx);
                idx = sort([idwv;idum;idps]);
            else
                caseidx = getClassInstances(lobj,'ModelType','Volumes');
                idps = caseRec(muicat,caseidx);
                idx = sort([idwv,idps]);
            end
            caserec = idx(caserec);
        end
%--------------------------------------------------------------------------
% Additional methods used to define tab content
%--------------------------------------------------------------------------
        function setYGORtab(obj,src)   
            %customise the layout of the YGOR tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
           
            %Header size and text
            S.HeadPos = [0.8,0.14]; %vertical position and height of header
            txt1 = 'The YGOR model uses a representative wave energy and a';
            txt2 = 'shoreline position data set (chainage at defined elevation).';
            txt3 = 'Select data sets and assign to buttons';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            S.Titles = {'Case','Datset','Variable','Model option'};            
            S.Style = {'popupmenu','popupmenu','popupmenu','popupmenu'};
            S.Order = {'Case','Dataset','Variable','Type'};
            S.Type  = {'Fit model','Simulation','Preload mean forcing'};
            S.Scaling = {};  %options for ScaleVariable - exclude option
            
            %Tab control button options - use default values
            
            %XYZ panel definition (if required)
            S.XYZnset = 2;                        %minimum number of buttons to use
            S.XYZmxvar = [1,1];                   %maximum number of dimensions per selection
            S.XYZpanel = [0.05,0.15,0.9,0.2];     %position for XYZ button panel
            S.XYZlabels = {'Waves','Position'};   %default button labels
            
            %Action button specifications - use default values

            obj.TabContent(itab) = S;             %update object
        end    
%%      
        function setBMVtab(obj,src)
            %customise the layout of the BMV tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
           
            %Header size and text
            S.HeadPos = [0.8,0.14]; %vertical position and height of header
            txt1 = 'The BMV model uses a nearshore wave data set and a';
            txt2 = 'timeseries of beach profile volumes and centroids.';
            txt3 = 'Select data sets and assign to buttons';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            S.Titles = {'Case','Datset','Variable','Model option'};            
            S.Style = {'popupmenu','popupmenu','popupmenu','popupmenu'};
            S.Order = {'Case','Dataset','Variable','Type'};
            S.Type  = {'Fit model','Simulation','Preload mean forcing'};
            S.Scaling = {};  %options for ScaleVariable - exclude option
            
            %Tab control button options - use default values
            
            %XYZ panel definition (if required)
            S.XYZnset = 2;                        %minimum number of buttons to use
            S.XYZmxvar = [1,1];                   %maximum number of dimensions per selection
            S.XYZpanel = [0.05,0.15,0.9,0.2];     %position for XYZ button panel
            S.XYZlabels = {'Waves','Volumes'};    %default button labels
            
            %Action button specifications - use default values
            
            obj.TabContent(itab) = S;             %update object
        end
    end     
end