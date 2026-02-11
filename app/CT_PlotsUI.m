classdef CT_PlotsUI < muiDataUI
%
%-------class help------------------------------------------------------
% NAME
%   CT_PlotsUI.m
% PURPOSE
%   Class implements the muiDataUI to access data for use in
%   plotting
% SEE ALSO
%   muiDataUI.m, muiPLotsUI, muiPlots.m
%
% Author: Ian Townend
% CoastalSEA (c) Dec 2020
%--------------------------------------------------------------------------
% 
    properties (Transient)
        %Abstract variables for muiDataUI---------------------------        
        % names of tabs providing different data accces options
        TabOptions = {'Timeseries','Profiles','Rose','2D','3D','2DT','3DT'};       
        % selections that force a call to setVariableLists
        updateSelections = {'Case','Dataset'};
        %Additional variables for application------------------------------
        Tabs2Use         %number of tabs to include  (set in getPlotGui)     
    end  
%%  
    methods (Access=protected)
        function obj = CT_PlotsUI(mobj)
            %initialise standard figure and menus
            guititle = 'Plotting';
            setDataUIfigure(obj,mobj,guititle);    %initialise figure     
        end
    end
%%    
    methods (Static)
        function obj = getPlotsUI(mobj)
            %this is the function call to initialise the UI and assigning
            %to a handle of the main model UI (mobj.mUI.PlotsUI) 
            %options for selection on each tab are defined in setTabContent
            if isempty(mobj.Cases.Catalogue.CaseID)
                warndlg('No data available to plot');
                obj = [];
                return;
            elseif isa(mobj.mUI.PlotsUI,'CT_PlotsUI')
                obj = mobj.mUI.PlotsUI;
                if isempty(obj.dataUI.Figure)
                    %initialise figure 
                    guititle = 'Select Data for Plotting';
                    setDataUIfigure(obj,mobj,guititle);    
                    setDataUItabs(obj,mobj); %add tabs 
                else
                    getdialog('Plot UI is open');
                end
            else
                obj = CT_PlotsUI(mobj);
                idx = find(~ismember(mobj.DataUItabs.Plot,obj.TabOptions));
                if ~isempty(idx)
                    txt1 = 'Unknown plot type defined in main UI for DataUItabs.Plot';
                    for j=1:length(idx)
                        txt1 = sprintf('%s\n"%s" not defined',txt1,mobj.DataUItabs.Plot{idx(j)});
                    end
                    warndlg(txt1)
                    obj = [];
                    return
                end
                obj.Tabs2Use = mobj.DataUItabs.Plot;
                setDataUItabs(obj,mobj); %add tabs                
            end                
        end
    end
%%
%--------------------------------------------------------------------------
% Abstract methods required by muiDataUI to define tab content
%--------------------------------------------------------------------------
    methods (Access=protected) 
        function setTabContent(obj,src)
            %setup default layout options for individual tabs
            %Abstract function required by muiDataUI
            itab = find(strcmp(obj.Tabs2Use,src.Tag));
            obj.TabContent(itab) = muiDataUI.defaultTabContent;
            
            %customise the layout of each tab. Overload the default
            %template with a function for the tab specific definition
            switch src.Tag
                case 'Timeseries'
                    setTS_tab(obj,src);
                case 'Profiles'
                    setProf_tab(obj,src);
                case 'Rose'
                    setRose_tab(obj,src);    
                case '2D'
                    set2D_tab(obj,src);
                case '3D'
                    set3D_tab(obj,src);  
                case '2DT'
                    set2DT_tab(obj,src);
                case '3DT'
                    set3DT_tab(obj,src);     
            end             
        end                
%%
        function setVariableLists(obj,src,mobj)
            %initialise the variable lists or values
            %Abstract function required by muiDataUI
            %called during initialisation by muiDataUI.setDataUItabs and 
            %whenever muiDataUI.updateCaseList is called
            muicat = mobj.Cases.Catalogue;
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            sel_uic = S.Selections;
            listrec = sel_uic{strcmp(S.Order,'Case')}.Value; 
            caserec = subCaseRec(obj,src,muicat,listrec);
            if isempty(caserec), return; end
            cobj = getCase(mobj.Cases,caserec);
            dsnames = fieldnames(cobj.Data);
            ids = sel_uic{strcmp(S.Order,'Dataset')}.Value;
            for i=1:length(sel_uic)                
                switch sel_uic{i}.Tag
                    case 'Case'                        
                        sel_uic{i}.String = subCaseList(obj,src,muicat);
                        sel_uic{i}.UserData = sel_uic{i}.Value; %used to track changes
                    case 'Dataset'
                        sel_uic{i}.String = dsnames;                        
                        if sel_uic{i}.Value>length(dsnames)
                            sel_uic{i}.Value = 1; %reset if case changed to
                            ids = 1;              %case with fewer datasets
                        end
                        sel_uic{i}.UserData = sel_uic{i}.Value; %used to track changes
                    case 'Variable' 
                        sel_uic{i}.String =  subVarList(obj,src,cobj,ids);
                        sel_uic{i}.Value = 1;
                    case 'Type'
                        sel_uic{i}.String = S.Type;
                    case 'Other'
                        sel_uic{i}.String = 'Not set';
                end
            end        
            obj.TabContent(itab).Selections = sel_uic;
        end
%%
        function caserec = subCaseRec(obj,src,muicat,listrec)
            %find correct caserec from listrec (id in sublist)
            cdesc = subCaseList(obj,src,muicat);
            if ~isempty(cdesc)
                casedesc = cdesc{listrec};
                caserec = find(strcmp(muicat.CaseDescription,casedesc)); 
            else
                caserec = [];
            end
        end
%%
        function cdesc = subCaseList(~,src,muicat)
            %modify list to a subset of the full case list
            classes = muicat.CaseClass;
            switch src.Tag
                case 'Timeseries'
                    idx = ~ismember(classes,'ctBeachProfileData');                    
                case 'Profiles'
                    idx = ismember(classes,'ctBeachProfileData');
                case 'Rose'
                    roseclasses = {'ctWaveData','ctWindData','ctWaveModel',...
                                   'ctWindWaveModel','WRM_WaveModel'...
                                   'muiTableImport','muiUserData','muiUserModel'};
                    idx = ismember(classes,roseclasses);
                otherwise
                    idx = true(length(classes),1);
            end
            cdesc = muicat.CaseDescription(idx);
        end
%%
        function varnum = subVarRec(obj,src,cobj,listnum,ids)
            %find correct varnum from subvarlist (id in sub-list)
            vdesc = subVarList(obj,src,cobj,ids);
            vardesc = vdesc{listnum};
            varnum = find(strcmp(varlist,vardesc));
        end
%%
        function [vdesc,varlist] = subVarList(~,src,cobj,ids)
            %modify list to a subset of the full case list
            ds = fieldnames(cobj.Data);
            varlist = cobj.Data.(ds{ids}).VariableDescriptions;            
             switch src.Tag                  
                case 'Profiles'
                    pvars = {'FeatureCode','Elevation'};
                    idx = ismember(varlist,pvars);
                    vdesc = varlist(idx);
                otherwise
                    vdesc = varlist;
             end
        end
%%        
        function useSelection(obj,src,mobj)  
            %make use of the selection made to create a plot of selected type
            %Abstract function required by muiDataUI
            if strcmp(src.String,'Save')   %save animation to file
                saveAnimation(obj,src,mobj);
            else
                %some tabs use sub-selections of Cases and/or Variables
                %however muiDataUI already checks for subselection and adjusts
                CT_Plots.getPlot(obj,mobj);
            end
        end     
    end
%%
%--------------------------------------------------------------------------
% Additional methods used to define tab content
%--------------------------------------------------------------------------
    methods (Access=private)
        function setTS_tab(obj,src)
            %customise the layout of the timeseries tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a timeseries plot (line, bar, etc) select Case, Dataset and Variable';
            txt2 = 'Assign to the TSvar button and adjust the variable and time ranges, if required';
            txt3 = 'Select the Type of plot and then use the New, Add or Delete buttons to continue';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);  
            %Specification of uicontrol for each selection variable  
            %Use default lists
            
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.1,0.14;0.3,0.14;0.5,0.14;0.7,0.14]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 1;                           %minimum number of buttons to use
            S.XYZmxvar = 1;                          %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.25,0.91,0.15];      %position for XYZ button panel
            S.XYZlabels = {'TSvar'};                 %default button labels
            
            %Action button specifications - use default
            %Action button specifications
            S.ActButNames = {'Refresh','RightYaxis'}; %names assigned selection struct
            S.ActButText = {char(174),'yL'};     %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.30];  %positions for action buttons   
            % action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setRightYaxis(src,evt)'};
            % tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                           'Swap from left to right axis'};         
            obj.TabContent(itab) = S;            %update object
        end
%%
        function setProf_tab(obj,src)
            %customise the layout of the profiles tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For XY profile plots, select  Case, Dataset and Variable and assign to the ''Profile 1'' button';
            txt2 = 'For XT-Z profile plots, assign profile selection to the ''Profile 1'' button and Time to the ''T or PrN'' button';
            txt3 = 'For multi-profile plots, assign alongshore Start and End profiles to the ''Profile 1'' and ''T or PrN'' buttons';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3); 
            
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.1,0.14;0.3,0.14;0.5,0.14;0.7,0.14]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 1;                           %minimum number of buttons to use
            S.XYZmxvar = [2,2];                      %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.25,0.91,0.2];       %position for XYZ button panel
            S.XYZlabels = {'Profile 1','T or PrN'};   %default button labels
            %Action button specifications - use default
            
            obj.TabContent(itab) = S;                %update object
        end
%%
        function setRose_tab(obj,src)
            %customise the layout of the rose tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a Polar or Rose plot select the Case, Dataset and Variable';
            txt2 = 'Assign to the ''Rose'' buttton';
            txt3 = 'Add and Delete buttons only work for Polar plots';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3); 
            
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.1,0.14;0.3,0.14;0.5,0.14;0.7,0.14]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 1;                           %minimum number of buttons to use
            S.XYZmxvar = 1;                          %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.25,0.91,0.15];      %position for XYZ button panel
            S.XYZlabels = {'Rose'};               %default button labels
            
            %Action button specifications - use default
            
            obj.TabContent(itab) = S;                %update object            
        end
%%
        function set2D_tab(obj,src)
            %customise the layout of the 2D tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a cartesian plot (line, bar, etc) select variable and dimension to use for X and Y axes';
            txt2 = 'For a Polar or Rose plot, switch the + button to O and set X to a direction variable';
            txt3 = 'Use the XY button to swap the X and Y variables without needing to reassign them';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            %Use default lists
            
            %Tab control button options
            S.TabButText = {'New','Add','Delete','Clear'}; %labels for tab button definition
            S.TabButPos = [0.1,0.14;0.3,0.14;0.5,0.14;0.7,0.14]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 2;                           %minimum number of buttons to use
            S.XYZmxvar = [1,1];                      %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.25,0.91,0.2];       %position for XYZ button panel
            S.XYZlabels = {'Var','X'};               %default button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','RightYaxis','Swap','Polar'}; %names assigned selection struct
            S.ActButText = {char(174),'yL','XY','+'};     %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.38;0.895,0.32;0.895,0.26];%positions for action buttons   
            % action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setRightYaxis(src,evt)',...
                            '@(src,evt)setXYorder(src,evt)',...
                            '@(src,evt)setPolar(src,evt)'};
            % tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                           'Swap from left to right axis',...
                           'Swap from X-Y to Y-X axes',...
                           'XY to Polar; X data in degrees or radians'};           
            obj.TabContent(itab) = S;                %update object
        end
%%
        function set3D_tab(obj,src)
            %customise the layout of the 3D tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            
            %Header size and text
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For a contour or surface plot, select a variable with at least 2 dimensions';
            txt2 = 'If the variable has more than 2 dimensions you will be prompted to select a sub-set';
            txt3 = 'Select dimensions to use for the X-Y axes';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            %Use default lists except for:
            S.Type = {'surf','contour','contourf','contour3','mesh','User'}; 
            
            %Tab control button options - use defaults
           
            %XYZ panel definition (if required)
            S.XYZnset = 3;                         %minimum number of buttons to use
            S.XYZmxvar = [2,1,1];                  %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.14,0.91,0.3];     %position for XYZ button panel
            S.XYZlabels = {'Var','X','Y'};         %button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','Polar'};   %names assigned selection struct
            S.ActButText = {char(174),'+'};        %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.27];    %positions for action buttons   
            %action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setPolar(src,evt)'};
            %tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                            'XY to Polar; X data in degrees'}; 
                        
            obj.TabContent(itab) = S;              %update object
        end
%%
        function set2DT_tab(obj,src)
            %customise the layout of the 2DT tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            
            %Header size and text
            S.HeadPos = [0.86,0.1];    %header vertical position and height
            txt1 = 'For an animation of a line, surface or volume';
            txt2 = 'Select the Variable, Time and the property to define the X-axis';
            txt3 = 'Ensure that Time and X selections are correctly matched to the dimensions of the variable';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            
            %Specification of uicontrol for each selection variable  
            %Use default lists 
            
            %Tab control button options
            S.TabButText = {'Run','Save','Clear'};  %labels for tab button definition
            S.TabButPos = [0.1,0.03;0.3,0.03;0.5,0.03]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 3;                          %minimum number of buttons to use
            S.XYZmxvar = [2,1,1];                   %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.14,0.91,0.30];     %position for XYZ button panel
            S.XYZlabels = {'Var','Time','X'}; %button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','Polar'};   %names assigned selection struct
            S.ActButText = {char(174),'+'};        %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.17];    %positions for action buttons   
            %action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setPolar(src,evt)'};
            %tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                            'XY to Polar; X data in degrees'};   
     
            obj.TabContent(itab) = S;               %update object
        end   
%%
        function set3DT_tab(obj,src)
            %customise the layout of the 2DT tab
            %overload defaults defined in muiDataUI.defaultTabContent
            itab = strcmp(obj.Tabs2Use,src.Tag);
            S = obj.TabContent(itab);
            
            %Header size and text
            S.HeadPos = [0.88,0.1];    %header vertical position and height
            txt1 = 'For an animation of a line, surface or volume, select a variable with at least 2 dimensions';
            txt2 = 'Select the Time and properties to define the X and Y axes';
            txt3 = 'Ensure that Time and XY selections are correctly matched to the dimensions of the variable';
            S.HeadText = sprintf('1 %s\n2 %s\n3 %s',txt1,txt2,txt3);
            %Specification of uicontrol for each selection variable  
            %Use default lists except for:
            S.Type = {'surf','contour','contourf','contour3','mesh','User'}; 
            
            %Tab control button options
            S.TabButText = {'Run','Save','Clear'};  %labels for tab button definition
            S.TabButPos = [0.1,0.03;0.3,0.03;0.5,0.03]; %default positions
            
            %XYZ panel definition (if required)
            S.XYZnset = 3;                          %minimum number of buttons to use
            S.XYZmxvar = [3,1,1,1];                 %maximum number of dimensions per selection
            S.XYZpanel = [0.04,0.14,0.91,0.40];     %position for XYZ button panel
            S.XYZlabels = {'Var','Time','X','Y'}; %button labels
            
            %Action button specifications
            S.ActButNames = {'Refresh','Polar'};   %names assigned selection struct
            S.ActButText = {char(174),'+'};        %labels for additional action buttons
            % Negative values in ActButPos indicate that a
            % button is alligned with a selection option numbered in the 
            % order given by S.Titles
            S.ActButPos = [0.86,-1;0.895,0.27];    %positions for action buttons   
            %action button callback function names
            S.ActButCall = {'@(src,evt)updateCaseList(obj,src,evt,mobj)',...
                            '@(src,evt)setPolar(src,evt)'};
            %tool tips for buttons             
            S.ActButTip = {'Refresh data list',...
                            'XY to Polar; X data in degrees'}; 
    
            obj.TabContent(itab) = S;               %update object            
        end       
    end
end